# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Bibliography Tracking Module

sbpy classes and functions can automatically report citations to the
user through a bibliography registry.
"""

__all__ = [
    'register',
    'reset',
    'status',
    'stop',
    'track',
    'Tracking',
    'cite',
    'show',
    'to_text',
    'to_bibtex',
    'to_aastex',
    'to_icarus',
    'to_mnras'
]

__doctest_requires__ = {
    "cite": ["ads"],
}

import warnings
from functools import wraps
from collections import OrderedDict, defaultdict
import atexit

try:
    import ads
except ImportError:
    ads = None

from astropy import log

from ..utils.decorators import requires


def register(task, citations):
    """Register a citation with the `sbpy` bibliography tracker.


    Parameters
    ----------
    task : string or function
        The source requesting a citation.

    citations : dict
        A dictionary of NASA Astrophysics Data System (ADS) bibcodes,
        DOIs, or free-form strings to cite.  The key is the aspect
        that requires citation.


    Examples
    --------
    >>> register('sbpy.thermal.neatm', {
    ...     'method': '1998Icar..131..291H',
    ...     'parameter: beaming': '2013Icar..226.1138F'
    ... })

    Citations may also be a list.

    >>> register('user task', {
    ...     'classification': ['1950BAN....11...91O',
    ...                        '1978M&P....19..305W']
    ... })

    """

    global _bibliography, _track

    if isinstance(task, str):
        source = task
    else:
        source = '.'.join((task.__module__, task.__qualname__))

    if _track:
        for key, citation in citations.items():
            c = [citation] if isinstance(citation, str) else list(citation)
            _bibliography[source][key].update(c)


def _bibliography_task_generator():
    """Generator for empty bibliography tasks."""
    return defaultdict(set)


def reset():
    """Reset `sbpy` bibliography tracking."""
    global _bibliography
    _bibliography = defaultdict(_bibliography_task_generator)


def stop():
    """Disable `sbpy` bibliography tracking."""
    global _track
    _track = False


def status():
    """Report `sbpy` bibliography tracking status.


    Returns
    -------
    status : bool
        `True` if bibliography tracking is enabled.

    """
    return _track


def track():
    """Enable `sbpy` bibliography tracking."""
    global _track
    _track = True
    register('sbpy', {'software: sbpy': [
             '2019JOSS....4.1426M']})


class Tracking:
    def __init__(self, reporter=None):
        self.reporter = reporter

    def __enter__(self):
        track()

    def __exit__(self, type, value, tb):
        stop()
        if self.reporter is not None:
            print(self.reporter())


def cite(citations):
    """Decorator that registers citations within ``sbpy``.

    This decorator is the primary mechanism for registering citations.
    As an alternative, `~register` may be called directly.


    Parameters
    ----------
    citations : dict
        A dictionary of NASA Astrophysics Data System (ADS) bibcodes,
        DOIs, or free-form strings to cite.  The key references the
        aspect that requires citation, e.g., `{'method':
        '1998Icar..131..291H'}`, or `{'parameter : beaming':
        '2013Icar..226.1138F'}`.


    Returns
    -------
    decorator : function
        Function decorator.


    Examples
    --------
    >>> @cite({'method': '1687pnpm.book.....N'})
    ... def force(mass, acceleration):
    ...     return mass * acceleration
    >>>
    >>> with Tracking(to_text):
    ...     print(force(1, 2))    # doctest: +REMOTE_DATA
    2
    sbpy:
      software: sbpy:
          Mommert, Kelley, de Val-Borro, Li et al. 2019, The Journal of Open Source Software, Vol 4, 38, 1426
    sbpy.bib.core.force:
      method:
          Newton 1687, Philosophiae Naturalis Principia Mathematica. Auctore Js. Newton
    """

    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            # only cite after successful call
            result = f(*args, **kwargs)
            register(f, citations)
            return result
        return wrapper
    return decorator


def _filter(filter):
    """Private function that filters by key."""
    if filter is None:
        return _bibliography

    filtered = defaultdict(OrderedDict)
    for task, ref in _bibliography.items():
        for key in ref:
            if key.startswith(filter):
                filtered[task][key] = ref[key]

    return filtered


def show(filter=None):
    """Show current bibliography.

    Bibcodes are not coverted to text.  Use this function when you do
    not have the NASA `~ads` module or internet access.


    Parameters
    ----------
    filter : string, optional
        Filter the bibliography by key, showing only those that start
        with this string.


    Returns
    -------
    bibliography : string

    """

    output = ''
    for task, ref in _filter(filter).items():
        output += '{:s}:\n'.format(task)
        for key, citations in ref.items():
            output += '  {:s}:\n'.format(key)
            for citation in citations:
                output += '    {:s}\n'.format(citation)

    return output


@requires("ads")
def to_text(filter=None):
    """Convert bibcodes to human readable text.


    Parameters
    ----------
    filter : string, optional
        Filter the bibliography by key, showing only those that start
        with this string.


    Returns
    -------
    text : string
        Text with human-readable information retrieved from ADS using
        the bibcodes.

    """

    output = ''
    for task, ref in _filter(filter).items():
        output += '{:s}:\n'.format(task)
        try:
            for key, value in ref.items():
                output += '  {:s}:\n'.format(key)
                for citation in value:
                    with warnings.catch_warnings():
                        # warnings.filterwarnings('error')
                        try:
                            # request needed fields to avoid lazy loading
                            paper = list(
                                ads.SearchQuery(
                                    bibcode=citation,
                                    fl=['first_author', 'author', 'volume',
                                        'pub', 'issue', 'page', 'year']
                                ))[0]
                        except (IndexError, Warning, RuntimeWarning) as e:
                            # if query failed,
                            output += '      {:s}\n'.format(citation)
                            continue

                    # format authors
                    if len(paper.author) > 4:
                        # more than 4 authors
                        author = '{:s} et al.'.format(
                            ', '.join([au.split(',')[0] for au in
                                       paper.author[:4]]))
                    elif len(paper.author) > 1:
                        # less than or equal to 3 authors
                        author = ', '.join([au.split(',')[0] for au in
                                            paper.author[:-1]])
                        author += ' & {:s}'.format(paper.author[-1].
                                                   split(',')[0])
                    else:
                        # single author
                        author = paper.first_author.split(',')[0]

                    # year, journal
                    output += '      {:s} {:s}, {:s}'.format(
                        author, paper.year, str(paper.pub))

                    # volume
                    if paper.volume is not None:
                        output += ', Vol {:s}'.format(str(paper.volume))

                    # issue
                    if paper.issue is not None:
                        output += ', {:s}'.format(str(paper.issue))

                    # page
                    if paper.page is not None:
                        if len(paper.page) == 2:
                            output += ', {:s}-{:s}'.format(
                                str(paper.page[0]), str(paper.page[1]))
                        else:
                            output += ', {:s}'.format(str(paper.page[0]))

                    output += '\n'

        except AttributeError:
            pass

    return output


@requires("ads")
def _to_format(format, filter=None):
    """Convert bibcodes to a range of different output formats.


    Parameters
    ----------
    format : string
        Output format: ``bibtex`` | ``aastex``  | ``icarus`` | ``mnras``.

    filter : string, optional
        Filter the bibliography by key, showing only those that start
        with this string.


    Returns
    -------
    text : string
        ADS entries for all the bibliographic items in the given
        format.  Uses a query to the export service to get the data
        for each reference.

    """

    output = ''
    for task, ref in _filter(filter).items():
        with warnings.catch_warnings():
            # warnings.filterwarnings('error')
            for key, val in ref.items():
                # This method avoids using multiple calls to the
                # API that may impact rate limits
                # https://github.com/adsabs/adsabs-dev-api/blob/master/Export_API.ipynb
                try:
                    query = ads.ExportQuery(list(val), format=format)
                    data = query.execute()
                    output += '% {:s}/{:s}:\n{:s}\n'.format(
                        task, key.replace(' ', '_'), data)
                except ads.exceptions.APIResponseError:
                    output += '% {:s}/{:s}:\n{:s}\n\n'.format(
                        task, key.replace(' ', '_'), ", ".join(list(val)))

    return output


def to_bibtex(filter=None):
    """Convert bibliography to BibTeX format.


    Parameters
    ----------
    filter : string, optional
        Filter the bibliography by key, showing only those that start
        with this string.


    Returns
    -------
    text : string
        ADS data in BibTeX format.

    """
    return _to_format('bibtex', filter=filter)


def to_aastex(filter=None):
    """Convert bibliography to AASTeX format.


    Parameters
    ----------
    filter : string, optional
        Filter the bibliography by key, showing only those that start
        with this string.


    Returns
    -------
    text : string
        ADS data in AASTeX format.

    """
    return _to_format('aastex', filter=filter)


def to_icarus(filter=None):
    """Convert bibliography to Icarus LATeX format.


    Parameters
    ----------
    filter : string, optional
        Filter the bibliography by key, showing only those that start
        with this string.


    Returns
    -------
    text : string
        ADS data in Icarus LATeX format.

    """
    return _to_format('icarus', filter=filter)


def to_mnras(filter=None):
    """Convert bibliography to MNRAS LATeX format.


    Parameters
    ----------
    filter : string, optional
        Filter the bibliography by key, showing only those that start
        with this string.


    Returns
    -------
    text : string
        ADS data in MNRAS LATeX format.

    """
    return _to_format('mnras', filter=filter)


@atexit.register
def _report_at_exit():
    if status():
        log.info('Thank you for using sbpy.  ' +
                 'Your session results were based on:\n' +
                 show())


stop()  # default is no bibliography tracking
reset()  # creates empty _bibliography
