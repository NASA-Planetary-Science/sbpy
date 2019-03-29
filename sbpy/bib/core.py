# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Bibliography Tracking Module

sbpy classes and functions can automatically report citations to the
user through a bibliography registry.
"""

__all__ = ['register', 'reset', 'status', 'stop', 'track', 'Tracking',
           'to_text', 'to_bibtex', 'to_aastex', 'to_icarus', 'to_mnras']

import warnings
from collections import OrderedDict
import atexit
from astropy import log


def register(task, citation):
    """Register a citation with the `sbpy` bibliography tracker.


    Parameters
    ----------
    task : string
        The name of the source module requesting a citation.

    citation : dict
        A dictionary of a single NASA Astrophysics Data System (ADS)
        bibcode to cite.  The key references the aspect that requires
        citation, e.g., `{'method': '1998Icar..131..291H'}`, or
        `{'beaming parameter': '2013Icar..226.1138F'}`.

    """
    global _bibliography, _track
    if _track:
        if task in _bibliography:
            for newsubtask, newcitation in citation.items():
                if newsubtask in _bibliography[task].keys():
                    _bibliography[task][newsubtask].update([newcitation])
                else:
                    _bibliography[task][newsubtask] = set(
                        [list(citation.values())[0]])
        else:
            thiscitation = OrderedDict()
            for newsubtask, newcitation in citation.items():
                thiscitation[newsubtask] = set([list(citation.values())[0]])
            _bibliography[task] = thiscitation


def reset():
    """Reset `sbpy` bibliography tracking."""
    global _bibliography
    _bibliography = OrderedDict()


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


class Tracking:

    def __enter__(self):
        track()

    def __exit__(self, type, value, tb):
        stop()


def to_text(filter=None):
    """convert bibcodes to human readable text

    Returns
    -------
    text : string
        Text with human-readable information retrieved from ADS using
        the bibcodes.

    """
    import ads

    output = ''
    for task, ref in _bibliography.items():
        output += '{:s}:\n'.format(task)
        try:
            if filter is None:
                pass
            else:
                ref = {i: ref[i] for i in ref if i == filter}
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
                        # more than 3 authors
                        author = '{:s} et al.'.format(
                            ', '.join([au.split(',')[0] for au in
                                       paper.author[:3]]))
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


def _to_format(format):
    """Convert bibcodes to a range of different output formats.


    Parameters
    ----------
    format : string
        Output format: ``bibtex`` | ``aastex``  | ``icarus`` | ``mnras``


    Returns
    -------
    text : string
        ADS entries for all the bibliographic items in the given
        format.  Uses a query to the export service to get the data
        for each reference.

    """
    import ads

    output = ''
    for task, ref in _bibliography.items():
        with warnings.catch_warnings():
            # warnings.filterwarnings('error')
            try:
                for key, val in ref.items():
                    # This method avoids using multiple calls to the
                    # API that may impact rate limits
                    # https://github.com/adsabs/adsabs-dev-api/blob/master/
                    # export.md
                    query = ads.ExportQuery(list(val), format=format)
                    data = query.execute()
                    output += '% {:s}/{:s}:\n{:s}\n'.format(task, key,
                                                            data)
            except ads.exceptions.APIResponseError as e:
                e = str(e)
                if '<title>' in e:
                    e = e[e.find('<title>')+7: e.find('</title>')]
                warnings.warn('cannot obtain ADS data for {:s}/{:s}: ({:s})'.
                              format(task, key, e),
                              RuntimeWarning)
                pass

    return output


def to_bibtex():
    """Convert bibliography to BibTeX format.


    Returns
    -------
    text : string
        ADS data in BibTeX format.

    """
    return _to_format('bibtex')


def to_aastex():
    """Convert bibliography to AASTeX format.


    Returns
    -------
    text : string
        ADS data in AASTeX format.

    """
    return _to_format('aastex')


def to_icarus():
    """Convert bibliography to Icarus LATeX format.


    Returns
    -------
    text : string
        ADS data in Icarus LATeX format.

    """
    return _to_format('icarus')


def to_mnras():
    """Convert bibliography to MNRAS LATeX format.


    Returns
    -------
    text : string
        ADS data in MNRAS LATeX format.

    """
    return _to_format('mnras')


@atexit.register
def _report_at_exit():
    if _track:
        log.info('Thank you for using sbpy.  ' +
                 'Your session results were based on:\n' +
                 to_text())


_track = False  # default is no bibliography tracking
_bibliography = OrderedDict()
