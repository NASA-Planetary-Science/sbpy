# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Bibliography Tracking Module

sbpy classes and functions can automatically report citations to the
user through a bibliography registry.
"""

__all__ = ['register', 'reset', 'status', 'stop', 'track', 'Tracking',
           'to_text', 'to_bibtex']

from collections import OrderedDict


def register(task, citation):
    """Register a citation with the `sbpy` bibliography when running a
    function

    Parameters
    ----------
    task : string
      The name of the source module requesting a citation.
    citation : dict
      A dictionary of a single NASA Astrophysics Data System (ADS) bibcode
      to cite.  The key references the aspect that requires citation,
      e.g., `{'method': '1998Icar..131..291H'}`, or
      `{'beaming parameter': '2013Icar..226.1138F'}`.

    """
    global _bibliography, _track
    if _track:
        if task in _bibliography:
            for newsubtask, newcitation in citation.items():
                if newsubtask in _bibliography[task].keys():
                    _bibliography[task][newsubtask].update([newcitation])
                else:
                    _bibliography[task][newsubtask] = set([list(citation.values())[0]])
        else:
            for newsubtask, newcitation in citation.items():
                citation[newsubtask] = set([list(citation.values())[0]])
            _bibliography[task] = citation


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


def to_text():
    """convert bibcodes to human readable text

    Returns
    -------
    text : string
      Text with human-readable information retrieved from ADS using the
      bibcodes.
    """
    import ads
    import warnings

    output = ''
    for task, ref in _bibliography.items():
        output += '{:s}:\n'.format(task)
        try:
            for key, value in ref.items():
                output += '  {:s}:\n'.format(key)
                for citation in value:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')
                        try:
                            # request needed fields to avoid lazy loading
                            paper = list(
                                ads.SearchQuery(
                                    bibcode=citation,
                                    fl=['first_author', 'author', 'year']
                                ))[0]
                        except (IndexError, Warning, RuntimeWarning) as e:
                            # if query failed,
                            output += '      {:s}\n'.format(citation)
                            continue

                    if len(paper.author) > 4:
                        author = '{:s} et al.'.format(
                            paper.first_author.split(',')[0])
                    elif len(paper.author) > 1:
                        author = ','.join([au.split(',')[0] for au in
                                           paper.author[:-1]])
                        author += ' & {:s}'.format(paper.author[-1].split(',')
                                                   [0])
                    else:
                        # single author
                        author = paper.first_author.split(',')[0]

                    output += '      {:s} {:s}, {:s}\n'.format(author,paper.year,
                                                         citation)
        except AttributeError:
            pass

    return output


def to_bibtex():
    """Convert bibcodes to LaTeX BibTeX

    Returns
    -------
    text : string
      ADS BibTex entries for all the bibliographic items.  Uses a
      query to the export service to get the BibTeX for each
      reference.

    """
    import ads

    output = ''
    for task, ref in _bibliography.items():
        try:
            for key, val in ref.items():
                # This method to get bibtex records avoids using
                # multiple calls to the API that may impact rate
                # limits
                # https://github.com/adsabs/adsabs-dev-api/blob/master/export.md
                query = ads.ExportQuery(val, format='bibtex')
                bibtex = query.execute()
                output += '% {:s}/{:s}:\n{:s}\n'.format(task, key, bibtex)
        except ads.exceptions.APIResponseError:
            pass

    return output


_track = False  # default is no bibliography tracking
_bibliography = OrderedDict()
