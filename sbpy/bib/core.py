# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Bibliography Tracking Module

sbpy classes and functions can automatically report citations to the
user through a bibliography registry.
"""

__all__ = ['register', 'reset', 'status', 'stop', 'track', 'Tracking',
           'to_text', 'to_bibtex']

from collections import OrderedDict


def register(task, citations):
    """Register a citation with the `sbpy` bibliography when running a
    function

    Parameters
    ----------
    task : string
      The name of the source module requesting a citation.
    citations : dict
      A dictionary of NASA Astrophysics Data System (ADS) bibcode(s)
      to cite.  The keys reference the aspect that requires citation,
      e.g., `{'method': '1998Icar..131..291H'}`, or
      `{'beaming parameter': '2013Icar..226.1138F'}`.

    """
    global _bibliography, _track
    if _track:
        if task in _bibliography:
            key = []
            value = set()
            key.append(list(citations.keys()))
            for i in range(len(key)):
                if list(citations.keys())[i] in _bibliography[task]: #if subtask already exists
                    subtasks = _bibliography[task]
                    if isinstance(list(subtasks.values())[i], list):
                        for j in range(len(list(subtasks.values())[i])):
                            value.add(list(subtasks.values())[i][j])
                    else:
                        value.add(list(subtasks.values())[i])
                        value.add(list(citations.values())[i])
                        subtasks[list(citations.keys())[i]] = value #append new citation to existing subtask
                else: #if subtask doesn't exists
                    subtasks = _bibliography[task]
                    value.add(list(citations.values())[i])
                    subtasks[list(citations.keys())[i]] = value
        else:
            _bibliography[task] = citations


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
            for key, val in ref.items():
                if isinstance(val, str):
                    val = [val]
                val = list(val)
                for i in range(len(val)):
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')
                        try:
                            # request needed fields to avoid lazy loading
                            paper = list(
                                ads.SearchQuery(
                                    bibcode=val[i],
                                    fl=['first_author', 'author', 'year']
                                ))[0]
                        except Warning as e:
                            # if query failed,
                            output += '  {:s}: {:s}\n'.format(key, val[i])
                            continue

                    if len(paper.author) > 4:
                        author = '{:s} et al.'.format(
                            paper.first_author.split(',')[0])
                    elif len(paper.author) > 1:
                        author = ','.join([au.split(',')[0] for au in
                                           paper.author[:-1]])
                        author += ' & {:s}'.format(paper.author[-1].split(',')[0])
                    else:
                        # single author
                        author = paper.first_author.split(',')[0]

                    output += '  {:s}: {:s} {:s}, {:s}\n'.format(key, author,
                                                                 paper.year, val[i])
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
