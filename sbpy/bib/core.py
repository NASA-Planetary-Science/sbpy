# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=================================
SBPy Bibliography Tracking Module
=================================

`sbpy` classes and functions can automatically report citations to the
user through a bibliography registry.  Use `track` to enable citation
tracking, and `citations` to report the references used.

Example
-------
>>> from sbpy import bib, data
>>> bib.track()
>>> eph = data.Ephem.from_horizons('433', epoch=None, observatory='500')
>>> print(bib.to_text())  # doctest: +REMOTE_DATA
sbpy.data.Ephem:
  implementation: Giorgini et al. 1996, 1996DPS....28.2504G

Bibliography tracking can be used in a context manager::

  >>> from sbpy import bib
  >>> from sbpy.data import Ephem
  >>> with bib.Tracking():
  >>>     eph = Ephem.from_horizons('encke', epoch=None, observatory='500')
  >>> bib.to_text()
  JPL Horizons:
    implementation: 1996DPS....28.2504G

Functions
---------
register : Register a citation.
status   : Bibliography tracking status.
stop     : Stop bibliography tracking.
track    : Start bibliography tracking.
reset    : clear bibliography.
to_text  : output bibliography in clear text
to_bibtex: output bibliography in bibtex format

Context managers
----------------
Tracking : Bibliography tracking context manager.

"""

__all__ = ['register', 'reset', 'status', 'stop', 'track', 'Tracking',
           'to_text', 'to_bibtex']

from collections import OrderedDict


def register(task, citations):
    """Register a citation with the `sbpy` bibliography.

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
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        # request needed fields to avoid lazy loading
                        paper = list(
                            ads.SearchQuery(
                                bibcode=val,
                                fl=['first_author', 'author', 'year']
                            ))[0]
                    except Warning as e:
                        # if query failed,
                        output += '  {:s}: {:s}\n'.format(key, val)
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
                                                             paper.year, val)
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
