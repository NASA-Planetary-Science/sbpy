# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
================================
SBPy Biliography Tracking Module
================================

`sbpy` classes and functions can automatically report citations to the
user through a biliography registry.  Use `track` to enable citation
tracking, and `citations` to report the references used.

Example
-------
>>> from sbpy import bib, data
>>> bib.track()
>>> eph = data.Ephem.from_horizons('encke')
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

"""

__all__ = ['register', 'reset', 'status', 'stop', 'track',
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


def to_text():
    """convert bibcodes to human readable text
    
    not yet implemented
    """
    output = ''
    for task, ref in _bibliography.items():
        output += '{:s}:\n'.format(task)
        try:
            for key, val in ref.items():
                # transform val to ads(val)
                output += '  {:s}: {:s}\n'.format(key, val)
        except:
            pass
        
    return output

def to_bibtex():
    """ convert bibcodes to LATeX bibtex
    
        not yet implemented
    """
    output = ''
    for task, ref in _bibliography.items():
        try:
            for key, val in ref.items():
                # transform val to ads.bibtex(val)
                output += '% {:s}/{:s}:\n{:s}\n'.format(task, key, val)
        except:
            pass

    return output

    
_track = False  # default is no bibliography tracking
_bibliography = OrderedDict()
