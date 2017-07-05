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
>>> citations = bib.report()
>>> citations.to_text()
JPL Horizons:
  implementation: 1996DPS....28.2504G

Classes
-------
Bib : A container for bibilography entries.

Functions
---------
register : Register a citation.
report   : Report registered citations.
status   : Bibliography tracking status.
stop     : Stop bibliography tracking.
track    : Start bibliography tracking.

"""

__all__ = ['Bib', 'register', 'report', 'status', 'stop', 'track']

from collections import OrderedDict

class Bib():
    """Bibliography class.

    References for specific tasks are provided as bibcode elements that can be queried at `ADS`_.

    .. _ADS: http://adsabs.harvard.edu
    """

    def __init__(self):
        self.bib = OrderedDict()

    def __setitem__(self, task, payload):
        """Set references for a specific task

        Parameters
        ----------
        task : str, mandatory
            name of the task generating the reference
        payload : dict, mandatory
            a dictionary with bibcodes for several aspects of the task, e.g., 
            one could cite one paper for the general `method` and one for the
            actual `implementation` used in SBPy
        """
        self.bib[task] = payload

    def __getitem__(self, ident):
        """Return the references for a specific task or for the n-th task""" 
        if isinstance(ident, str):
            return self.bib[ident]
        elif isinstance(ident, int):
            key = list(self.bib.keys())[idx]
            return {'task': key, 'references': self.bib[key]}
        
    def __iter__(self):
        return self.bib
            
    def __len__(self):
        return len(self.bib)
    
    @property
    def tasks(self):
        return list(self.bib.keys())
    
    def to_text(self):
        """convert bibcodes to human readable text

        not yet implemented
        """
        output = ''
        for ref in self.bib: #__iter__():
            output += '{:s}:\n'.format(ref)
            for key, val in self.__getitem__(ref).items():
                # transform val to ads(val)
                output += '  {:s}: {:s}\n'.format(key, val)
        return output

    def to_bibtex(self):
        """ convert bibcodes to LATeX bibtex

        not yet implemented
        """
        output = ''
        for ref in self.bib: #__iter__():
            for key, val in self.__getitem__(ref).items():
                # transform val to ads.bibtex(val)
                output += '% {:s}/{:s}:\n{:s}\n'.format(ref, key, val)
        return output

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

def report():
    """Report tracked `sbpy` citations."""
    return _bibliography

def reset():
    """Reset `sbpy` bibliography tracking."""
    global _bibliography
    _bibliography = Bib()

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

_track = False  # default is no bibliography tracking
_bibliography = Bib()
