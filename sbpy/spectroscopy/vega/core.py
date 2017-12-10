# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=====================
SBPy Vega Core Module
=====================
"""

from astropy.utils.state import ScienceState
from ..core import SpectralStandard

__doctest_requires__ = {'Sun': 'synphot'}

__all__ = [
    'Vega',
    'default_vega'
]

class Vega(SpectralStandard):
    """Vega spectrum.

    Parameters
    ----------
    wave : `~astropy.units.Quantity`
      The spectral wavelengths.

    fluxd : `~astropy.units.Quantity`
      The solar flux densities, at 1 au.

    description : string, optional
      A brief description of the source spectrum.

    bibcode : string, optional
      Bibliography code for `sbpy.bib.register`.

    meta : dict, optional
      Any additional meta data, passed on to
      `~synphot.SourceSpectrum`.


    Attributes
    ----------
    wave        - Wavelengths of the source spectrum.
    fluxd       - Source spectrum.
    description - Brief description of the source spectrum.
    meta        - Meta data.


    Examples
    --------
    Create Vega from a file::
      >>> vega = Vega.from_file('filename')                     # doctest: +SKIP

    Evaluate Vega at 1 μm::
      >>> print(vega(1 * u.um))                   # doctest: +SKIP

    """
    pass

class default_vega(ScienceState):
    """The default Vega spectrum to use.

    To change it::

      >>> from sbpy.spectroscopy.vega import default_vega
      >>> with default_vega(Vega.from_file(filename))  # doctest: +SKIP
      ...     # Vega from filename in effect

    """
    _value = 'Bohlin2014'

    @staticmethod
    def get_vega_from_string(arg):
        """Return a Vega instance from a string."""

        import sys
        from . import sources
        
        try:
            parameters = getattr(sources, arg).copy()
            vega = Vega.from_file(**parameters)
        except AttributeError:
            msg = 'Unknown Vega spectrum "{}".  Valid spectra:\n{}'.format(arg, sources.keys())
            raise ValueError(msg)

        return vega

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return cls.get_vega_from_string(value)
        elif isinstance(value, Vega):
            return value
        else:
            raise TypeError("default_vega must be a string or Vega instance.")
