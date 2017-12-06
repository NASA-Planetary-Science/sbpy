# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=====================
SBPy Vega Core Module
=====================

Classes
-------
Vega
Bohlin2014    - Vega spectrum source classes.

default_vega  - Controls the default Vega spectrum in `sbpy`.

"""

from astropy.utils.state import ScienceState
import synphot  # synphot is required

__all__ = [
#    'Bohlin2014',
    'Vega',
    'default_vega'
]

available_spectra = ['Bohlin2014']

class Vega:
    """Vega spectrum.

    Parameters
    ----------
    source : `~synphot.SourceSpectrum`
      The source spectrum for Vega, flux density / spectral irradiance
      as seen from Earth.

    description : string, optional
      A description of the source spectrum.

    bibcode : string, optional
      The bibliography code for this spectrum.

    Attributes
    ----------
    description - Description of the source spectrum.
    wave        - Wavelengths of the source spectrum.
    fluxd       - The source spectrum.


    Examples
    --------

    """
    
    def __init__(self, source, description=None, bibcode=None):
        self._source = source
        self._description = description
        self._bibcode = bibcode

    @classmethod
    def from_file(cls, filename, wave_unit='um', flux_unit='W/(m2 um)',
                  **kwargs):
        """Load the source spectrum from a file.

        Parameters
        ----------
        filename : string
          The name of the file.  The file must be compatible with
          `~synphot.SourceSpectrum.from_file`.

        wave_unit, flux_unit : str or `~astropy.units.core.Unit`
          Wavelength and flux units.

        **kwargs
          Passed to `Vega` initialization.

        """
        
        source = synphot.SourceSpectrum.from_file(
            filename, wave_unit=wave_unit, flux_unit=flux_unit)
        return cls(source, **kwargs)
        
    @property
    def description(self):
        """Description of the source spectrum."""
        return self._description

    @property
    def wave(self):
        """Wavelengths of the source spectrum."""
        return self.source.waveset

    @property
    def fluxd(self):
        """The source spectrum."""
        return self.source(self.source.waveset, flux_unit='W / (m2 um)')

    @property
    def source(self):
        if self._bibcode is not None:
            bib.register('spectroscopy.vega', {description, self._bibcode})
        return self._source

 
#class Bohlin2014(Vega):
#    """Vega spectrum from Bohlin 2014.
#
#    References
#    ----------
#    Bohlin, R. C. 2014, Astronomical Journal, 147, 127.
#
#    """
#
#    def __init__(self):
#        from .. import bib
#        
#        description = 'Vega spectrum from Bohlin (2014).'
#        sp = synphot.SourceSpectrum.from_file(fn)
#
#        bib.register('sbpy.data.sun.Sun',
#                     {'Vega spectrum': '2014AJ....147..127B'})
#        
#        super().__init__(self, sp, description)

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
        
        try:
            vega = getattr(sys.modules[__name__], arg)
        except AttributeError:
            msg = 'Unknown Vega spectrum "{}".  Valid spectra:\n{}'.format(arg, available_spectra)
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
