# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
SBPy Sun Core Module
====================
"""

import numpy as np
from astropy.utils.state import ScienceState
from astropy.utils.data import _is_url
import synphot
from ... import bib

__all__ = ['Sun', 'default_sun']

class Sun:
    """Solar spectrum class.

    Parameters
    ----------
    source : `~synphot.SourceSpectrum`
      The source spectrum, flux density / spectral irradiance at 1 au.

    description : string, optional
      A description of the source spectrum.

    bibcode : string, optional
      Bibliography code for `sbpy.bib.register`.

    Attributes
    ----------
    description - Description of the source spectrum.
    wave        - Wavelengths of the source spectrum.
    fluxd       - The source spectrum.


    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import astropy.constants as const
    >>> from synphot import SourceSpectrum, BlackBody1D
    >>> from sbpy.spectroscopy.sun import Sun
    >>> scale = (np.pi * const.R_sun**2 / const.au**2).decompose().value
    >>> source = SourceSpectrum(BlackBody1D, temperature=5770 * u.K) * scale
    >>> sun = Sun(source)

    """
    
    def __init__(self, source, description=None, bibcode=None):
        self._source = source
        self._description = description
        self._bibcode = bibcode

    @classmethod
    def from_array(cls, wave, fluxd, meta=None, **kwargs):
        """Load the source spectrum from an array.

        Parameters
        ----------
        wave : `~astropy.units.Quantity`
          The spectral wavelengths.

        fluxd : `~astropy.units.Quantity`
          The solar flux densities, at 1 au.

        meta : dict, optional
          `synphot.SourceSpectrum` meta data.

        **kwargs
          Passed to `Sun` initialization.

        """
        
        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=wave, lookup_table=spec[2],
            meta=meta)
        
        return cls(source, **kwargs)

    @classmethod
    def from_file(cls, filename, wave_unit='um', flux_unit='W/(m2 um)',
                  cache=True, **kwargs):
        """Load the source spectrum from a file.

        Parameters
        ----------
        filename : string
          The name of the file.  The file must be compatible with
          `~synphot.SourceSpectrum.from_file`.  If it ends with
          '.fits' or '.fit', it is read as a FITS file, otherwise,
          text is assumed.

        wave_unit, flux_unit : str or `~astropy.units.core.Unit`
          Wavelength and flux units.

        cache : bool, optional
          If `True`, cache the contents of URLs.

        **kwargs
          Passed to `Sun` initialization.

        """

        from astropy.utils.data import download_file
        from synphot.specio import read_fits_spec, read_ascii_spec
        
        # URL cache because synphot.SourceSpectrum.from_file does not
        if _is_url(filename):
            if filename.lower().endswith(('.fits', '.fit')):
                read_spec = read_fits_spec
            else:
                read_spec = read_ascii_spec

            fn = download_file(filename, cache=True)
            spec = read_spec(fn, wave_unit=wave_unit, flux_unit=flux_unit)
            source = synphot.SourceSpectrum(
                synphot.Empirical1D, points=spec[1], lookup_table=spec[2],
                meta={'header': spec[0]})
        else:
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
        return self.source(self._source.waveset, flux_unit='W / (m2 um)')

    @property
    def source(self):
        if self._bibcode is not None:
            bib.register('spectroscopy.sun', {self._description, self._bibcode})
        return self._source
            
    def rebin(self, wave, unit='W / (m2 um)'):
        """Rebinned source spectrum.

        Parameters
        ----------
        wave : `~astropy.units.Quantity`
          Return solar spectrum at these wavelengths.

        unit : string, `~astropy.units.Unit`, optional
          Spectral units of the output: flux density, 'vegamag',
          'ABmag', or 'STmag'.

        Returns
        -------
        fluxd : `~astropy.units.Quantity`
          Spectrum of the Sun.  If multiple wavelengths are requested,
          the solar spectrum will be binned to match `wave`.
          Otherwise, the spectrum at that wavelength will be returned.

        """

        if np.size(wave) > 1:
            # Method adapted from http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
            specele = synphot.SpectralElement(synphot.ConstFlux1D(1))
            obs = synphot.Observation(self.source, specele, binset=wave,
                                      force='extrap')
            fluxd = obs.sample_binned(wave, unit)
        else:
            fluxd = self.source(wave, unit)

        return fluxd

    def filt(self, bp, unit='W / (m2 um)', vegaspec=None, **kwargs):
        """Solar spectrum observed through a filter.

        Parameters
        ----------
        bp : string or `~synphot.SpectralElement`
          The name of a filter, or a transmission spectrum as a
          `~synphot.SpectralElement`.  See notes for built-in filter
          names.

        unit : string, `~astropy.units.Unit`, optional
          Spectral units of the output: flux density, 'vegamag',
          'ABmag', or 'STmag'.

        vegaspec : `~synphot.SourceSpectrum`
          Use this spectrum for Vega when `unit == 'vegamag'`,
          otherwise the default `sbpy` spectrum will be used.

        **kwargs
          Additional keyword arguments for
          `~synphot.observation.Observation`.

        Returns
        -------
        wave : `~astropy.units.Quantity`
          Effective wavelength.

        fluxd : `~astropy.units.Quantity`
          The flux density.


        Notes
        -----

        Filter reference data is from STScI's Calibration Reference
        Data System.

        * ``'bessel_j'`` (Bessel *J*)
        * ``'bessel_h'`` (Bessel *H*)
        * ``'bessel_k'`` (Bessel *K*)
        * ``'cousins_r'`` (Cousins *R*)
        * ``'cousins_i'`` (Cousins *I*)
        * ``'johnson_u'`` (Johnson *U*)
        * ``'johnson_b'`` (Johnson *B*)
        * ``'johnson_v'`` (Johnson *V*)
        * ``'johnson_r'`` (Johnson *R*)
        * ``'johnson_i'`` (Johnson *I*)
        * ``'johnson_j'`` (Johnson *J*)
        * ``'johnson_k'`` (Johnson *K*)

        """
        
        assert isinstance(bp, (str, synphot.SpectralElement))

        if isinstance(bp, str):
            bp = get_bandpass(bp)

        obs = synphot.Observation(self.source, bp, **kwargs)

        if unit == 'vegamag' and vegaspec is None:
            from ..vega import default_vega
            vegaspec = default_vega.get().source

        wave = obs.effective_wavelength()
        fluxd = obs.effstim(unit, vegaspec=vegaspec)

        return wave, fluxd

class default_sun(ScienceState):
    """The default solar spectrum to use.

    To retrieve the current default::

      >>> sun = default_sun.get()

    To change it::

      >>> from sbpy.spectroscopy.sun import default_sun
      >>> with default_sun.set('E490_2014'):
      ...     # E490_2014 in effect
      ...     pass

    Or, you may use a string::

      >>> with default_sun.set('E490_2014LR'):
      ...     # E490_2014LR in effect
      ...     pass

    """
    _value = 'E490_2014'

    @staticmethod
    def get_sun_from_string(arg):
        """Return a Sun instance from a string."""
        import os
        import sys
        from . import sources
        
        try:
            parameters = getattr(sources, arg).copy()
            filename = parameters.pop('file')

            if not _is_url(filename):
                path = os.path.dirname(__file__)
                filename = os.sep.join((path, filename))
                
            sun = Sun.from_file(filename, **parameters)
        except AttributeError:
            msg = 'Unknown solar spectrum "{}".  Valid spectra:\n{}'.format(arg, available_spectra)
            raise ValueError(msg)

        return sun

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return cls.get_sun_from_string(value)
        elif isinstance(value, Sun):
            return value
        else:
            raise TypeError("default_sun must be a string or Sun instance.")

