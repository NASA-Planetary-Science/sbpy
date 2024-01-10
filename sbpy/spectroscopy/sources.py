# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy Spectroscopy Sources Module

Spectrophotometric classes that encapsulate synphot.SpectralSource and
synphot.Observation in order to generate sbpy spectra and photometry.

Requires synphot.

"""

__all__ = [
    'BlackbodySource', 'Reddening'
]

import numpy as np
from abc import ABC
import astropy.units as u
from astropy.utils.data import download_file, _is_url
from astropy.utils.decorators import deprecated

try:
    import synphot
    from synphot import SpectralElement, BaseUnitlessSpectrum
except ImportError:
    synphot = None

    class SpectralElement:
        pass

    class BaseUnitlessSpectrum:
        pass

from ..exceptions import SbpyException
from ..utils.decorators import requires

__doctest_requires__ = {
    'SpectralSource': ['synphot'],
    'BlackbodySource': ['synphot'],
}


class SinglePointSpectrumError(SbpyException):
    """Single point provided, but multiple values expected."""


@deprecated("v0.5.0")
class SynphotRequired(SbpyException):
    pass


class SpectralSource(ABC):
    """Abstract base class for SBPy spectral sources.

    Requires `~synphot`.


    Parameters
    ----------
    source : `~synphot.SourceSpectrum`
        The source spectrum.

    description : string, optional
        A brief description of the source spectrum.

    Attributes
    ----------
    wave        - Wavelengths of the source spectrum.
    fluxd       - Source spectrum.
    description - Brief description of the source spectrum.
    meta        - Meta data from ``source``, if any.

    """

    @requires("synphot")
    def __init__(self, source, description=None):
        self._source = source
        self._description = description

    @classmethod
    @requires("synphot")
    def from_array(cls, wave, fluxd, meta=None, **kwargs):
        """Create standard from arrays.


        Parameters
        ----------
        wave : `~astropy.units.Quantity`
            The spectral wavelengths.

        fluxd : `~astropy.units.Quantity`
            The spectral flux densities.

        meta : dict, optional
            Meta data.

        **kwargs
            Passed to object initialization.

        """

        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=wave, lookup_table=fluxd,
            meta=meta)

        return cls(source, **kwargs)

    @classmethod
    @requires("synphot")
    def from_file(cls, filename, wave_unit=None, flux_unit=None,
                  cache=True, **kwargs):
        """Load the source spectrum from a file.

        NaNs are dropped.


        Parameters
        ----------
        filename : string
            The name of the file.  See
            `~synphot.SourceSpectrum.from_file` for details.

        wave_unit, flux_unit : str or `~astropy.units.Unit`, optional
            Wavelength and flux units in the file.

        cache : bool, optional
            If ``True``, cache the contents of URLs.

        **kwargs
            Passed to object initialization.

        """

        if filename.lower().endswith(('.fits', '.fit', '.fz')):
            read_spec = synphot.specio.read_fits_spec
        else:
            read_spec = synphot.specio.read_ascii_spec

        # URL cache because synphot.SourceSpectrum.from_file does not
        if _is_url(filename):
            fn = download_file(filename, cache=True)
        else:
            fn = filename

        spec = read_spec(fn, wave_unit=wave_unit, flux_unit=flux_unit)
        i = np.isfinite(spec[1] * spec[2])
        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=spec[1][i], lookup_table=spec[2][i],
            meta={'header': spec[0]})

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
        return self._source

    @property
    def meta(self):
        self._source.meta

    def __call__(self, wave_or_freq, unit=None):
        """Evaluate/interpolate the source spectrum.


        Parameters
        ----------
        wave_or_freq : `~astropy.units.Quantity`
            Requested wavelengths or frequencies of the resulting
            spectrum.

        unit : string, `~astropy.units.Unit`, optional
            Spectral units of the output (flux density).  If ``None``,
            the default depends on ``wave_or_freq``: W/(m2 μm) for
            wavelengths, Jy for frequencies.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`
            The spectrum evaluated/interpolated to the requested
            wavelengths or frequencies.

        """

        from .. import units as sbu  # avoid circular dependency

        if unit is not None:
            unit = u.Unit(unit)
        elif wave_or_freq.unit.is_equivalent('m'):
            unit = u.Unit('W/(m2 um)')
        else:
            unit = u.Jy

        if unit.is_equivalent(sbu.VEGA):
            fluxd = self.source(wave_or_freq, 'W/(m2 um)').to(
                unit, sbu.spectral_density_vega(wave_or_freq))
        else:
            fluxd = self.source(wave_or_freq, unit)

        return fluxd

    def observe(self, wfb, unit=None, interpolate=False, **kwargs):
        """Observe source as through filters or spectrometer.

        Calls ``observe_bandpass``, ``observe_spectrum``, or
        ``self()``, as appropriate.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`
            Wavelengths, frequencies, or bandpasses.  May also be a
            list of ``SpectralElement``s.

        unit : string, `~astropy.units.Unit`, optional
            Spectral flux density units for the output.  If ``None``,
            the default depends on ``wfb``: W/(m2 μm) for wavelengths
            or bandpasses, Jy for frequencies.

        interpolate : bool, optional
            For wavelengths/frequencies, set to ``True`` for
            interpolation instead of rebinning.  Use this when the
            spectral resolution of the source is close to that of the
            requested wavelengths.

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`
            The spectrum rebinned.


        Raises
        ------
        SinglePointSpectrumError - If requested wavelengths or
            frequencies has only one value.


        Notes
        -----
        Method for spectra adapted from AstroBetter post by Jessica Lu:
        https://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/

        """

        if isinstance(wfb, (list, tuple, SpectralElement)):
            lambda_eff, fluxd = self.observe_bandpass(
                wfb, unit=unit, **kwargs)
        elif isinstance(wfb, u.Quantity):
            if interpolate:
                fluxd = self(wfb, unit=unit)
            else:
                fluxd = self.observe_spectrum(wfb, unit=unit, **kwargs)
        else:
            raise TypeError('Unsupported type for `wfb` type: {}'
                            .format(type(wfb)))

        return fluxd

    def observe_bandpass(self, bp, unit=None, **kwargs):
        """Observe through a bandpass.


        Parameters
        ----------
        bp : `~synphot.SpectralElement`, list, or tuple
            Bandpass.

        unit : string, `~astropy.units.Unit`, optional
            Spectral flux density units for the output.  The default
            is W/(m2 μm).

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        lambda_eff : `~astropy.units.Quantity`
            Effective wavelength(s) of the observation(s).

        fluxd : `~astropy.units.Quantity`
            The spectrum rebinned.

        """

        from .. import units as sbu  # avoid circular dependency

        # promote single bandpasses to a list, but preserve number of
        # dimensions
        if isinstance(bp, (SpectralElement, str)):
            ndim = 0
            bp = [bp]
        else:
            ndim = np.ndim(bp)

        if unit is None:
            unit = u.Unit('W/(m2 um)')
        else:
            unit = u.Unit(unit)

        fluxd = np.ones(len(bp)) * unit
        for i in range(len(bp)):
            obs = synphot.Observation(self.source, bp[i], **kwargs)
            lambda_eff = obs.effective_wavelength()
            lambda_pivot = obs.bandpass.pivot()
            _fluxd = obs.effstim('W/(m2 um)')

            if unit.is_equivalent(sbu.VEGAmag):
                fluxd[i] = _fluxd.to(unit, sbu.spectral_density_vega(bp[i]))
            else:
                fluxd[i] = _fluxd.to(unit, u.spectral_density(lambda_pivot))

        if np.ndim(fluxd) != ndim:
            fluxd = fluxd.squeeze()

        return lambda_eff, fluxd

    def observe_spectrum(self, wave_or_freq, unit=None, **kwargs):
        """Observe source as through a spectrometer.

        .. Important:: This method works best when the requested
           spectral resolution is lower than the spectral resolution
           of the internal data.  If the requested
           wavelengths/frequencies are exactly the same as the
           internal spectrum, then the internal spectrum will be
           returned without binning.  This special case does not work
           for subsets of the wavelengths.


        Parameters
        ----------
        wave_or_freq : `~astropy.units.Quantity`
            Wavelengths or frequencies of the spectrum.  Spectral bins
            will be centered at these values.  The length must be
            larger than 1.

        unit : string, `~astropy.units.Unit`, optional
            Spectral flux density units for the output.  If ``None``,
            the default is W/(m2 μm) for wavelengths, Jy for
            frequencies.

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`
            The spectrum rebinned.


        Raises
        ------
        SinglePointSpectrumError - If requested wavelengths or
            frequencies has only one value.


        Notes
        -----
        Method for spectra adapted from AstroBetter post by Jessica Lu:
        https://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/

        """

        from .. import units as sbu  # avoid circular dependency

        if np.size(wave_or_freq) == 1:
            raise SinglePointSpectrumError(
                'Multiple wavelengths or frequencies required for '
                'observe_spectrum.  Instead consider interpolation '
                'with {}().'
                .format(self.__class__.__name__))

        if unit is None:
            if wave_or_freq.unit.is_equivalent('m'):
                unit = u.Unit('W/(m2 um)')
            else:
                unit = u.Jy
        else:
            unit = u.Unit(unit)

        specele = synphot.SpectralElement(synphot.ConstFlux1D, amplitude=1)

        # Specele is defined over all wavelengths, but most spectral
        # standards are not.  force='taper' will affect retrieving
        # flux densities at the edges of the spectrum, but is
        # preferred to avoid wild extrapolation.
        kwargs['force'] = kwargs.get('force', 'taper')

        obs = synphot.Observation(
            self.source, specele, binset=wave_or_freq, **kwargs)

        if unit.is_equivalent(sbu.VEGAmag):
            fluxd = obs.sample_binned(flux_unit='W/(m2 um)').to(
                unit, sbu.spectral_density_vega(wave_or_freq))
        else:
            fluxd = obs.sample_binned(flux_unit=unit)

        return fluxd

    def color_index(self, wfb, unit):
        """Color index (magnitudes) and effective wavelengths.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity` or tuple of `~synphot.SectralElement`
            Two wavelengths, frequencies, or bandpasses.

        unit : string or `~astropy.units.MagUnit`
            Units for the calculation, e.g., ``astropy.units.ABmag`` or
            ``sbpy.units.VEGAmag``.


        Returns
        -------
        eff_wave : `~astropy.units.Quantity`
            Effective wavelengths for each ``wfb``.

        ci : `~astropy.units.Quantity`
            Color index, ``m_0 - m_1``, where 0 and 1 are element
            indexes for ``wfb``.

        """

        eff_wave = []
        m = np.zeros(2) * u.Unit(unit)
        for i in range(2):
            if isinstance(wfb[i], u.Quantity):
                if wfb[i].unit.is_equivalent(u.Hz):
                    eff_wave.append(wfb[i].to(u.um, u.spectral()))
                else:
                    eff_wave.append(wfb[i])
                m[i] = self(eff_wave[i], unit=unit)
            elif isinstance(wfb[i], (list, tuple, SpectralElement)):
                w, m[i] = self.observe_bandpass(wfb[i], unit=unit)
                eff_wave.append(w)
            else:
                raise TypeError('Unsupported type for `wfb` type: {}'
                                .format(type(wfb[i])))

        ci = m[0] - m[1]

        return u.Quantity(eff_wave), ci

    def redden(self, S):
        """Redden the spectrum.

        Parameters
        ----------
        S : `~SpectralGradient`
            The spectral gradient to redden.

        Returns
        -------
        spec : `~SpectralSource`
            Reddened spectrum

        """
        from copy import deepcopy
        r = Reddening(S)
        red_spec = deepcopy(self)
        red_spec._source = red_spec.source * r
        if red_spec.description is not None:
            red_spec._description = '{} reddened by {} at {}'.format(
                red_spec.description, S, S.wave0)
        return red_spec


class BlackbodySource(SpectralSource):
    """Blackbody sphere.

    Spectral flux densities are calculated from ``pi * B(T)``, where
    ``B`` is the Planck function.


    Parameters
    ----------
    T : `~astropy.units.Quantity`, required
        Temperature in Kelvin.

    """

    @requires("synphot")
    def __init__(self, T=None):
        super().__init__(None, description='πB(T)')

        if T is None:
            raise TypeError('T is required.')

        self._T = u.Quantity(T, u.K)
        self._source = synphot.SourceSpectrum(
            synphot.BlackBody1D, temperature=self._T.value) * np.pi

    def __repr__(self):
        return '<BlackbodySource: T={}>'.format(self._T)

    @property
    def T(self):
        return self._T


class Reddening(BaseUnitlessSpectrum):
    """Class to handle simple linear reddening.

    Parameters
    ----------
    S : `~SpectralGradient`
        The spectral gradient to redden.

    """

    @u.quantity_input(S=u.percent / u.um)
    @requires("synphot")
    def __init__(self, S):

        if getattr(S, 'wave0', None) is None:
            raise ValueError("Normalization wavelength in `S` (.wave0) is "
                             "required by not available.")

        wv = [1, 2] * S.wave0
        df = (S.wave0 * S).to('').value

        super().__init__(
            synphot.Empirical1D, points=wv, lookup_table=[1, 1+df],
            fill_value=None)
