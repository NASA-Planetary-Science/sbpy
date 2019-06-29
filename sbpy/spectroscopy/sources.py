# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy Spectroscopy Sources Module

Spectrophotometric classes that encasulate synphot.SpectralSource and
synphot.Observation in order to generate sbpy spectra and photometry.

Requires synphot.

"""

__all__ = [
    'BlackbodySource'
]

from abc import ABC
import numpy as np
import astropy.units as u
from astropy.utils.data import download_file, _is_url

try:
    import synphot
except ImportError:
    synphot = None

from ..exceptions import SbpyException
from ..bib import register


class SinglePointSpectrumError(SbpyException):
    '''Single point provided, but multiple values expected.'''


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

    def __init__(self, source, description=None):
        if synphot is None:
            raise ImportError(
                'synphot required for {}.'.format(cls.__name__))

        self._source = source
        self._description = description

    @classmethod
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
        spectral_source = cls(None, **kwargs)
        spectral_source._source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=wave, lookup_table=fluxd,
            meta=meta)

        return spectral_source

    @classmethod
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

        spectral_source = cls(None, **kwargs)

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
        spectral_source._source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=spec[1][i], lookup_table=spec[2][i],
            meta={'header': spec[0]})

        return spectral_source

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
        """Evaluate or interpolate the source spectrum.


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
            The spectrum evaluated at or interpolated to the requested
            wavelengths or frequencies.

        """

        from .. import units as sbu  # avoid circular dependency

        if unit is None:
            if wave_or_freq.unit.is_equivalent('m'):
                unit = u.Unit('W/(m2 um)')
            else:
                unit = u.Jy
        else:
            unit = u.Unit(unit)

        if unit.is_equivalent(sbu.VEGA):
            fluxd = self.source(wave_or_freq, 'W/(m2 um)').to(
                unit, sbu.spectral_density_vega(wave_or_freq))
        else:
            fluxd = self.source(wave_or_freq, unit)

        return fluxd

    def observe(self, wfb, unit=None, **kwargs):
        """Observe source as through filters or spectrometer.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`
            Wavelengths, frequencies, or bandpasses.  May also be a
            list of ``SpectralElement``s.

        unit : string, `~astropy.units.Unit`, optional
            Spectral units of the output (flux density).  If ``None``,
            the default depends on ``wfb``: W/(m2 μm) for wavelengths
            or bandpasses, Jy for frequencies.

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
        Method adapted from AstroBetter post by Jessica Lu:
        http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/

        """

        if synphot is None:
            raise ImportError(
                'synphot required for {}.observe.'.format(self.__name__))

        from .. import units as sbu  # avoid circular dependency

        # promote single bandpasses to a list, but preserve number of
        # dimensions
        if isinstance(wfb, (synphot.SpectralElement, str)):
            ndim = 0
            wfb = [wfb]
        else:
            ndim = np.ndim(wfb)

        if isinstance(wfb, (tuple, list)):
            if unit is None:
                unit = u.Unit('W/(m2 um)')
            else:
                unit = u.Unit(unit)

            fluxd = np.ones(len(wfb)) * unit
            for i in range(len(wfb)):
                fluxd[i] = self.filt(wfb[i], unit=unit, **kwargs)[2]
        else:
            if np.size(wfb) == 1:
                raise SinglePointSpectrumError(
                    'Multiple wavelengths or frequencies required for '
                    'observe.  Consider interpolation with {}() instead.'
                    .format(self.__class__.__name__))

            if unit is None:
                if wfb.unit.is_equivalent('m'):
                    unit = u.Unit('W/(m2 um)')
                else:
                    unit = u.Jy
            else:
                unit = u.Unit(unit)

            specele = synphot.SpectralElement(synphot.ConstFlux1D(1))

            # Use force='taper' to prevent PartialOverlap execption.
            # Specele is defined over all wavelengths, but most
            # spectral standards are not.
            kwargs['force'] = kwargs.get('force', 'taper')

            obs = synphot.Observation(
                self.source, specele, binset=wfb, **kwargs)

            if unit.is_equivalent(sbu.VEGAmag):
                fluxd = obs.sample_binned(flux_unit='W/(m2 um)').to(
                    unit, sbu.spectral_density_vega(wfb))
            else:
                fluxd = obs.sample_binned(flux_unit=unit)

        if np.ndim(fluxd) != ndim:
            # likely need a squeeze
            fluxd = fluxd.squeeze()

        return fluxd

    def filt(self, bp, unit='W / (m2 um)', **kwargs):
        """Observe the source through a single filter.

        Parameters
        ----------
        bp: string or `~synphot.SpectralElement`
            The name of a filter, or a transmission spectrum as a
            `~synphot.SpectralElement`.  See Notes for built-in filter
            names.

        unit: string, `~astropy.units.Unit`, optional
            Spectral flux density units of the output.

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        lambda_eff: `~astropy.units.Quantity`
            Effective wavelength.

        lambda_pivot: `~astropy.units.Quantity`
            Pivot wavelength.

        fluxd: `~astropy.units.Quantity`
            Spectral flux density.


        Notes
        -----

        Filter reference data is from STScI's Calibration Reference
        Data System.

        * ``'bessel_j'`` (Bessel * J*)
        * ``'bessel_h'`` (Bessel * H*)
        * ``'bessel_k'`` (Bessel * K*)
        * ``'cousins_r'`` (Cousins * R*)
        * ``'cousins_i'`` (Cousins * I*)
        * ``'johnson_u'`` (Johnson * U*)
        * ``'johnson_b'`` (Johnson * B*)
        * ``'johnson_v'`` (Johnson * V*)
        * ``'johnson_r'`` (Johnson * R*)
        * ``'johnson_i'`` (Johnson * I*)
        * ``'johnson_j'`` (Johnson * J*)
        * ``'johnson_k'`` (Johnson * K*)

        """

        if synphot is None:
            raise ImportError(
                'synphot required for {}.filt.'.format(cls.__name__))

        from .. import units as sbu  # avoid circular dependency

        if isinstance(bp, str):
            bp = synphot.SpectralElement.from_filter(bp)

        obs = synphot.Observation(self.source, bp, **kwargs)
        lambda_eff = obs.effective_wavelength()
        lambda_pivot = obs.bandpass.pivot()
        _unit = u.Unit(unit)

        if _unit.is_equivalent(sbu.VEGAmag):
            fluxd = obs.effstim('W/(m2 um)').to(
                _unit, sbu.spectral_density_vega(bp))
        else:
            fluxd = obs.effstim(flux_unit=_unit)

        return lambda_eff, lambda_pivot, fluxd

    def color_index(self, wfb, unit, equivalencies=[], equiv_func=None):
        """Color index (magnitudes) and effective wavelengths.


        Parameters
        ----------
        wfb : two-element `~astropy.units.Quantity` or tuple
            Wavelengths, frequencies, or bandpasses of the
            measurement.  See :func:`~observe`.

        unit : string or `~astropy.units.MagUnit`
            Units for the output, e.g., ``astropy.units.ABmag`` or
            ``sbpy.units.VEGAmag``.

        equivalencies : list, optional
            List of unit equivalencies for magnitude-flux-density
            conversion.

        equiv_func : callable, optional
            A function that takes the effective wavelength of the
            observation and generates an equivalency list, e.g., use
            :func:`~sbpy.units.spectral_density_vega` for Vega-based
            magnitude conversions.
            :func:`~sbpy.units.spectral_density` is automatically
            included.


        Returns
        -------
        eff_wave : `~astropy.units.Quantity`
            Effective wavelengths for each ``wfb``.

        ci : `~astropy.units.Quantity`
            Color index, ``m_0 - m_1``, where 0 and 1 are element
            indexes for ``wfb``.

        """

        eff_wave = np.zeros(2) * u.um
        m = np.zeros(2) * u.Unit(unit)
        for i in range(2):
            if isinstance(wfb[i], u.Quantity):
                eff_wave[i] = wfb[i].to(u.um, u.spectral())
                f = self(eff_wave[i])
            else:
                eff_wave[i], f = self.filt(wfb[i])

            eqv = u.spectral_density(eff_wave[i])
            if equiv_func is not None:
                eqv.extend(equiv_func(eff_wave[i]))

            m[i] = f.to(m.unit, equivalencies + eqv)
        ci = m[0] - m[1]

        return eff_wave, ci


class BlackbodySource(SpectralSource):
    """Blackbody sphere.

    Spectral flux densities are calculated from ``pi * B(T)``, where
    ``B`` is the Planck function.


    Parameters
    ----------
    T : `~astropy.units.Quantity`, required
        Temperature in Kelvin.

    """

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


class SpectralStandard(SpectralSource, ABC):
    """Abstract base class for SBPy spectral standards.


    Parameters
    ----------
    source : `~synphot.SourceSpectrum`
        The source spectrum.

    description : string, optional
        A brief description of the source spectrum.

    bibcode : string, optional
        Bibliography code for `sbpy.bib.register`.


    Attributes
    ----------
    wave        - Wavelengths of the source spectrum.
    fluxd       - Source spectrum.
    description - Brief description of the source spectrum.
    meta        - Meta data from `source`, if any.

    """

    def __init__(self, source, description=None, bibcode=None):
        self._source = source
        self._description = description
        self._bibcode = bibcode

    @property
    def source(self):
        if self._bibcode is not None:
            register('spectroscopy', {self._description: self._bibcode})
        return self._source

    def filt(self, bp, unit='W / (m2 um)', source_fluxd={}, **kwargs):
        """Observe the standard through a single filter.

        Uses the `sbpy` calibration system.  If `bp` has been set in
        `~sbpy.calib.solar_fluxd`, then those values will be returned.
        Otherwise, the solar spectrum defined by this object will be
        filtered with the provided bandpass.


        Parameters
        ----------
        bp: string or `~synphot.SpectralElement`
            The name of a filter, or a transmission spectrum as a
            `~synphot.SpectralElement`.  See
            `sbpy.spectroscopy.sources.SpectralSource` for built-in
            filter names.

        unit: string, `~astropy.units.Unit`, optional
            Spectral flux density units of the output.

        source_fluxd : dict
            Specific flux densities and wavelengths to use for filter
            bandpasses specified as strings in `bp`.  If the name of
            the filter is BP, then the keys would be:

                BP : flux density
                BP_lambda_eff : effective wavelength
                BP_lambda_pivot : pivot wavelength

            Multiple bandpasses may be specified in the same
            dictionary.

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        lambda_eff: `~astropy.units.Quantity`
            Effective wavelength.  ``None`` if it cannot be calculated.

        lambda_pivot: `~astropy.units.Quantity`
            Pivot wavelength.  ``None`` if it cannot be calculated.

        fluxd: `~astropy.units.Quantity`
            Spectral flux density.

        """

        if bp in source_fluxd:
            fluxd = source_fluxd[bp]
            lambda_eff = source_fluxd.get(bp + '_lambda_eff')
            lambda_pivot = source_fluxd.get(bp + '_lambda_pivot')
            if lambda_pivot is None:
                equiv = None
            else:
                equiv = [u.spectral_density(lambda_pivot)
                         + u.spectral_density_vega(lambda_pivot)]
            try:
                fluxd = fluxd.to(unit, equiv)
            except u.UnitConversionError as e:
                raise type(e)(
                    '{}  Is "{}_lambda_pivot" required and'
                    ' was it provided?'.format(e.message, bp))
        else:
            lambda_eff, lambda_pivot, fluxd = (
                super().filt(bp, unit=unit, **kwargs))

        return lambda_eff, lambda_pivot, fluxd
