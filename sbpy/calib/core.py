# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = [
    'solar_spectrum',
    'vega_spectrum',
    'solar_fluxd',
    'vega_fluxd',
    'Sun',
    'Vega',
    'UndefinedSourceError'
]

import os
from abc import ABC
import numpy as np
from astropy.utils.state import ScienceState
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
import astropy.units as u
from ..spectroscopy.sources import SpectralSource
from ..exceptions import SbpyException, OptionalPackageUnavailable
from .. import bib
from . import solar_sources, vega_sources

try:
    import synphot
    from synphot import SpectralElement
except ImportError:
    synphot = None

    class SpectralElement:
        pass

__doctest_requires__ = {'Sun': 'synphot'}


class UndefinedSourceError(SbpyException):
    "SpectralStandard was initialized without a source, but it was accessed."


class SpectralStandard(SpectralSource, ABC):
    """Abstract base class for SBPy spectral standards.


    Parameters
    ----------
    spectrum_state : `~astropy.utils.state.ScienceState`
        Context manager for this source's spectrum.

    fluxd_state : `~astropy.utils.state.ScienceState`
        Context manager for this source's calibration by filter.

    source : `~synphot.SourceSpectrum` or ``None``
        The source spectrum or ``None`` if unspecified.

    description : string, optional
        A brief description of the source spectrum.

    bibcode : string, optional
        Bibliography code for citation (see `sbpy.bib.register`).


    Attributes
    ----------
    wave        - Wavelengths of the source spectrum.
    fluxd       - Source spectrum.
    description - Brief description of the source spectrum.
    meta        - Meta data from `source`, if any.


    Notes
    -----

    """

    def __init__(self, source, description=None, bibcode=None):
        self._source = source
        self._description = description
        self._bibcode = bibcode
        self._bibtask = '.'.join((self.__module__, self.__class__.__name__))

    def __repr__(self):
        if self.description is None:
            return '<{}>'.format(self.__class__.__name__)
        else:
            return '<{}: {}>'.format(self.__class__.__name__,
                                     self.description)

    @classmethod
    def from_builtin(cls, name):
        """Spectrum from a built-in `sbpy` source.

        Parameters
        ----------
        name : string
            The name of the spectrum.  See `show_builtin()` for
            available sources.

        """

        from astropy.utils.data import _is_url

        try:
            parameters = getattr(cls._sources, name).copy()
        except AttributeError:
            msg = 'Unknown spectrum "{}".  Valid spectra:\n'.format(
                name) + cls.show_builtin(print=False)
            raise UndefinedSourceError(msg)

        if not _is_url(parameters['filename']):
            # find in the module's location
            parameters['filename'] = get_pkg_data_filename(
                os.path.join('data', parameters['filename']))

        return cls.from_file(**parameters)

    @classmethod
    def from_default(cls):
        """Initialize new spectral standard from current default.

        The spectrum will be ``None`` if `synphot` is not available.

        """
        if synphot is None:
            warnings.warn(OptionalPackageUnavailable(
                'synphot is not installed, returning an empty spectral'
                ' standard.'))
            standard = cls(None)
        else:
            standard = cls._spectrum_state.get()
        return standard

    @classmethod
    def show_builtin(cls, print=True):
        """List built-in spectra."""
        from builtins import print as print_func
        rows = []
        for name in cls._sources.available:
            source = getattr(cls._sources, name)
            rows.append((name, source['description']))
        tab = Table(rows=rows, names=('name', 'description')).pformat(
            max_lines=-1, max_width=-1)
        tab = '\n'.join(tab)

        if print:
            print_func(tab)
        else:
            return tab

    @property
    def source(self):
        if self._source is None:
            raise UndefinedSourceError('The source is not defined.')

        if self._bibcode is not None:
            bib.register(self._bibtask, {self._description: self._bibcode})

        return self._source

    def observe(self, wfb, unit=None, **kwargs):
        """Observe as through filters or spectrometer.

        Calls `observe_bandpass`, `observe_spectrum`, or
        `observe_filter` as appropriate.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, string
            Wavelengths, frequencies, or bandpasses.  Bandpasses may
            be a filter name (string).  May also be a list of
            ``SpectralElement``s or strings.

        unit : string, `~astropy.units.Unit`, optional
            Units of the output (spectral flux density).

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`

        """

        if isinstance(wfb, (list, tuple)):
            fluxd = []
            for i in range(len(wfb)):
                fluxd.append(self.observe(wfb[i], unit=unit, **kwargs))
            fluxd = u.Quantity(fluxd)
        elif isinstance(wfb, str):
            lambda_eff, lambda_pivot, fluxd = self.observe_filter(
                wfb, unit=unit)
        else:
            fluxd = super().observe(wfb, unit=unit, **kwargs)

        return fluxd

    def observe_filter(self, filt, unit=None):
        """Flux density through this filter.

        Does not use the spectrum, but instead the flux density
        calibration manager.  If the name of the filter is BP, then
        the expected keys are:

            BP : flux density
            BP(lambda eff) : effective wavelength, optional
            BP(lambda pivot) : pivot wavelength, optional


        Parameters
        ----------
        filt : string
            Name of the filter.

        unit : string, `~astropy.units.Unit`, optional
            Spectral flux density units for the output.


        Returns
        -------
        lambda_eff: `~astropy.units.Quantity`
            Effective wavelength.  ``None`` if it is not provided.

        lambda_pivot: `~astropy.units.Quantity`
            Pivot wavelength.  ``None`` if it is not provided.

         fluxd : `~astropy.units.Quantity`
            Spectral flux density.


        Raises
        ------
        ``KeyError`` if the filter is not defined.

        """

        unit = 'W/(m2 um)' if unit is None else unit

        source_fluxd = self._fluxd_state.get()
        fluxd = source_fluxd[filt]

        lambda_eff = source_fluxd.get(filt + '(lambda eff)')
        lambda_pivot = source_fluxd.get(filt + '(lambda pivot)')

        # convert to requested units, may need lambda_pivot
        if lambda_pivot is None:
            equiv = None
        else:
            equiv = u.spectral_density(lambda_pivot)

        try:
            fluxd = fluxd.to(unit, equiv)
        except u.UnitConversionError as e:
            raise type(e)(
                '{}  Is "{}(lambda pivot)" required and'
                ' was it provided?'.format(e, filt))

        return lambda_eff, lambda_pivot, fluxd

    def color_index(self, wfb, unit):
        """Color index (magnitudes) and effective wavelengths.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, or tuple/list of `~synphot.SectralElement`, string
            Two wavelengths, frequencies, or bandpasses.

        unit : string or `~astropy.units.MagUnit`
            Units for the output, e.g., ``astropy.units.ABmag`` or
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
            elif isinstance(wfb[i], str):
                w, pivot, m[i] = self.observe_filter(
                    wfb[i], unit=unit)
                eff_wave.append(w)
            else:
                raise TypeError('Unsupported type for `wfb` type: {}'
                                .format(type(wfb[i])))

        ci = m[0] - m[1]

        return u.Quantity(eff_wave), ci


class solar_spectrum(ScienceState):
    """Get/set the `sbpy` default solar spectrum.

    To retrieve the current default:

    >>> from sbpy.calib import solar_spectrum
    >>> sun = solar_spectrum.get()

    To change it:

    >>> with solar_spectrum.set('E490_2014'):
    ...     # E490_2014 in effect
    ...     pass

    Or, you may use a string:

    >>> with solar_spectrum.set('E490_2014LR'):
    ...     # E490_2014LR in effect
    ...     pass

    """
    _value = 'E490_2014'

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return Sun.from_builtin(value)
        elif isinstance(value, Sun):
            return value
        else:
            raise TypeError(
                "solar_spectrum must be a string or Sun instance.")


class vega_spectrum(ScienceState):
    """Get/set the `sbpy` default Vega spectrum.

    To retrieve the current default:

    >>> from sbpy.calib import vega_spectrum
    >>> vega = vega_spectrum.get()

    To change it:
    >>> with vega_spectrum.set(Vega.from_file(filename))  # doctest: +SKIP
    ...     # Vega from filename in effect
    ...     pass

    """
    _value = 'Bohlin2014'

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return Vega.from_builtin(value)
        elif isinstance(value, Vega):
            return value
        else:
            raise TypeError(
                "vega_spectrum must be a string or Vega instance.")


class solar_fluxd(ScienceState):
    """Get/set the `sbpy` solar flux density.

    To set the current values:

    >>> from sbpy.calib import solar_fluxd
    >>> import sbpy.units as sbu
    >>> with solar_fluxd.set({'V': -26.76 * sbu.VEGAmag}):
    ...     pass

    The units must be flux density per unit wavelength or per unit
    frequency.

    To retrieve the current values:

    >>> import sbpy.units as sbu
    >>> with solar_fluxd.set({'V': -26.76 * sbu.VEGAmag}):
    ...     S = solar_fluxd.get()
    ...     print(S['V'])
    -26.76 mag(VEGA)

    Multiple values are allowed:

    >>> import astropy.units as u
    >>> solar_fluxd.set({
    ...     'PS1_g': -26.54 * u.ABmag,
    ...     'PS1_r': -26.93 * u.ABmag,
    ...     'PS1_i': -27.05 * u.ABmag
    ... })  # doctest: +IGNORE_OUTPUT

    When wavelength is required, specify ``bandpass(lambda eff)`` for
    effective wavelength and/or ``bandpass(lambda pivot)`` for pivot
    wavelength:

    >>> import sbpy.units as sbu
    >>> solar_fluxd.set({
    ...     'V': -26.76 * sbu.VEGAmag,
    ...     'V(lambda eff)': 548 * u.nm,
    ...     'V(lambda pivot)': 551 * u.nm
    ... })  # doctest: +IGNORE_OUTPUT

    """
    _value = {}  # default is disabled


class vega_fluxd(ScienceState):
    """Get/set the `sbpy` Vega flux density.

    To set the current values:

    >>> from sbpy.calib import vega_fluxd
    >>> import astropy.units as u
    >>> with vega_fluxd.set({'V': 3674 * u.Jy}):
    ...     pass

    The units must be flux density per unit wavelength or per unit
    frequency.

    To retrieve the current values:

    >>> import astropy.units as u
    >>> with vega_fluxd.set({'V': 3674 * u.Jy}):
    ...     S = vega_fluxd.get()
    ...     print(S['V'])
    3674.0 Jy

    Multiple values are allowed:

    >>> import astropy.units as u
    >>> vega_fluxd.set({
    ...     'PS1_g': 4026 * u.Jy,
    ...     'PS1_r': 3252 * u.Jy,
    ...     'PS1_i': 2656 * u.Jy
    ... })  # doctest: +IGNORE_OUTPUT

    When wavelength is required, specify ``bandpass(lambda eff)`` for
    effective wavelength and/or ``bandpass(lambda pivot)`` for pivot
    wavelength:

    >>> import astropy.units as u
    >>> vega_fluxd.set({
    ...     'V': 3674 * u.Jy,
    ...     'V(lambda eff)': 548 * u.nm,
    ...     'V(lambda pivot)': 551 * u.nm
    ... })  # doctest: +IGNORE_OUTPUT

    """
    _value = {}  # default is disabled


class Sun(SpectralStandard):
    """Solar spectrum.

    Most functionality requires `synphot`.  The exception is retrieval
    of flux densities by filter name via the ``solar_fluxd`` context
    manager.


    Parameters
    ----------
    wave : `~astropy.units.Quantity`
        Spectral wavelengths.

    fluxd : `~astropy.units.Quantity`
        Solar spectral flux density at 1 au.

    description : string, optional
        Brief description of the source spectrum.

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
    Get the default solar spectrum:

    >>> sun = Sun.from_default()

    Create solar standard from `synphot.SourceSpectrum`:

    >>> import astropy.constants as const
    >>> import synphot
    >>> source = (synphot.SourceSpectrum(synphot.BlackBody1D, temperature=5770)
    ...           (3.14159 * const.R_sun**2 / const.au**2).decompose())
    >>> sun = Sun(source, description='5770 K blackbody Sun')

    Create solar standard from an array:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> import astropy.constants as const
    >>> from astropy.modeling.blackbody import blackbody_lambda
    >>> wave = np.logspace(-1, 2) * u.um
    >>> fluxd = (blackbody_lambda(wave, 5770 * u.K) * np.pi * u.sr
    ...          (const.R_sun**2 / const.au**2).decompose())
    >>> sun = Sun.from_array(wave, fluxd, description='5770 K blackbody Sun')

    Create solar standard from a file:

    >>> sun = Sun.from_file('filename')        # doctest: +SKIP

    Interpolate to 0.62 μm:
    >>> sun(0.62 * u.um)                       # doctest: +FLOAT_CMP
    <Quantity 1720.5108871 W / (m2 um)>

    Observe as through a spectrometer:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> sun = Sun.from_default()
    >>> wave = np.linspace(1, 2.5) * u.um
    >>> fluxd = sun.observe(wave)              # doctest: +IGNORE_OUTPUT

    Observe as through a filter:

    >>> from sbpy.utils import get_bandpass
    >>> sun = Sun.from_default()
    >>> v = get_bandpass('Johnson V')
    >>> sun.observe(v)               # doctest: +FLOAT_CMP
    <Quantity [1839.93273227] W / (m2 um)>

    Observe through a filter, using `sbpy`'s filter calibration system:

    >>> from sbpy.calib import solar_fluxd
    >>> import sbpy.units as sbu
    >>> solar_fluxd.set({
    ...     'V': -26.76 * sbu.VEGAmag,
    ...     'V(lambda eff)': 548 * u.nm,
    ...     'V(lambda pivot)': 551 * u.nm
    ... })
    >>> sun = Sun.from_default()
    >>> print(sun.observe('V'))
    -26.76 VEGAmag

    """

    _sources = solar_sources
    _spectrum_state = solar_spectrum
    _fluxd_state = solar_fluxd


class Vega(SpectralStandard):
    """Vega spectrum.

    Parameters
    ----------
    wave : `~astropy.units.Quantity`
        Spectral wavelengths.

    fluxd : `~astropy.units.Quantity`
        Spectral flux density.

    description : string, optional
        Brief description of the source spectrum.

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
    Get the default Vega spectrum:
    >>> vega = Vega.from_default()               # doctest: +IGNORE_OUTPUT

    Create Vega from a file:
    >>> vega = Vega.from_file('filename')        # doctest: +SKIP

    Evaluate Vega at 1 μm (interpolation):
    >>> print(vega(1 * u.um))                    # doctest: +FLOAT_CMP
    6.326492514857613e-09 W / (m2 um)

    Observe Vega through as if through a spectrometer:
    >>> import numpy as np
    >>> wave = np.linspace(0.4, 0.6) * u.um
    >>> spec = vega.observe(wave)

    Observe Vega through a filter:
    >>> import sbpy.utils
    >>> V = sbpy.utils.get_bandpass('Johnson V')
    >>> fluxd = vega.observe(V)

    User provided calibration:
    >>> from sbpy.calib import vega_fluxd
    >>> vega = Vega(None)
    >>> with vega_fluxd.set({'V': 3674 * u.Jy,
    ...                      'V(lambda pivot)': 5511 * u.AA}):
    ...     print(vega.observe('V', unit='Jy'))
    3674.0 Jy

    """

    _sources = vega_sources
    _spectrum_state = vega_spectrum
    _fluxd_state = vega_fluxd
