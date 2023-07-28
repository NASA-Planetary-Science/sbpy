# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy units core module

This package defines Vega-based magnitude systems, and a units
equivalency function for converting to/from flux density.

To use these units in the top-level `astropy.units` namespace::

    >>> import astropy.units as u
    >>> import sbpy.units
    >>> sbpy.units.enable()    # doctest: +IGNORE_OUTPUT
    >>> u.Unit('mag(VEGA)')
    Unit("mag(VEGA)")

"""

__all__ = [
    'hundred_nm',
    'spectral_density_vega',
    'enable',
    'VEGA',
    'VEGAmag',
    'JM',
    'JMmag',
    'reflectance',
    'projected_size',
]

__doctest_requires__ = {
    "spectral_density_vega": ["synphot"]
}

from warnings import warn
import numpy as np
import astropy.units as u
import astropy.constants as const
from ..calib import (
    Vega, Sun, vega_fluxd, FilterLookupError, UndefinedSourceError
)
from .. import data as sbd
from ..spectroscopy.sources import SinglePointSpectrumError
from ..exceptions import RequiredPackageUnavailable, OptionalPackageUnavailable, SbpyWarning


VEGA = u.def_unit(['VEGA', 'VEGAflux'],
                  doc='Spectral flux density of Vega.')

VEGAmag = u.MagUnit(VEGA)
VEGAmag.__doc__ = "Vega-based magnitude: Vega is 0 mag at all wavelengths"

JM = u.def_unit(['JM', 'JMflux'], represents=VEGA * 10**(0.4 * 0.03),
                doc=('Johnson-Morgan magnitude system flux density '
                     'zeropoint (Johnson et al 1966; Bessell & Murphy '
                     '2012).'))

JMmag = u.MagUnit(JM)
JMmag.__doc__ = ("Johnson-Morgan magnitude system: Vega is 0.03 mag at "
                 "all wavelengths (Johnson et al 1966; Bessell & Murphy "
                 "2012).")

hundred_nm = u.def_unit('100 nm', represents=100 * u.nm,
                        doc=('Convenience unit for expressing spectral '
                             'gradients.'))


def enable():
    """Enable `sbpy` units in the top-level `astropy.units` namespace.

    Allows the use in `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    May be used with the ``with`` statement to enable them temporarily.

    """
    import inspect
    return u.add_enabled_units(inspect.getmodule(enable))


def spectral_density_vega(wfb):
    """Flux density equivalencies with Vega-based magnitude systems.

    Uses the default `sbpy` Vega spectrum, or `~sbpy.calib.vega_fluxd`.

    Vega is assumed to have an apparent magnitude of 0 in the
    ``VEGAmag`` system, and 0.03 in the Johnson-Morgan, ``JMmag``,
    system [Joh66, BM12]_.


    Parameters
    ----------
    wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass of the corresponding flux
        density being converted.  See
        :func:`~synphot.SpectralElement.from_filter()` for possible
        bandpass names.


    Returns
    -------
    equiv : list
        List of equivalencies.


    Examples
    --------
    Monochromatic flux density:
    >>> import astropy.units as u
    >>> from sbpy.units import spectral_density_vega, VEGAmag
    >>> m = 0 * VEGAmag
    >>> fluxd = m.to(u.Jy, spectral_density_vega(5557.5 * u.AA))
    >>> print(fluxd)   # doctest: +FLOAT_CMP
    3544.75836649349 Jy

    Use your favorite bandpass and zeropoint (Willmer 2018):
    >>> from sbpy.calib import vega_fluxd
    >>> cal = {'SDSS_r': 3255 * u.Jy, 'SDSS_r_lambda_pivot': 0.6176 * u.um}
    >>> with vega_fluxd.set(cal):
    ...     fluxd = m.to(u.Jy, spectral_density_vega('SDSS_r'))
    ...     print(fluxd)    # doctest: +FLOAT_CMP
    3255.0 Jy

    Use your favorite bandpass and zeropoint (Willmer 2018):
    >>> from sbpy.calib import vega_fluxd
    >>> cal = {'SDSS_r': 3255 * u.Jy, 'SDSS_r_lambda_pivot': 0.6176 * u.um}
    >>> with vega_fluxd.set(cal):
    ...     fluxd = m.to(u.Jy, spectral_density_vega('SDSS_r'))
    ...     print(fluxd)    # doctest: +FLOAT_CMP
    3255.0 Jy

    References
    ----------
    [Joh66] Johnson et al. 1966, Commun. Lunar Planet. Lab. 4, 99

    [BM12] Bessell & Murphy 2012, PASP 124, 140-157

    """

    vega = Vega.from_default()

    # warn rather than raise exceptions so that code that uses
    # spectral_density_vega when it doesn't need it will still run.
    equiv = []
    for unit in ['W/(m2 Hz)', 'W/(m2 um)']:
        try:
            try:
                fluxd0 = vega.observe(wfb, unit=unit)
            except SinglePointSpectrumError:
                fluxd0 = vega(wfb, unit=unit)

            # pass fluxd0 as an optional argument to dereference it,
            # otherwise both equivalencies will use the fluxd0 for
            # W/(m2 um)
            equiv.append((
                fluxd0.unit, VEGA,
                lambda f_phys, fluxd0=fluxd0.value: f_phys / fluxd0,
                lambda f_vega, fluxd0=fluxd0.value: f_vega * fluxd0
            ))
        except RequiredPackageUnavailable:
            warn(OptionalPackageUnavailable(
                'synphot is required for Vega-based magnitude conversions'
                ' with {}'.format(wfb)))
        except UndefinedSourceError:
            pass
        except u.UnitConversionError as e:
            warn(SbpyWarning(str(e)))

    return equiv


@u.quantity_input(cross_section='km2', reflectance='1/sr')
def reflectance(wfb, cross_section=None, reflectance=None, **kwargs):
    """Reflectance related equivalencies.

    Supports conversion from/to reflectance and scattering
    cross-section to/from total flux or magnitude at 1 au for both
    heliocentric and observer distances.  Uses `sbpy`'s photometric
    calibration system: `~sbpy.calib.solar_spectrum` and
    `~sbpy.calib.solar_fluxd`.

    Spectral flux density equivalencies for Vega are automatically
    used, if possible.  Dimensionless logarithmic units are also supported
    if the corresponding solar value is set by `~sbpy.calib.solar_fluxd.set`.


    Parameters
    ----------
    wfb : `astropy.units.Quantity`, `synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass corresponding to the flux
        density being converted.

    cross_section : `astropy.units.Qauntity`, optional
        Total scattering cross-section.  One of `cross_section` or
        `reflectance` is required.

    reflectance : `astropy.units.Quantity`, optional
        Average reflectance.  One of `cross_section` or `reflectance`
        is required.

    **kwargs
        Keyword arguments for `~Sun.observe()`.


    Returns
    -------
    equiv : list
        List of equivalencies


    Examples
    --------
    Convertion between scattering cross-section and reflectance
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from sbpy.units import reflectance, VEGAmag, spectral_density_vega
    >>> from sbpy.calib import solar_fluxd, vega_fluxd
    >>>
    >>> solar_fluxd.set({'V': -26.77471503 * VEGAmag})
    ...                                             # doctest: +IGNORE_OUTPUT
    >>> vega_fluxd.set({'V': 3.5885e-08 * u.Unit('W / (m2 um)')})
    ...                                             # doctest: +IGNORE_OUTPUT
    >>> mag = 3.4 * VEGAmag
    >>> cross_sec = np.pi * (460 * u.km)**2
    >>> ref = mag.to('1/sr', reflectance('V', cross_section=cross_sec))
    >>> print('{0:.4f}'.format(ref))
    0.0287 1 / sr
    >>> mag1 = ref.to(VEGAmag, reflectance('V', cross_section=cross_sec))
    >>> print('{0:.2f}'.format(mag1))
    3.40 mag(VEGA)

    >>> # Convertion between magnitude and scattering cross-section
    >>> ref = 0.0287 / u.sr
    >>> cross_sec = mag.to('km2', reflectance('V', reflectance=ref))
    >>> radius = np.sqrt(cross_sec/np.pi)
    >>> print('{0:.2f}'.format(radius))
    459.69 km
    >>> mag2 = cross_sec.to(VEGAmag, reflectance('V', reflectance=ref))
    >>> print('{0:.2f}'.format(mag2))
    3.40 mag(VEGA)

    """

    # Solar flux density at 1 au in different units
    f_sun = []
    sun = Sun.from_default()
    for unit in ('W/(m2 um)', 'W/(m2 Hz)', VEGA):
        try:
            f_sun.append(sun.observe(wfb, unit=unit, **kwargs))
        except SinglePointSpectrumError:
            f_sun.append(sun(wfb, unit=unit))
        except (u.UnitConversionError, FilterLookupError):
            pass
    if len(f_sun) == 0:
        try:
            f_sun.append(sun.observe(wfb, **kwargs))
        except (SinglePointSpectrumError, u.UnitConversionError, FilterLookupError):
            pass

    # pass fluxd0 as an optional argument to dereference it,
    # otherwise both equivalencies will use the fluxd0 for
    # the last item in f_sun
    equiv = []
    if cross_section is not None:
        xsec = cross_section.to('au2').value
        for fluxd0 in f_sun:
            if fluxd0.unit in [u.mag, u.dB, u.dex]:
                equiv.append((
                    fluxd0.unit, u.sr**-1,
                    lambda mag, mag0=fluxd0.value: u.Quantity(mag - mag0,
                        fluxd0.unit).to('', u.logarithmic()).value / xsec,
                    lambda ref, mag0=fluxd0.value: u.Quantity(ref * xsec).to(
                        fluxd0.unit, u.logarithmic()).value + mag0
                    ))
            else:
                equiv.append((
                    fluxd0.unit, u.sr**-1,
                    lambda fluxd, fluxd0=fluxd0.value: fluxd / (fluxd0 * xsec),
                    lambda ref, fluxd0=fluxd0.value: ref * fluxd0 * xsec
                ))
    elif reflectance is not None:
        ref = reflectance.to('1/sr').value
        au2km = (const.au.to('km')**2).value
        for fluxd0 in f_sun:
            if fluxd0.unit in [u.mag, u.dB, u.dex]:
                equiv.append((
                    fluxd0.unit, u.km**2,
                    lambda mag, mag0=fluxd0.value: u.Quantity(mag - mag0,
                        fluxd0.unit).to('', u.logarithmic()).value / ref *
                        au2km,
                    lambda xsec, mag0=fluxd0.value: u.Quantity(ref *
                        xsec / au2km).to(fluxd0.unit, u.logarithmic()).value
                        + mag0
                    ))
            else:
                equiv.append((
                    fluxd0.unit, u.km**2,
                    lambda fluxd, fluxd0=fluxd0.value: (
                        fluxd / (fluxd0 * ref) * au2km),
                    lambda xsec, fluxd0=fluxd0.value: (
                        fluxd0 * ref * xsec / au2km)
                ))
    return equiv


@sbd.dataclass_input
@sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
def projected_size(eph: sbd.Ephem):
    """Angular size and projected linear size equivalency.

    Based on the tangent formula:

        length = delta * tan(angle)


    Parameters
    ----------
    eph : `~astropy.units.Quantity`, dict-like, `sbpy.data.Ephem`
        Distance to the object as a quantity or the field `'delta'`.


    Returns
    -------
    eqiv : list
        List of equivalencies


    Examples
    --------
    >>> import astropy.units as u
    >>> import sbpy.units as sbu
    >>> (1 * u.arcsec).to('km', sbu.projected_size(1 * u.au))
    ... # doctest: +FLOAT_CMP
    <Quantity [725.27094381] km>
    >>> (725.27 * u.km).to('arcsec', sbu.projected_size(1 * u.au))
    ... # doctest: +FLOAT_CMP
    <Quantity [1.00] arcsec>

    """

    delta = eph['delta'].to('m').value

    equiv = [(
        u.rad, u.m,
        lambda angle, delta=delta: delta * np.tan(angle),
        lambda length, delta=delta: np.arctan(length / delta)
    )]

    return equiv
