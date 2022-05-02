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
    'albedo_unit',
    'dimensionless_albedo',
    'projected_size',
]

from warnings import warn
import numpy as np
import astropy.units as u
import astropy.constants as const
from ..exceptions import SbpyWarning
from ..calib import (
    Vega, Sun, vega_fluxd, FilterLookupError, UndefinedSourceError
)
from .. import data as sbd
from ..spectroscopy.sources import SinglePointSpectrumError, SynphotRequired
from ..exceptions import OptionalPackageUnavailable, SbpyWarning


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

# various reflectance and albedo units
albedo_unit = u.def_unit(['albedo'], represents=u.dimensionless_unscaled,
                    doc=('Integrated reflectance of a planetary body at '
                         'arbitrary phase angle.  Albedo is the product '
                         'of geometric albedo and disk-integrated phase '
                         'function.'))


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
        except SynphotRequired:
            warn(OptionalPackageUnavailable(
                'synphot is required for Vega-based magnitude conversions'
                ' with {}'.format(wfb)))
        except UndefinedSourceError:
            pass
        except u.UnitConversionError as e:
            warn(SbpyWarning(str(e)))

    return equiv


@u.quantity_input(flux=['W / (m2 nm)', 'W / (m2 Hz)', VEGA, u.AB, u.ST],
                  cross_section='km2', albedo=albedo_unit, rh=u.au,
                  delta=u.au)
def dimensionless_albedo(wfb, flux=None, cross_section=None, albedo=None,
                         rh=1*u.au, delta=1*u.au, **kwargs):
    """Dimensionless albedo related equivalencies.

    Supports conversion between disk-integrated albedo, scattering
    cross-section, and the total flux or magnitude.  Uses `sbpy`'s
    photometric calibration system: `~sbpy.calib.solar_spectrum` and
    `~sbpy.calib.solar_fluxd`.

    To perform the conversion, one of the three keywords `flux`,
    `cross_section`, or `albedo` must be provided, depending on the units
    of the from and to quantities.

    Albedo is defined as the ratio of the disk-integrated brightness
    of a solar system planetary body at an arbitrary phase angle to that
    of a perfect Lambert disk of the same radius and at the same distance
    as the body, but illuminated and observed perpendicularly.

    Based on this definition, albedo is essentially the product of
    geometric albedo (pv) and the disk-integrated phase function (Phi(a),
    where a is phase angle):

        albedo = pv * Phi(a)

    Spectral flux density equivalencies for Vega are automatically
    used, if possible.  Dimensionless logarithmic units are also supported
    if the corresponding solar value is set by `~sbpy.calib.solar_fluxd.set`.


    Parameters
    ----------
    wfb : `astropy.units.Quantity`, `synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass corresponding to the flux
        density being converted.

    flux : `astropy.units.Quantity`, `astropy.units.Magnitude`, optional
        Total flux or magnitude of a planetary body.

    cross_section : `astropy.units.Qauntity`, optional
        Total scattering cross-section.  One of `cross_section` or
        `reflectance` is required.

    albedo : `astropy.units.Quantity`, optional
        Albedo as defined earlier

    rh : `astropy.units.Quantity`, optional
        Heliocentric distance of the target.

    delta : `astropy.units.Quantity`, optional
        Observer distance of the target.

    **kwargs
        Keyword arguments for `~Sun.observe()`.


    Returns
    -------
    equiv : list
        List of equivalencies


    Examples
    --------
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from sbpy.units import dimensionless_albedo
    >>> from sbpy.units import VEGAmag, spectral_density_vega, albedo_unit
    >>> from sbpy.calib import solar_fluxd
    >>>
    >>> solar_fluxd.set({'V': -26.77471503 * VEGAmag})
    ...                                             # doctest: +IGNORE_OUTPUT

    >>> # calculate albedo from magnitude
    >>> mag = 3.4 * VEGAmag
    >>> xsec = np.pi * (460 * u.km)**2
    >>> alb = mag.to(albedo_unit, dimensionless_albedo('V', cross_section=xsec))
    >>> print('{0:.4f}'.format(alb))
    0.0900 albedo

    >>> # calculate size from magnitude
    >>> alb = 0.09 * albedo_unit
    >>> xsec = mag.to('km2', dimensionless_albedo('V', albedo=alb))
    >>> radius = np.sqrt(xsec/np.pi)
    >>> print('{0:.2f}'.format(radius))
    460.11 km

    >>> # calculate albedo from size
    >>> radius = 460 * u.km
    >>> xsec = np.pi * radius**2
    >>> alb = xsec.to(albedo_unit, dimensionless_albedo('V', flux=mag))
    >>> print('{0:.4f}'.format(alb))
    0.0900 albedo
    """

    if [flux, cross_section, albedo].count(None) != 2:
        raise ValueError('One and only one of `flux`, `cross_section`, '
            'or `albedo` should be specified.')

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
        xsec = cross_section.to_value('au2')
        for fluxd0 in f_sun:
            if fluxd0.unit in [u.mag, u.dB, u.dex]:
                equiv.append((
                    fluxd0.unit, albedo_unit,
                    lambda mag, mag0=fluxd0.value: u.Quantity(mag - mag0,
                        fluxd0.unit).to_value('', u.logarithmic()) \
                        / xsec * np.pi,
                    lambda ref, mag0=fluxd0.value: u.Quantity(ref * xsec \
                        / np.pi).to_value(fluxd0.unit, u.logarithmic()) + mag0
                    ))
            else:
                equiv.append((
                    fluxd0.unit, albedo_unit,
                    lambda fluxd, fluxd0=fluxd0.value: fluxd / (fluxd0 * xsec \
                        / np.pi),
                    lambda ref, fluxd0=fluxd0.value: ref * fluxd0 * xsec / np.pi
                ))
    elif albedo is not None:
        ref = albedo.value / np.pi
        au2km = const.au.to_value('km')**2
        for fluxd0 in f_sun:
            if fluxd0.unit in [u.mag, u.dB, u.dex]:
                equiv.append((
                    fluxd0.unit, u.km**2,
                    lambda mag, mag0=fluxd0.value: u.Quantity(mag - mag0,
                        fluxd0.unit).to_value('', u.logarithmic()) / ref *
                        au2km,
                    lambda xsec, mag0=fluxd0.value: u.Quantity(ref *
                        xsec / au2km).to_value(fluxd0.unit, u.logarithmic())
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
    elif flux is not None:
        au2km = const.au.to_value('km')**2
        try:
            f_sun = sun.observe(wfb, unit=flux.unit, **kwargs)
        except SinglePointSpectrumError:
            f_sun = sun(wfb, unit=unit)
        except (u.UnitConversionError, FilterLookupError):
            return equiv
        if isinstance(flux, u.Magnitude):
            f_over_fsun = (flux - f_sun).to_value('', u.logarithmic())
        else:
            f_over_fsun = flux / f_sun
        equiv.append((
            albedo_unit, u.km**2,
            lambda alb: f_over_fsun / alb * au2km * np.pi,
            lambda xsec: f_over_fsun / xsec * au2km * np.pi,
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
