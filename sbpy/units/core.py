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
    'reflectance'
]

from warnings import warn
import astropy.units as u
from ..exceptions import SbpyWarning
from ..calib import Vega, Sun, vega_fluxd
from ..spectroscopy.sources import SinglePointSpectrumError


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

    Requires `~synphot`.

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
    eqv : list
        List of equivalencies.


    Examples
    --------
    >>> import astropy.units as u
    >>> from sbpy.units import spectral_density_vega, VEGAmag
    >>> m = 0 * VEGAmag
    >>> fluxd = m.to(u.Jy, spectral_density_vega(5500 * u.AA))
    >>> print(fluxd)   # doctest: +FLOAT_CMP
    3578.9571538333985 Jy

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

    # warn rather than raise an exception so that code that uses
    # spectral_density_vega when it doesn't need it will still run.
    equiv = []
    for unit in ['W/(m2 Hz)', 'W/(m2 um)']:
        try:
            if isinstance(wfb, u.Quantity):
                fluxd = vega(wfb, unit=unit)
            else:
                fluxd = vega.observe(wfb, unit=unit)

            equiv.append([
                fluxd.unit, VEGA,
                lambda x: x / fluxd.value,
                lambda x: x * fluxd.value
            ])
        except SynphotRequired:
            warn(SbpyWarning('synphot is required for Vega-based '
                             'magnitude conversions with : {}'
                             .format(wfb)))

    return equiv


@u.quantity_input(cross_section='km2', reflectance='1/sr',
                  f_sun=['W/(m2 um)', 'W/(m2 Hz)'])
def reflectance(wfb, cross_section=None, reflectance=None):
    """Reflectance related equivalencies.

    Supports conversion from/to reflectance and scattering
    cross-section to/from total flux or magnitude at 1 au for both
    heliocentric and observer distances.  Uses `sbpy`'s photometric
    calibration system: `~sbpy.calib.solar_spectrum` and
    `~sbpy.calib.solar_fluxd`.

    Spectral flux density equivalencies for Vega are automatically
    used, if possible.


    Parameters
    ----------
    wfb : `astropy.units.Quantity`, `synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass corresponding to the flux
        density being converted.  See `~sbpy.calib.solar_fluxd` or
        :func:`~synphot.SpectralElement.from_filter()` for possible
        bandpass names.

    cross_section : `astropy.units.Qauntity`, optional
        Total scattering cross-section.  One of `cross_section` or
        `reflectance` is required.

    reflectance : `astropy.units.Quantity`, optional
        Average reflectance.  One of `cross_section` or `reflectance`
        is required.


    Returns
    -------
    eqv : list
        List of equivalencies


    Examples
    --------
    Convertion between scattering cross-section and reflectance
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from sbpy.units import reflectance, VEGAmag, spectral_density_vega
    >>> from sbpy.calib import solar_fluxd, vega_fluxd
    >>>
    >>> solar_fluxd.set({'V': -26.77 * VEGAmag,
    ...                  'V_lambda_pivot': 5511 * u.AA})
    >>> vega_fluxd.set({'V': 3.5885e-08 * u.Unit('W / (m2 um)'),
    ...                 'V_lambda_pivot': 5511 * u.AA})
    >>>
    >>> mag = 3.4 * VEGAmag
    >>> cross_sec = np.pi * (460 * u.km)**2
    >>> ref = mag.to('1/sr', reflectance('V', cross_section=cross_sec))
    >>> print('{0:.4f}'.format(ref))
    0.0287 1 / sr
    >>> mag1 = ref.to(VEGAmag, reflectance('V', cross_section=cross_sec))
    >>> print('{0:.2f}'.format(mag1))
    3.40 mag(VEGA)

    >>> # Convertion bretween magnitude and scattering cross-section
    >>> ref = 0.0287 / u.sr
    >>> cross_sec = mag.to('km2', reflectance('V', reflectance=ref))
    >>> radius = np.sqrt(cross_sec/np.pi)
    >>> print('{0:.2f}'.format(radius))
    459.69 km
    >>> mag2 = cross_sec.to(VEGAmag, reflectance('V', reflectance=ref))
    >>> print('{0:.2f}'.format(mag2))
    3.40 mag(VEGA)

    """

    sun = Sun.from_default()
    try:
        f_sun_lam = sun.observe(wfb, unit='W/(m2 um)')
    except u.UnitConversionError:
        f_sun_lam = None

    try:
        f_sun_nu = sun.observe(wfb, unit='W/(m2 Hz)')
    except u.UnitConversionError:
        f_sun_nu = None

    vega_equiv = spectral_density_vega(wfb)

    # M_sun = -26.77471503 * VEGAmag  # Solar magnitude in 'Johnson-V'
    # f_sun_phys = 5.12726792e+10  # Solar flux in VEGA
    # f_sun_lam = 1839.93273227  # Solar flux in 'Johnson-V' in W/(m2 um)
    # f_sun_nu = 1.86599755e-12  # Solar flux in 'Johnson-V' in W/(m2 Hz)

    equiv = []
    if cross_section is not None:
        xsec = cross_section.to('au2').value
        if f_sun_lam is not None:
            equiv.append((u.Unit('W/(m2 um)'),
                          u.sr**-1,
                          lambda flux: flux/(f_sun_lam*xsec),
                          lambda ref: ref*f_sun_lam*xsec))
        if f_sun_nu is not None:
            equiv.append((u.Unit('W/(m2 Hz)'),
                          u.sr**-1,
                          lambda flux: flux/(f_sun_nu*xsec),
                          lambda ref: ref*f_sun_nu*xsec))
    elif reflectance is not None:
        ref = reflectance.to('1/sr').value
        if f_sun_lam is not None:
            equiv.append((u.Unit('W/(m2 um)'),
                          u.km**2,
                          lambda flux: flux/(f_sun_lam*ref)*u.au.to('km')**2,
                          lambda xsec: f_sun_lam*ref*xsec*u.km.to('au')**2))
        if f_sun_nu is not None:
            equiv.append((u.Unit('W/(m2 Hz)'),
                          u.km**2,
                          lambda flux: flux/(f_sun_nu*ref)*u.au.to('km')**2,
                          lambda xsec: f_sun_nu*ref*xsec*u.km.to('au')**2))
    return equiv
