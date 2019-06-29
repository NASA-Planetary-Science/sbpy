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
from astropy.utils.exceptions import AstropyWarning
from ..calib.vega import Vega
from ..calib.sun import Sun
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

    Uses the default `sbpy` Vega spectrum.

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
    >>> fluxd.value   # doctest: +FLOAT_CMP
    3578.9571538333985


    References
    ----------
    [Joh66] Johnson et al. 1966, Commun. Lunar Planet. Lab. 4, 99

    [BM12] Bessell & Murphy 2012, PASP 124, 140-157

    """

    # warn rather than raise an exception so that code that uses
    # spectral_density_vega when it doesn't need it will still run.
    try:
        import synphot
    except ImportError:
        warn(AstropyWarning('synphot required for Vega-based magnitude'
                            ' conversions.'))
        return []

    vega = Vega.from_default()
    if isinstance(wfb, u.Quantity):
        wav = wfb
        fnu0 = vega(wfb, unit='W/(m2 Hz)')
        flam0 = vega(wfb, unit='W/(m2 um)')
    elif isinstance(wfb, (synphot.SpectralElement, str)):
        fnu0 = vega.filt(wfb, unit='W/(m2 Hz)')[2]
        flam0 = vega.filt(wfb, unit='W/(m2 um)')[2]

    return [
        (fnu0.unit, VEGA, lambda x: x / fnu0.value,
         lambda x: x * fnu0.value),
        (flam0.unit, VEGA, lambda x: x / flam0.value,
         lambda x: x * flam0.value)
    ]


@u.quantity_input(cross_section='km2', reflectance='1/sr',
                  f_sun=['W/(m2 um)', 'W/(m2 Hz)'])
def reflectance(cross_section=None, reflectance=None, wfb=None, M_sun=None,
                f_sun=None):
    """Reflectance related equivalencies.

    Supports conversion from/to reflectance and scattering cross-section
    to/from total flux or magnitude at 1 au for both heliocentric and observer
    distances.

    The default bandpass is Johnson V and the solar spectrum model is
    E490-00a (2014) reference solar spectrum (Table 3), doi:10.1520/E0490.  If
    the user does not specify a wavelength, frequency, or bandpass via
    ``wfb``, then `~sbpy.units.VEGAmag` and spectral flux density
    equivalencies will be provided for the V-band.

    Users wanting conversions with magnitude systems incompatible with these
    units will need to define ``M_sun`` in that system or add equivalencies to
    convert from the magnitude system to physical units or Vega-based
    magnitudes.

    The V-band solar magnitude and flux for the default solar spectrum model
    are:

        M_sun('johnson_v') = -26.7747 VEGAmag
        f_sun('johnson_v') = 1839.93 W / (m2 µm)
        f_sun('johnson_v') = 1.86600e-12 W / (m2 Hz)

    If other wavelength/frequency or bandpass is used, then it has to be passed
    via parameter ``wfb``.

    Parameters
    ----------
    cross_section : `astropy.units.Qauntity`
        Total scattering cross-section
    reflectance : `astropy.units.Quantity`
        Average reflectance
    wfb : `astropy.units.Quantity`, `synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass of the corresponding flux
        density being converted.  See
        :func:`~synphot.SpectralElement.from_filter()` for possible
        bandpass names.  ``wfb`` overrides ``f_sun`` and ``M_sun`` if
        present at the same time.
    f_sun : `astropy.units.Quantity`
        Solar flux in a unit convertible to the unit of the quantity to be
        converted.  If ``wfb`` is not provided, then ``f_sun`` will be used in
        the conversion.  ``f_sun`` overrides ``M_sun`` if both present.
    M_sun : `astropy.units.Quantity`
        Solar magnitude in the same magnitude system as the quantity to be
        converted.  If neither ``wfb`` nor ``f_sun`` is provided, then
        ``M_sun`` will be used in the conversion.  If ``M_sun`` is not
        present either, then the default V-band solar magnitude/flux will be
        used.

    Returns
    -------
    eqv : list
        List of equivalencies

    Examples
    --------
    >>> # Convertion between scattering cross-section and reflectance
    >>> # Note that these examples assumes V-band magnitude and uses the
    >>> # Built-in default V-magnitude of the Sun -26.77 VEGAmag.
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from sbpy.units import reflectance, VEGAmag
    >>> mag = 3.4 * VEGAmag
    >>> cross_sec = np.pi * (460 * u.km)**2
    >>> ref = mag.to('1/sr', reflectance(cross_section=cross_sec))
    >>> print('{0:.4f}'.format(ref))
    0.0287 1 / sr
    >>> mag1 = ref.to(VEGAmag, reflectance(cross_section=cross_sec))
    >>> print('{0:.2f}'.format(mag1))
    3.40 mag(VEGA)

    >>> # Convertion bretween magnitude and scattering cross-section
    >>> ref = 0.0287 / u.sr
    >>> cross_sec = mag.to('km2', reflectance(reflectance=ref))
    >>> radius = np.sqrt(cross_sec/np.pi)
    >>> print('{0:.2f}'.format(radius))
    459.69 km
    >>> mag2 = cross_sec.to(VEGAmag, reflectance(reflectance=ref))
    >>> print('{0:.2f}'.format(mag2))
    3.40 mag(VEGA)
    """
    if wfb is not None:
        sun = Sun.from_default()
        try:
            f_sun_lam = sun.observe(wfb, unit='W/(m2 um)')
            f_sun_nu = sun.observe(wfb, unit='W/(m2 Hz)')
        except SinglePointSpectrumError:
            warn(AstropyWarning('Solar spectrum is interpolated.'))
            f_sun_lam = sun(wfb, unit='W/(m2 um)')
            f_sun_nu = sun(wfb, unit='W/(m2 Hz)')
        f_sun_phys = f_sun_lam.to(VEGA, spectral_density_vega(wfb))
        M_sun = f_sun_phys.to(VEGAmag)
        f_sun_phys = f_sun_phys.value
        f_sun_lam = f_sun_lam.value
        f_sun_nu = f_sun_nu.value
    elif f_sun is not None:
        if f_sun.unit.is_equivalent('W/(m2 um)'):
            f_sun_lam = f_sun.to('W/(m2 um)').value
            f_sun_nu = None
        else:
            f_sun_nu = f_sun.to('W/(m2 Hz)').value
            f_sun_lam = None
        M_sun = None
    elif M_sun is not None:
        if not hasattr(M_sun, 'unit'):
            raise TypeError("Argument 'M_sun' has no 'unit' attribute.  "
                            "You may want to pass in an astropy Quantity instead.")
        if not hasattr(M_sun.unit, 'physical_unit'):
            raise u.UnitsError("Argument 'M_sun' must be in a magnitude "
                               "system based on a physical unit")
        f_sun_lam = None
        f_sun_nu = None
        f_sun_phys = M_sun.to(M_sun.unit.physical_unit).value
    else:
        M_sun = -26.77471503 * VEGAmag  # Solar magnitude in 'Johnson-V'
        f_sun_phys = 5.12726792e+10  # Solar flux in VEGA
        f_sun_lam = 1839.93273227  # Solar flux in 'Johnson-V' in W/(m2 um)
        f_sun_nu = 1.86599755e-12  # Solar flux in 'Johnson-V' in W/(m2 Hz)
    equiv = []
    if cross_section is not None:
        xsec = cross_section.to('au2').value
        if M_sun is not None:
            equiv.append((M_sun.unit.physical_unit,
                          u.sr**-1,
                          lambda flux: flux/(f_sun_phys*xsec),
                          lambda ref: ref*f_sun_phys*xsec))
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
        if M_sun is not None:
            equiv.append((M_sun.unit.physical_unit,
                          u.km**2,
                          lambda flux: flux/(f_sun_phys*ref)*u.au.to('km')**2,
                          lambda xsec: f_sun_phys*ref*xsec*u.km.to('au')**2))
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
