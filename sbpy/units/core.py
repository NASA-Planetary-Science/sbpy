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
from ..spectroscopy.vega import Vega
from ..spectroscopy.sun import Sun


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
        fnu0 = vega.filt(wfb, unit='W/(m2 Hz)')[1]
        flam0 = vega.filt(wfb, unit='W/(m2 um)')[1]

    return [
        (fnu0.unit, VEGA, lambda x: x / fnu0.value,
         lambda x: x * fnu0.value),
        (flam0.unit, VEGA, lambda x: x / flam0.value,
         lambda x: x * flam0.value)
    ]


@u.quantity_input(cross_section='km2', reflectance='1/sr')
def reflectance(cross_section=None, reflectance=None, wfb=None, M_sun=None):
    """Reflectance related equivalencies.

    Support conversion between reflectance, scattering cross-section, and total
    magnitude.

    The default magnitude system is ``VEGAmag``, where the apparent magnitude
    of Vega is assumed to be 0 in all and any wavelengths and bands.  If other
    magnitude system is used, then it is implicitly inferred from the solar
    magnitude passed by keyword 'M_sun'.

    Parameters
    ----------
    cross_section : `astropy.units.Qauntity`
        Total scattering cross-section
    reflectance : `astropy.units.Quantity`
        Average reflectance
    wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass of the corresponding flux
        density being converted.  See
        :func:`~synphot.SpectralElement.from_filter()` for possible
        bandpass names.  If provided, then this parameter overrides `M_sun`.
    M_sun : `~astropy.units.Quantity`
        Solar magnitude in the same magnitude system as the quantity to be
        converted.  If `wfb` is not provided, then `M_sun` will be used
        in the conversion.  If neither `wfb` nor `M_sun` is present, then
        the default V-band solar magnitude, -26.775 VEGAmag will be used.

    Returns
    -------
    eqv : list
        List of equivalencies

    Examples
    --------
    >>> # Convertion between scattering cross-section and reflectance
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
        try:
            import synphot
        except ImportError:
            warn(AstropyWarning('synphot required for Vega-based magnitude'
                            ' conversions.'))
            return []
        sun = Sun.from_default()
        if isinstance(wfb, u.Quantity):
            f_sun = sun(wfb, unit='W/(m2 um)')
        elif isinstance(wfb, (synphot.SpectralElement, str)):
            f_sun = sun.filt(wfb, unit='W/(m2 um)')[1]
        else:
            raise ValueError('unrecognized type for `wfb`')
        f_sun_phys = f_sun.to(VEGA, spectral_density_vega(wfb))
        M_sun = f_sun_phys.to(VEGAmag)
        f_sun = f_sun.value
    else:
        if M_sun is None:
            M_sun = -26.77471503 * VEGAmag
            f_sun = 1839.93273227
        else:
            if not hasattr(M_sun.unit, 'physical_unit'):
                raise u.UnitTypeError('the magnitude system must be based on a physical unit')
            f_sun = None # M_sun.to('W/(m2 um)', mag_equiv('johnson_v')).value
        f_sun_phys = M_sun.to(M_sun.unit.physical_unit)
    equiv = []
    if cross_section is not None:
        xsec = cross_section.to('au2').value
        equiv.append((M_sun.unit.physical_unit,
                      u.sr**-1,
                      lambda flux: flux/(f_sun_phys*xsec),
                      lambda ref: ref*f_sun_phys*xsec))
        if f_sun is not None:
            equiv.append((u.Unit('W/(m2 um)'),
                          u.sr**-1,
                          lambda flux: flux/(f_sun*xsec),
                          lambda ref: ref*f_sun*xsec))
    if reflectance is not None:
        ref = reflectance.to('1/sr').value
        equiv.append((M_sun.unit.physical_unit,
                      u.km**2,
                      lambda flux: flux/(f_sun_phys*ref)*u.au.to('km')**2,
                      lambda xsec: f_sun_phys*ref*xsec*u.km.to('au')**2))
        if f_sun is not None:
            equiv.append((u.Unit('W/(m2 um)'),
                          u.km**2,
                          lambda flux: flux/(f_sun*ref)*u.au.to('km')**2,
                          lambda xsec: f_sun*ref*xsec*u.km.to('au')**2))
    return equiv
