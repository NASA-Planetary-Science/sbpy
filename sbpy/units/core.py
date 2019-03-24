# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy units core module

This package defines Vega-based magnitude systems, and a units
equivalency function for converting to/from flux density.

To use these units in the top-level `astropy.units` namespace::

    >>> import sbpy.units
    >>> sbpy.units.enable()    # doctest: +SKIP

"""

__all__ = [
    'spectral_density_vega',
    'enable',
    'VEGA',
    'VEGAmag',
    'JM',
    'JMmag',
    'magnitude_reflectance'
]

from warnings import warn
import numpy as np
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from ..spectroscopy.vega import Vega


VEGA = u.def_unit(['VEGA', 'VEGAflux'], doc='Spectral flux density of Vega.')

VEGAmag = u.MagUnit(VEGA)
VEGAmag.__doc__ = "Vega-based magnitude: Vega is 0 mag at all wavelengths"

JM = u.def_unit(['JM', 'JMflux'], represents=VEGA * 10**(0.4 * 0.03),
                doc=('Johnson-Morgan magnitude system flux density '
                     'zeropoint (Johnson et al. 1966; Bessell & Murphy '
                     '2012).'))

JMmag = u.MagUnit(JM)
JMmag.__doc__ = ("Johnson-Morgan magnitude system: Vega is 0.03 mag at "
                 "all wavelengths (Johnson et al. 1966; Bessell & Murphy "
                 "2012).")


def enable():
    """Enable `sbpy` units in the top-level `~astropy.units` namespace.

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


def magnitude_reflectance(xsec, wfb=None, M_sun=None):
    """Magnitude - reflectance equivalencies for asteroids

    Parameters
    ----------
    xsec : `~astropy.units.Quantity`
        Cross-sectional area of the object of interest
    wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass of the corresponding flux
        density being converted.  See
        :func:`~synphot.SpectralElement.from_filter()` for possible
        bandpass names.  If provided, then this parameter overrides `M_sun`.
    M_sun : `~astropy.units.Quantity`
        Solar magnitude.  If `wfb` is not provided, then `M_sun` will be used
        in the conversion.  If either `wfb` or `M_sun` is not present, then
        the default V-band solar magnitude in V-band, -26.775 will be used.

    Returns
    -------
    eqv : list
        List of equivalencies

    Examples
    --------
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from sbpy.units import magnitude_reflectance
    >>> m = 3.4 * u.mag
    >>> xsec = np.pi * (460 * u.km)**2
    >>> ref = m.to('1/sr', magnitude_reflectance(xsec))
    >>> print(f'{ref:.4f}')
    0.0287 1 / sr
    >>>
    >>> m1 = ref.to(u.mag, magnitude_reflectance(xsec))
    >>> print(f'{m1:.2f}')
    3.40 mag
    """
    if wfb is not None:
        sun = Sun.from_default()
        if isinstance(wfb, u.Quantity):
            M_sun = sun(wfb).to(VEGAmag, spectral_density_vega(wfb)).value
        elif isinstance(wfb, (synphot.SpectralElement, str)):
            M_sun = sun.filt(wfb)[1].to(VEGAmag, spectral_density_vega(wfb)).value
        else:
            raise ValueError('unrecognized type for `wfb`')
    else:
        if M_sun is None:
            M_sun = -26.775
    return [(u.mag, u.sr**-1, lambda mag: 10**((M_sun-mag)*0.4)/(xsec/u.au**2).decompose().value, lambda ref: M_sun-2.5*np.log10(ref*(xsec/u.au**2).decompose().value))]


def magnitude_radius(ref, wfb=None, M_sun=None):
    """Magnitude - radius equivalencies for asteroids

    Parameters
    ----------
    ref : `~astropy.units.Quantity`
        Reflectance of the object of interest
    wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, string
        Wavelength, frequency, or a bandpass of the corresponding flux
        density being converted.  See
        :func:`~synphot.SpectralElement.from_filter()` for possible
        bandpass names.  If provided, then this parameter overrides `M_sun`.
    M_sun : `~astropy.units.Quantity`
        Solar magnitude.  If `wfb` is not provided, then `M_sun` will be used
        in the conversion.  If either `wfb` or `M_sun` is not present, then
        the default V-band solar magnitude in V-band, -26.775 will be used.

    Returns
    -------
    eqv : list
        List of equivalencies

    Examples
    --------
    >>> from astropy import units as u
    >>> from sbpy.units import magnitude_reflectance
    >>> m = 2.19 * u.mag
    >>> ref = 0.09 / u.sr
    >>> radius = m.to('km', magnitude_reflectance(ref)) # doctest: +SKIP
    >>> print(radius) # doctest: +SKIP
    """
    pass
