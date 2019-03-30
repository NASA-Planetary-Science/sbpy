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
    'JMmag'
]

from warnings import warn
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from ..spectroscopy.vega import Vega


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
