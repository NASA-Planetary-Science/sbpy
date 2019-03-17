# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy units core module
"""

__all__ = [
    'spectral_density_vega',
    'fluxd_vega',
    'vega_mag',
    'j66_mag'
]

import warnings
from fractions import Fraction
import numpy as np
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from ..spectroscopy.vega import Vega


fluxd_vega = u.def_unit('Vega', doc='Spectral flux density of Vega.')
vega_mag = u.mag(fluxd_vega)
j66_mag = u.mag(fluxd_vega * 10**(0.4 * 0.03))


def spectral_density_vega(wfb):
    """Flux density equivalencies with Vega-based magnitude systems.

    Uses the default `sbpy` Vega spectrum.

    Vega is assumed to have an apparent magnitude of 0 in the VEGAMAG
    system (``vega_mag``), and 0.03 in the Johnson-Morgan system
    (``j66_mag``) [Joh66, BM12]_.


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
    >>> from sbpy.units import spectral_density_vega, vega_mag
    >>> m = 0 * vega_mag
    >>> fluxd = m.to(u.Jy, spectral_density_vega(5500 * u.AA))
    >>> fluxd.value   # doctest: +FLOAT_CMP
    3578.9571538333985


    References
    ----------
    [Joh66] Johnson et al. 1966, Commun. Lunar Planet. Lab. 4, 99

    [BM12] Bessell & Murphy 2012, PASP 124, 140-157

    """

    try:
        import synphot
        from synphot.units import VEGAMAG
    except ImportError:
        raise ImportError('synphot required for Vega-based magnitude'
                          ' conversions.')

    vega = Vega.from_default()
    if isinstance(wfb, u.Quantity):
        wav = wfb
        fnu0 = vega(wfb, unit='W/(m2 Hz)')
        flam0 = vega(wfb, unit='W/(m2 um)')
    elif isinstance(wfb, (synphot.SpectralElement, str)):
        fnu0 = vega.filt(wfb, unit='W/(m2 Hz)')[1]
        flam0 = vega.filt(wfb, unit='W/(m2 um)')[1]

    def vegamag_forward(fluxd_vega):
        """nan/inf are returned as -99 mag.

        nan/inf handling consistent with
        :func:`~synphot.units.spectral_density_vega`.

        """
        def f(fluxd):
            m = -2.5 * np.log10(fluxd / fluxd_vega.value)
            m = np.ma.MaskedArray(m, mask=~np.isfinite(m)).filled(-99)
            if np.ndim(m) == 0:
                m = float(m)
            return m
        return f

    def vegamag_backward(fluxd_vega):
        def b(m):
            return fluxd_vega.value * 10**(-0.4 * m)
        return b

    return [
        (fnu0.unit, fluxd_vega, lambda x: x / fnu0.value,
         lambda x: x * fnu0.value),
        (flam0.unit, fluxd_vega, lambda x: x / flam0.value,
         lambda x: x * flam0.value)
    ]
