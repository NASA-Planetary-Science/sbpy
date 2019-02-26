# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy photometric calibration

"""

__all__ = [
    'fluxd_to_mag',
    'mag_to_fluxd'
]

import sys
import numpy as np
import astropy.units as u

from ..units import VegaMag
from ..spectroscopy.vega import Vega
from ..spectroscopy.sun import Sun

try:
    import synphot
except ImportError:
    synphot = None


def fluxd_to_mag(fluxd, to_unit, bandpass=None, wave=None):
    """Convert flux density into magnitude.

    Vega is assumed to have an apparent magnitude of 0.03 in the
    Vega-based magnitude system (Bessell & Murphy 2012, PASP 124,
    140-157).


    Parameters
    ----------
    fluxd : `~astropy.units.Quantity`
        Flux density to convert.

    to_unit : `~astropy.units.MagUnit`
        Convert into this magnitude system: VegaMag, ABmag, or STmag.
        One of ``bandpass`` or ``wave`` is required for ``VegaMag`` or
        when conversion requires a transformation between flux density
        per unit frequency and per unit wavelength.

    bandpass : string or `~synphot.SpectralElement`, optional
        Name or transmission profile of the bandpass.

    wave : `~astropy.units.Quantity`, optional
        Effective wavelength of the flux density.


    Returns
    -------
    m : `~astropy.units.Quantity`
        Equivalent magnitude.

    """

    _to_unit = validate_mag(to_unit)

    if _to_unit == VegaMag:
        vega = Vega.from_default()
        if bandpass is not None:
            zp = vega.filt(bandpass, unit=fluxd.unit)[1]
        elif wave is not None:
            zp = vega(wave)
        else:
            raise ValueError('One of bandpass or wave is required for'
                             ' conversion to VegaMag.')
        # Vega == +0.03 mag: zeropoint is brighter than Vega
        zp *= 1.0280
    elif _to_unit in (u.ABmag, u.STmag):
        zp = 0 * _to_unit

    if bandpass is not None:
        args = (u.spectral_density(bandpass.pivot()),)
    else:
        args = (u.spectral_density(wave),)

    try:
        fluxd0 = zp.to(fluxd.unit, *args)
    except u.UnitConversionError as e:
        raise (type(e)(str(e) +
                       '.  Did you provide a bandpass or wavelength?')
               .with_traceback(sys.exc_info()[2]))

    return -2.5 * np.log10(fluxd / fluxd0).value * to_unit


def mag_to_fluxd(mag, to_unit, bandpass=None, wave=None):
    """Convert magnitude into flux density.

    Vega is assumed to have an apparent magnitude of 0.03 in the
    Vega-based magnitude system (Bessell & Murphy 2012, PASP 124,
    140-157).


    Parameters
    ----------
    mag : `~astropy.units.Quantity`
        Magnitude to convert: VegaMag, ABmag, or STmag.  One of
        ``bandpass`` or ``wave`` is required for ``VegaMag`` or when
        conversion requires a transformation between flux density per
        unit frequency and per unit wavelength.

    to_unit : `~astropy.units.MagUnit`
        Convert into this flux density unit.

    bandpass : string or `~synphot.SpectralElement`, optional
        Name or transmission profile of the bandpass.

    wave : `~astropy.units.Quantity`, optional
        Effective wavelength of the flux density.


    Returns
    -------
    fluxd : `~astropy.units.Quantity`
        Equivalent flux density.

    """

    _mag = validate_mag(mag)

    if _mag.unit == VegaMag:
        vega = Vega.from_default()
        if bandpass is not None:
            wave, zp = Vega.filt(bandpass, unit=fluxd.unit)
        elif wave is not None:
            zp = Vega(wave)
        else:
            raise ValueError('One of bandpass or wave is required for'
                             ' conversion from VegaMag.')
        zp *= 1.0280  # Vega == 0.03 mag
    elif to_unit in (u.ABmag, u.STmag):
        zp = 0 * _mag.unit
    else:
        raise ValueError('Unsupported unit: {}'.format(unit))

    fluxd0 = zp.to(to_unit, u.spectral_density(wave))

    return fluxd0 * 10**(-0.4 * _mag.value)


def validate_mag(m):
    """Verify that provided unit may be used with sbpy.

    Converts synphot's VEGAMAG into sbpy's VegaMag.


    Parameters
    ----------
    m : `~astropy.unit.Unit`, `~astropy.unit.Quantity`
        Quantity or unit to test.  May use units of ABmag, STmag,
        sbpy's VegaMag, or synphot's VEGAMAG.


    Returns
    -------
    new_mag : `~astropy.unit.Unit`, `~astropy.unit.Quantity`
        Validated magnitude unit/quantity (depending on input), ready
        for use with sbpy.

    """

    unit = m.unit if isinstance(m, (u.Quantity, u.Magnitude)) else m

    new_unit = None
    if unit.is_equivalent((VegaMag, u.ABmag, u.STmag)):
        new_unit = unit
    elif synphot:
        if unit == synphot.units.VEGAMAG:
            # use sbpy's VegaMag instead
            new_unit = VegaMag

    if new_unit is None:
        raise ValueError('units are not magnitudes, found ' + str(unit))

    if isinstance(m, (u.Quantity, u.Magnitude)):
        return m.value * new_unit
    else:
        return new_unit
