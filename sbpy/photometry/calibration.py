# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy photometric calibration

"""

__all__ = [
    'convert_mag'
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


def convert_mag(from_value, to_unit, bandpass=None, wave=None):
    """Convert to/from magnitudes.

    Vega is assumed to have an apparent magnitude of 0.03 in the
    Vega-based magnitude system (Bessell & Murphy 2012, PASP 124,
    140-157).

    One of `from_value` or `to_unit` must have units of magnitude.

    If a convertion between flux density per unit wavelength and per
    unit frequency is required, then one of `bandpass` or `wave` must
    be given.


    Parameters
    ----------
    from_value : `~astropy.units.Quantity`
        Value to convert.

    to_unit : `~astropy.units.Unit` or `~astropy.units.MagUnit`
        Unit or magnitude system to convert to.

    bandpass : string or `~synphot.SpectralElement`, optional
        Name or transmission profile of the bandpass.

    wave : `~astropy.units.Quantity`, optional
        Effective wavelength of the flux density.


    Returns
    -------
    v : `~astropy.units.Quantity` or `~astropy.units.Magnitude`
        The converted value.


    Examples
    --------
    >>> import astropy.units as u
    >>> from sbpy.photometry import convert_mag
    >>> fluxd = convert_mag(0 * u.ABmag, u.Jy)
    >>> fluxd.value   # doctest: +FLOAT_CMP
    3630.78

    """

    try:
        _from = validate_mag(from_value)
    except ValueError:
        _from = from_value

    try:
        _to = validate_mag(to_unit)
    except ValueError:
        _to = to_unit

    kwargs = dict(bandpass=bandpass, wave=wave)
    if isinstance(_from, u.Magnitude) and isinstance(_to, u.MagUnit):
        if _from.unit == u.ABmag:
            fluxd = mag_to_fluxd(_from, u.Jy, **kwargs)
        else:
            fluxd = mag_to_fluxd(_from, 'W/(m2 um)', **kwargs)

        v = fluxd_to_mag(fluxd, _to, **kwargs)
    elif isinstance(_from, u.Magnitude):
        v = mag_to_fluxd(_from, _to, **kwargs)
    else:
        v = fluxd_to_mag(_from, _to, **kwargs)

    return v


def fluxd_to_mag(fluxd, to_unit, bandpass=None, wave=None):
    """Convert flux density into magnitude.

    Vega is assumed to have an apparent magnitude of 0.03 in the
    Vega-based magnitude system(Bessell & Murphy 2012, PASP 124,
    140-157).

    If a convertion between flux density per unit wavelength and per
    unit frequency is required, then one of `bandpass` or `wave` must
    be given.


    Parameters
    ----------
    fluxd: `~astropy.units.Quantity`
        Flux density to convert.

    to_unit: `~astropy.units.MagUnit`
        Convert into this magnitude system: VegaMag, ABmag, or STmag.

    bandpass: string or `~synphot.SpectralElement`, optional
        Name or transmission profile of the bandpass.

    wave: `~astropy.units.Quantity`, optional
        Effective wavelength of the flux density.

    Returns
    -------
    m: `~astropy.units.Quantity`
        Equivalent magnitude.


    See also
    --------
    convert_mag : convert to/from magnitudes


    Examples
    --------
    >>> import astropy.units as u
    >>> from sbpy.photometry.calibration import fluxd_to_mag
    >>> m = fluxd_to_mag(3630.78 * u.Jy, u.ABmag)
    >>> m.value   # doctest: +FLOAT_CMP
    1.6e-7

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
    else:
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
    Vega-based magnitude system(Bessell & Murphy 2012, PASP 124,
    140-157).

    If a convertion between flux density per unit wavelength and per
    unit frequency is required, then one of `bandpass` or `wave` must
    be given.


    Parameters
    ----------
    mag: `~astropy.units.Quantity`
        Magnitude to convert: VegaMag, ABmag, or STmag.

    to_unit: `~astropy.units.MagUnit`
        Convert into this flux density unit.

    bandpass: string or `~synphot.SpectralElement`, optional
        Name or transmission profile of the bandpass.

    wave: `~astropy.units.Quantity`, optional
        Effective wavelength of the flux density.


    Returns
    -------
    fluxd: `~astropy.units.Quantity`
        Equivalent flux density.


    See also
    --------
    convert_mag : convert to/from magnitudes


    Examples
    --------
    >>> import astropy.units as u
    >>> from sbpy.photometry.calibration import mag_to_fluxd
    >>> fluxd = mag_to_fluxd(0 * u.ABmag, u.Jy)
    >>> fluxd.value   # doctest: +FLOAT_CMP
    3630.78

    """

    _mag = validate_mag(mag)

    if _mag.unit == VegaMag:
        vega = Vega.from_default()
        if bandpass is not None:
            wave, zp = vega.filt(bandpass, unit=to_unit)
        elif wave is not None:
            zp = vega(wave)
        else:
            raise ValueError('One of bandpass or wave is required for'
                             ' conversion to VegaMag.')
        zp *= 1.0280  # Vega == 0.03 mag
    else:
        zp = 0 * _mag.unit

    if bandpass is not None:
        args = (u.spectral_density(bandpass.pivot()),)
    else:
        args = (u.spectral_density(wave),)

    try:
        fluxd0 = zp.to(to_unit, *args)
    except u.UnitConversionError as e:
        raise (type(e)(str(e) +
                       '.  Did you provide a bandpass or wavelength?')
               .with_traceback(sys.exc_info()[2]))

    return fluxd0 * 10**(-0.4 * _mag.value)


def validate_mag(m):
    """Verify that provided unit may be used with sbpy.

    Converts synphot's VEGAMAG into sbpy's VegaMag.


    Parameters
    ----------
    m: `~astropy.unit.Magnitude`, `~astropy.unit.Quantity`
        Quantity or unit to test.  May use units of ABmag, STmag,
        sbpy's VegaMag, or synphot's VEGAMAG.


    Returns
    -------
    new_mag: `~astropy.unit.Unit`, `~astropy.unit.Quantity`
        Validated magnitude unit/quantity(depending on input), ready
        for use with sbpy.

    """

    if isinstance(m, (u.Quantity, u.Magnitude)):
        unit = m.unit
    else:
        unit = u.Unit(m)

    new_unit = None
    if unit in (VegaMag, u.ABmag, u.STmag):
        new_unit = unit
    elif synphot:
        if unit == synphot.units.VEGAMAG:
            # use sbpy's VegaMag instead
            new_unit = VegaMag

    if unit == u.m:
        pass  # stop
    if new_unit is None:
        raise ValueError('units are not magnitudes, found ' + str(unit))

    if isinstance(m, (u.Quantity, u.Magnitude)):
        return m.value * new_unit
    else:
        return new_unit
