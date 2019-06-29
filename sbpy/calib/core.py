# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = [
    'solar_spectrum',
    'vega_spectrum',
    'solar_fluxd',
    'vega_fluxd'
]

from astropy.utils.state import ScienceState
from .sun import Sun
from .vega import Vega


class solar_spectrum(ScienceState):
    """Get/set the `sbpy` default solar spectrum.

    To retrieve the current default:

    >>> from sbpy.calib import solar_spectrum
    >>> sun = solar_spectrum.get()

    To change it:

    >>> with solar_spectrum.set('E490_2014'):
    ...     # E490_2014 in effect
    ...     pass

    Or, you may use a string:

    >>> with solar_spectrum.set('E490_2014LR'):
    ...     # E490_2014LR in effect
    ...     pass

    """
    _value = 'E490_2014'

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return Sun.from_builtin(value)
        elif isinstance(value, Sun):
            return value
        else:
            raise TypeError(
                "solar_spectrum must be a string or Sun instance.")


class vega_spectrum(ScienceState):
    """Get/set the `sbpy` default Vega spectrum.

    To retrieve the current default:

    >>> from sbpy.calib import vega_spectrum
    >>> vega = vega_spectrum.get()

    To change it:
    >>> with vega_spectrum.set(Vega.from_file(filename))  # doctest: +SKIP
    ...     # Vega from filename in effect
    ...     pass

    """
    _value = 'Bohlin2014'

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return Vega.from_builtin(value)
        elif isinstance(value, Vega):
            return value
        else:
            raise TypeError(
                "vega_spectrum must be a string or Vega instance.")


class solar_fluxd(ScienceState):
    """Get/set the `sbpy` solar flux density.

    To set the current values:

    >>> from sbpy.calib import solar_fluxd
    >>> import sbpy.units as sbu
    >>> with solar_fluxd.set({'V': -26.76 * sbu.VEGAmag}):
    ...     pass

    The units must be flux density per unit wavelength or per unit
    frequency.

    To retrieve the current values:

    >>> import sbpy.units as sbu
    >>> with solar_fluxd.set({'V': -26.76 * sbu.VEGAmag}):
    ...     S = solar_fluxd.get()
    ...     print(S['V'])
    -26.76 mag(VEGA)

    Multiple values are allowed:

    >>> import astropy.units as u
    >>> solar_fluxd.set({
    ...     'PS1_g': -26.54 * u.ABmag,
    ...     'PS1_r': -26.93 * u.ABmag,
    ...     'PS1_i': -27.05 * u.ABmag
    ... })  # doctest: +IGNORE_OUTPUT

    When wavelength is required, specify ``bandpass_lambda_eff`` for
    effective wavelength and/or ``bandpass_lambda_pivot`` for pivot
    wavelength:

    >>> import sbpy.units as sbu
    >>> solar_fluxd.set({
    ...     'V': -26.76 * sbu.VEGAmag,
    ...     'V_lambda_eff': 548 * u.nm,
    ...     'V_lambda_pivot': 551 * u.nm
    ... })  # doctest: +IGNORE_OUTPUT

    """
    _value = {}  # default is disabled


class vega_fluxd(ScienceState):
    """Get/set the `sbpy` Vega flux density.

    To set the current values:

    >>> from sbpy.calib import vega_fluxd
    >>> import astropy.units as u
    >>> with vega_fluxd.set({'V': 3674 * u.Jy}):
    ...     pass

    The units must be flux density per unit wavelength or per unit
    frequency.

    To retrieve the current values:

    >>> import astropy.units as u
    >>> with vega_fluxd.set({'V': 3674 * u.Jy}):
    ...     S = vega_fluxd.get()
    ...     print(S['V'])
    3674.0 Jy

    Multiple values are allowed:

    >>> import astropy.units as u
    >>> vega_fluxd.set({
    ...     'PS1_g': 4026 * u.Jy,
    ...     'PS1_r': 3252 * u.Jy,
    ...     'PS1_i': 2656 * u.Jy
    ... })  # doctest: +IGNORE_OUTPUT

    When wavelength is required, specify ``bandpass_lambda_eff`` for
    effective wavelength and/or ``bandpass_lambda_pivot`` for pivot
    wavelength:

    >>> import astropy.units as u
    >>> vega_fluxd.set({
    ...     'V': 3674 * u.Jy,
    ...     'V_lambda_eff': 548 * u.nm,
    ...     'V_lambda_pivot': 551 * u.nm
    ... })  # doctest: +IGNORE_OUTPUT

    """
    _value = {}  # default is disabled
