# Thermal models

import numpy as np
import astropy.units as u
from sbpy.bib import cite
from .core import NonRotThermalModel, FastRotThermalModel
from ..data import dataclass_input, quantity_to_dataclass


__all__ = ['STM', 'NEATM', 'FRM']


class STM(NonRotThermalModel):
    """Standard thermal model (STM)

    References
    ----------
    Lebofsky, L.A., et al., 1986.  A refined "standard" thermal model for
        asteroids based on observatons of 1 Ceres and 2 Pallas.  Icarus 68,
        239-251.
    """

    @cite({'method': '1986Icar...68..239L'})
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.beaming = 0.756

    @staticmethod
    @u.quantity_input(phase=u.deg, phase_slope=u.mag/u.deg)
    def _phase_corr(phase, phase_slope):
        return u.Magnitude((phase * phase_slope).to_value('mag')).physical

    @u.quantity_input(phase=u.deg, phase_slope=u.mag/u.deg)
    @dataclass_input(eph=Ephem)
    def fluxd(self, wave_freq, eph, phase_slope=0.01*u.mag/u.deg, **kwargs):
        """Calculate total flux density.

        Parameters
        ----------
        wave_freq : u.Quantity
            Wavelength or frequency of observations
        eph : `~sbpy.data.Ephem`, dict_like, number, or
            `~astropy.units.Quantity`
            If `~sbpy.data.Ephem` or dict_like, ephemerides of the object that
            can include the observer distance and phase angle via keywords
            `delta` and `phase`, respectively.  If observer distance is not
            found, then it will be assumed to be 1 au.  If phase angle is not
            found, then it will be assumed to be 0 deg.  The default unit for
            observer distance is au, and the default unit for phase angle is
            degrees.
        """
        delta = u.Quantity(eph['delta'] if 'delta' in eph else 1., u.au)
        phase = u.Quantity(eph['phase'] if 'phase' in eph else 0., u.deg)
        scl = self._phase_corr(phase, phase_slope)
        sublon = 0. * u.deg
        sublat = 0. * u.deg
        return super().fluxd(wave_freq, delta, sublon, sublat, **kwargs) * scl


class FRM(FastRotThermalModel):
    """Fast rotating model (FRM)

    References
    ----------
    Lebofsky, L.A., Spencer, J.R., 1989.  Radiometry and thermal modeling of
        asteroids.  In: Asteroids II, p. 128-147.
    """

    @cite({'method': '1989aste.conf..128L'})
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @quantity_to_dataclass(delta=(Ephem, 'delta'))
    def fluxd(self, wave_freq, delta, **kwargs):
        """Calculate total flux density.

        Parameters
        ----------
        wave_freq : u.Quantity
            Wavelength or frequency of observations
        delta : `~sbpy.data.Ephem`, dict_like, number, or
            `~astropy.units.Quantity`
            If `~sbpy.data.Ephem` or dict_like, ephemerides of the object that
            can include the observer distance via keywords `delta`.  If
            observer distance is not found, then it will be assumed to be
            1 au.  The default unit for observer distance is au.
        """
        delta = u.Quantity(delta['delta'] if 'delta' in delta else 1., u.au)
        sublon = 0. * u.deg
        sublat = 0. * u.deg
        return super().fluxd(wave_freq, delta, sublon, sublat, **kwargs)


class NEATM(NonRotThermalModel):
    """Near-Earth asteroid thermal model (NEATM)

    References
    ----------
    Harris, A.W., 1998.  A Thermal Model for Near-Earth Asteroids.
        Icarus 131, 291-301
    """

    @cite({'method': '1998Icar..131..291H'})
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @u.quantity_input(phase=u.deg)
    def fluxd(self, wave_freq, eph, **kwargs):
        """Calculate total flux density.

        Parameters
        ----------
        wave_freq : u.Quantity
            Wavelength or frequency of observations
        eph : `~sbpy.data.Ephem`, dict_like, number, or
            `~astropy.units.Quantity`
            If `~sbpy.data.Ephem` or dict_like, ephemerides of the object that
            can include the observer distance and phase angle via keywords
            `delta` and `phase`, respectively.  If observer distance is not
            found, then it will be assumed to be 1 au.  If phase angle is not
            found, then it will be assumed to be 0 deg.  The default unit for
            observer distance is au, and the default unit for phase angle is
            degrees.
        """
        delta = u.Quantity(eph['delta'] if 'delta' in eph else 1., u.au)
        sublon = u.Quantity(eph['phase'] if 'phase' in eph else 0., u.deg)
        sublat = 0. * u.deg
        return super().fluxd(wave_freq, delta, sublon, sublat, **kwargs)
