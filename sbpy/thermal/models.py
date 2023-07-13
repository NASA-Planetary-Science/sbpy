# Thermal models

import numpy as np
import astropy.units as u
from sbpy.bib import cite
from .core import NonRotThermalModel, FastRotThermalModel


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
        """
        Parameters
        ----------
        rh : u.Quantity
            Heliocentric distance
        R : u.Quantity
            Radius of asteroid
        albedo : float, u.Quantity
            Bolometric Bond albedo
        emissivity : float, u.Quantity
            Emissivity of surface
        beaming : float, u.Quantity
            Beaming parameter
        """
        kwargs.setdefault('beaming', 0.756)
        super().__init__(*args, **kwargs)

    @staticmethod
    @u.quantity_input(phase=u.deg, phase_slope=u.mag/u.deg)
    def _phase_corr(phase, phase_slope):
        return u.Magnitude((phase * phase_slope).to_value('mag')).physical

    @u.quantity_input(phase=u.deg, phase_slope=u.mag/u.deg)
    def fluxd(self, wave_freq, delta, phase=0*u.deg,
            phase_slope=0.01*u.mag/u.deg, **kwargs):
        """Calculate total flux density.
        """
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

    def fluxd(self, wave_freq, delta, **kwargs):
        """Calculate total flux density."""
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
        """Initialization

        rh : u.Quantity
            Heliocentric distance
        R : u.Quantity
            Radius of asteroid
        albedo : float, u.Quantity
            Bolometric Bond albedo
        emissivity : float, u.Quantity
            Emissivity of surface
        beaming : float, u.Quantity
            Beaming parameter
        """
        super().__init__(*args, **kwargs)

    @u.quantity_input(phase=u.deg)
    def fluxd(self, wave_freq, delta, phase=0*u.deg, **kwargs):
        """Calculate total flux density.
        """
        sublon = phase
        sublat = 0. * u.deg
        return super().fluxd(wave_freq, delta, sublon, sublat, **kwargs)
