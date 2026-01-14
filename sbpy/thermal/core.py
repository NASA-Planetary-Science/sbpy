# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Thermal Module

created on June 27, 2017
"""

__all__ = ["ThermalClass", "NonRotThermalModel", "FastRotThermalModel"]

import abc
import numpy as np
from numpy.linalg import norm
from scipy.integrate import dblquad
import astropy.units as u
import astropy.constants as const
from astropy.modeling.models import BlackBody
from ..data import Phys, Obs, Ephem, dataclass_input, quantity_to_dataclass


class ThermalClass(abc.ABC):
    """Abstract base class for thermal models.

    This class implements the basic calculation for thermal models,
    such as integration of total flux based on a temperature distribution.
    """

    @u.quantity_input(rh=u.km, R=u.km)
    def __init__(self, rh, R, albedo=0.1, emissivity=1.0, beaming=1.0):
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
        self.rh = rh
        self.R = R
        self.albedo = albedo
        self.emissivity = emissivity
        self.beaming = beaming

    @abc.abstractmethod
    def T(self, lon, lat):
        """Temperature on the surface of an object.

        Needs to be overridden in subclasses.  This function needs to be able
        to return a valid quantity for the full range of lon and lat, i.e.,
        include the night side of an object.

        lon : u.Quantity
            Longitude
        lat : u.Quantity
            Latitude
        """
        pass

    @u.quantity_input(wave_freq=u.m, equivalencies=u.spectral())
    def _int_func(self, lon, lat, m, unit, wave_freq):
        """Integral function for `fluxd`.

        Parameters
        ----------
        lon : float
            Longitude in radiance
        lat : float
            Latitude in fradiance
        m : numpy array of shape (3, 3)
            Transformation matrix to convert a vector in the frame to perform
            integration to body-fixed frame.  This matrix can be calculated
            with private method `_transfer_to_bodyframe`.  The integration
            is performed in a frame where the sub-observer point is defined at
            lon = 0, lat = 0.
        unit : str or astropy.units.Unit
            Unit of the integral function.
        wave_freq : u.Quantity
            Wavelength or frequency of calculation

        Returns
        -------
        float : Integral function to calculate total flux density.
        """
        _, lon1, lat1 = xyz2sph(m.dot(sph2xyz(lon, lat)))
        T = self.T(lon1 * u.rad, lat1 * u.rad)
        if np.isclose(T, 0 * u.K):
            return 0.0
        else:
            # the integral term needs to include a correction for latitude
            # with cos(lat), and a Lambertian emission term cos(lat) + cos(lon)
            coslat = np.cos(lat)
            coslon = np.cos(lon)
            f = BlackBody(T)(wave_freq) * coslat * coslat * coslon
            return f.to_value(unit, u.spectral_density(wave_freq))

    @staticmethod
    @u.quantity_input(sublon=u.deg, sublat=u.deg)
    def _transfer_to_bodyframe(sublon, sublat):
        """Calculate transformation matrix.

        The numerical integration to calculate total flux density is performed
        in a reference frame where the sub-observer point is at
        (lon, lat) = (0, 0).  This matrix supports the transformation from
        this frame to the body-fixed frame to facilitate the calculation of
        surface temperature.
        """
        coslat = np.cos(sublat).value
        if abs(coslat) < np.finfo(type(coslat)).resolution:
            if sublat.value > 0:
                m = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
            else:
                m = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
        else:
            m = twovec(
                [sublon.to_value("deg"), sublat.to_value("deg")], 0, [0, 90], 2
            ).T
        return m

    @u.quantity_input(
        wave_freq=u.m, delta=u.m, lon=u.deg, lat=u.deg, equivalencies=u.spectral()
    )
    def fluxd(
        self,
        wave_freq,
        delta,
        sublon,
        sublat,
        unit="W m-2 um-1",
        error=False,
        epsrel=1e-3,
        **kwargs
    ):
        """Model thermal flux density of an object.

        Parameters
        ----------
        eph : `sbpy.data.Ephem` instance, mandatory
            provide object ephemerides and flux measurements

        Examples
        --------
        >>> from sbpy.thermal import STM
        >>> stmfit = STM.fit(eph) # doctest: +SKIP

        not yet implemented

        """


class NonRotThermalModel(ThermalClass):
    """Non-rotating object temperature distribution, i.e., STM, NEATM"""

    @property
    def Tss(self):
        f_sun = const.L_sun / (4 * np.pi * self.rh**2)
        return (
            (
                (1 - self.albedo)
                * f_sun
                / (self.beaming * self.emissivity * const.sigma_sb)
            )
            ** 0.25
        ).decompose()

    @u.quantity_input(lon=u.deg, lat=u.deg)
    def T(self, lon, lat):
        """Surface temperature at specific (lat, lon)

        lon : u.Quantity in units equivalent to deg
            Longitude
        lat : u.Quantity in units equivalent to deg
            Latitude

        Returns
        -------
        u.Quantity : Surface temperature.
        """
        coslon = np.cos(lon)
        coslat = np.cos(lat)
        prec = np.finfo(coslat.value).resolution
        if (abs(coslon) < prec) or (abs(coslat) < prec) or (coslon < 0):
            return 0 * u.K
        else:
            return self.Tss * (coslon * coslat) ** 0.25


class FastRotThermalModel(ThermalClass):
    """Fast-rotating object temperature distribution, i.e., FRM"""

    @property
    def Tss(self):
        f_sun = const.L_sun / (4 * np.pi * self.rh**2)
        return (
            ((1 - self.albedo) * f_sun / (np.pi * self.emissivity * const.sigma_sb))
            ** 0.25
        ).decompose()

    @u.quantity_input(lon=u.deg, lat=u.deg)
    def T(self, lon, lat):
        """Surface temperature at specific (lat, lon)

        lon : u.Quantity in units equivalent to deg
            Longitude
        lat : u.Quantity in units equivalent to deg
            Latitude

        Returns
        -------
        u.Quantity : Surface temperature.
        """
        coslat = np.cos(lat)
        return self.Tss * coslat**0.25
