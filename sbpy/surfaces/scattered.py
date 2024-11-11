# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

import astropy.units as u

from .surface import Surface
from .lambertian import LambertianSurface
from ..calib import Sun
from ..units.typing import SpectralQuantity, SpectralRadianceQuantity, UnitLike


class ScatteredSunlight(Surface):
    """Abstract base class to observe light scattered by a surface illuminated by the Sun."""

    @u.quantity_input
    def radiance(
        self,
        wave_freq: SpectralQuantity,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
        rh: u.physical.length,
        unit: UnitLike = "W/(m2 sr um)",
    ) -> u.Quantity:
        """Observed radiance from a surface.


        Parameters
        ----------
        wave_freq : `astropy.units.Quantity`
            Wavelength or frequency at which to evaluate the Sun.  Arrays are
            evaluated with `sbpy.calib.core.Sun.observe()`.

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.

        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.

        rh : `~astropy.units.Quantity`
            Heliocentric distance.

        unit : `~astropy.units.Unit`, optional
            Unit of the return value.


        Returns
        -------
        radiance : `~astropy.units.Quantity`
            Observed radiance.

        """

        sun = Sun.from_default()
        if wave_freq.size == 1:
            F_i = sun(wave_freq)
        else:
            F_i = sun.observe(wave_freq)

        F_i /= rh.to_value("au") ** 2
        return (F_i * self.reflectance(i, e, phi) / u.sr).to(unit)

    @u.quantity_input
    def radiance_from_vectors(
        self,
        wave_freq: SpectralQuantity,
        n: np.ndarray[3],
        rs: u.physical.length,
        ro: u.physical.length,
        unit: UnitLike = "W/(m2 sr um)",
    ) -> u.Quantity:
        """Observed radiance from a surface.


        Parameters
        ----------
        wave_freq : `astropy.units.Quantity`
            Wavelength or frequency at which to evaluate the Sun.  Arrays are
            evaluated with `sbpy.calib.core.Sun.observe()`.

        n : `numpy.ndarray`
            Surface normal vector.

        rs : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `numpy.ndarray`
            Radial vector from the surface to the observer.

        rh : `~astropy.units.Quantity`
            Heliocentric distance.

        unit : `~astropy.units.Unit`, optional
            Unit of the return value.


        Returns
        -------
        radiance : `~astropy.units.Quantity`
            Observed radiance.

        """
        rh = np.linalg.norm(rs).to("au")
        i, e, phi = self._vectors_to_angles(n, rs, ro)
        return self.radiance(wave_freq, i, e, phi, rh, unit=unit)

    __doc__ += Surface.__doc__[Surface.__doc__.index("\n") + 1 :]


class LambertianSurfaceScatteredSunlight(LambertianSurface, ScatteredSunlight):
    """Sunlight scattered from a Lambertian surface.

    The surface is assumed to be illuminated by the Sun.

    """

    __doc__ += Surface.__doc__[Surface.__doc__.index("\n") :]
