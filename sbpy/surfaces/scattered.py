# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

import astropy.units as u

from .surface import Surface
from .lambertian import LambertianSurface
from ..calib import Sun
from ..units.typing import SpectralQuantity, SpectralFluxDensityQuantity, UnitLike


class ScatteredLight(Surface):
    """Abstract base class to observe light scattered by a surface."""

    @u.quantity_input
    def radiance(
        self,
        F_i: SpectralFluxDensityQuantity,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        """Observed light reflected from a surface."""
        return F_i * self.reflectance(i, e, phi) / u.sr

    radiance.__doc__ += Surface.radiance.__doc__[
        Surface.radiance.__doc__.index("\n") :  # noqa: E203
    ]


class ScatteredSunlight(ScatteredLight):
    """Abstract base class to observe sunlight scattered by a surface."""

    @u.quantity_input
    def scattered_sunlight(
        self,
        wave_freq: SpectralQuantity,
        rh: u.physical.length,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
        unit: UnitLike = "W/(m2 sr um)",
    ) -> u.Quantity:
        """Radiance from sunlight scattered by a surface.


        Parameters
        ----------
        wave_freq : `astropy.units.Quantity`
            Wavelength or frequency at which to evaluate the Sun.  Arrays are
            evaluated with `sbpy.calib.core.Sun.observe()`.

        rh : `~astropy.units.Quantity`
            Heliocentric distance.

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.

        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.

        unit : `~astropy.units.Unit`, optional
            Unit of the return value.


        Returns
        -------
        radiance : `~astropy.units.Quantity`
            Observed radiance.

        """

        sun = Sun.from_default()
        flux_density_unit = u.Unit(unit) * u.sr
        if wave_freq.size == 1:
            F_i = sun(wave_freq, unit=flux_density_unit)
        else:
            F_i = sun.observe(wave_freq, unit=flux_density_unit)

        F_i /= rh.to_value("au") ** 2
        return self.radiance(F_i, i, e, phi).to(unit)

    @u.quantity_input
    def scattered_sunlight_from_vectors(
        self,
        wave_freq: SpectralQuantity,
        n: np.ndarray[3],
        rs: u.physical.length,
        ro: u.physical.length,
        unit: UnitLike = "W/(m2 sr um)",
    ) -> u.Quantity:
        """Observed sunlight reflected from a surface.


        Parameters
        ----------
        wave_freq : `astropy.units.Quantity`
            Wavelength or frequency at which to evaluate the Sun.  Arrays are
            evaluated with `sbpy.calib.core.Sun.observe()`.

        n : `numpy.ndarray`
            Surface normal vector.

        rs : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `~astropy.units.Quantity`
            Radial vector from the surface to the observer.

        unit : `~astropy.units.Unit`, optional
            Unit of the return value.


        Returns
        -------
        radiance : `~astropy.units.Quantity`
            Observed radiance.

        """
        rh = np.linalg.norm(rs).to("au")
        i, e, phi = self._vectors_to_angles(n, rs, ro)
        return self.scattered_sunlight(wave_freq, rh, i, e, phi, unit=unit)

    __doc__ += Surface.__doc__[Surface.__doc__.index("\n") + 1 :]  # noqa: E203


class LambertianSurfaceScatteredSunlight(LambertianSurface, ScatteredSunlight):
    """Sunlight scattered from a Lambertian surface.

    The surface is assumed to be illuminated by the Sun.

    """

    __doc__ += Surface.__doc__[Surface.__doc__.index("\n") :]  # noqa: E203
