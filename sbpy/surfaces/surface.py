# Licensed under a 3-clause BSD style license - see LICENSE.rst

from abc import ABC, abstractmethod

import numpy as np
from astropy import units as u

from ..data.phys import Phys
from ..data.decorators import dataclass_input
from ..units.typing import SpectralFluxDensityQuantity


def min_zero_cos(a: u.physical.angle) -> u.Quantity:
    """Use to ensure that cos(>=90 deg) equals 0."""

    # handle scalars separately
    if a.ndim == 0 and u.isclose(np.abs(a), 90 * u.deg):
        return u.Quantity(0)

    x = np.cos(a)
    x[u.isclose(np.abs(a), 90 * u.deg)] = 0

    return np.maximum(x, 0)


class Surface(ABC):
    """Abstract base class for all small-body surfaces."""

    @staticmethod
    @abstractmethod
    def absorption(
        F_i: SpectralFluxDensityQuantity,
        epsilon: float,
        i: u.physical.angle,
    ) -> u.Quantity:
        r"""Absorption of directional, incident light.

        The surface is illuminated by incident flux density, :math:`F_i`, at an
        angle of :math:`i`, measured from the surface normal direction.


        Parameters
        ----------

        F_i : `astropy.units.Quantity`
            Incident light (spectral flux density / spectral irradiance).

        epsilon : float
            Emissivity of the surface.

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.


        Returns
        -------
        F_a : `~astropy.units.Quantity`
            Absorbed spectral flux density.

        """

    @staticmethod
    @abstractmethod
    def emission(
        X_e: SpectralFluxDensityQuantity,
        epsilon: float,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        r"""Emission of light from a surface, as seen by a distant observer.

        The surface is observed at an angle of :math:`e`, measured from the
        surface normal direction, and at a solar phase angle of :math:`phi`.


        Parameters
        ----------
        X_e : `astropy.units.Quantity`
            Emitted spectral radiance.

        epsilon : float
            Emissivity of the surface.

        e : `~astropy.units.Quantity`
            Observed angle from normal.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.


        Returns
        -------
        F_e : `~astropy.units.Quantity`
            Spectral flux density / spectral irradiance received by the
            observer.

        """

    @staticmethod
    @abstractmethod
    def reflectance(
        F_i: SpectralFluxDensityQuantity,
        albedo: float,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        r"""Bidirectional reflectance.

        The surface is illuminated by incident flux density (irradiance),
        :math:`F_i`, at an angle of :math:`i`, and emitted toward an angle of
        :math:`e`, measured from the surface normal direction.  :math:`\phi` is
        the source-target-observer (phase) angle.  Both the source and the
        emitted light are assumed to be collimated.


        Parameters
        ----------
        F_i : `astropy.units.Quantity`
            Incident light (spectral flux density).

        albedo : float
            Surface albedo.

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.

        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.


        Returns
        -------
        F_r : `~astropy.units.Quantity`
            Spectral flux density / spectral irradiance received by the observer.

        """

    @staticmethod
    def _vectors_to_angles(
        n: np.ndarray,
        rs: u.physical.length,
        ro: np.ndarray,
    ) -> tuple:
        n_hat = n / np.linalg.norm(n)
        rs_hat = rs / np.linalg.norm(rs)
        ro_hat = ro / np.linalg.norm(ro)

        i = u.Quantity(np.arccos(np.dot(n_hat, rs_hat)), "rad")
        e = u.Quantity(np.arccos(np.dot(n_hat, ro_hat)), "rad")
        phi = u.Quantity(np.arccos(np.dot(rs_hat, ro_hat)), "rad")

        return i, e, phi

    @u.quantity_input
    def reflectance_from_vectors(
        self,
        F_i: SpectralFluxDensityQuantity,
        n: np.ndarray,
        rs: u.physical.length,
        ro: np.ndarray,
    ) -> u.Quantity:
        """Vector based alternative to reflectance().

        Input vectors do not need to be normalized.


        Parameters
        ----------
        F_i : `astropy.units.Quantity`
            Incident light (spectral flux density).

        n : `numpy.ndarray`
            Surface normal vector.

        rs : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `numpy.ndarray`
            Radial vector from the surface to the observer.


        Returns
        -------
        F_r : `~astropy.units.Quantity`
            Spectral flux density / spectral irradiance received by the observer.

        """

        return self.reflectance(F_i, *self._vectors_to_angles(n, rs, ro))
