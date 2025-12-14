# Licensed under a 3-clause BSD style license - see LICENSE.rst

from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from astropy import units as u

from ..data.decorators import dataclass_input
from ..units.typing import SpectralFluxDensityQuantity, SpectralRadianceQuantity


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
        i: u.physical.angle,
    ) -> SpectralFluxDensityQuantity:
        r"""Absorption of directional, incident light.

        The surface is illuminated at an angle of :math:`i`, measured from the
        surface normal direction.


        Parameters
        ----------
        i : `~astropy.units.Quantity`
            Angle from normal of incident light.


        Returns
        -------
        a : `~astropy.units.Quantity`
            Fraction of incident light absorbed.

        """

    @staticmethod
    @abstractmethod
    def emission(
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> SpectralRadianceQuantity:
        r"""Emission of directional light from a surface.

        The surface is observed at an angle of :math:`e`, measured from the
        surface normal direction.  Anisotropic emission is characterized by the
        angle `phi`.


        Parameters
        ----------

        e : `~astropy.units.Quantity`
            Observed angle from normal.

        phi : `~astropy.units.Quantity`
            Angle to account for anisotropic emission.


        Returns
        -------
        F_e : `~astropy.units.Quantity`
            Spectral radiance / specific intensity received by the observer.

        """

    @staticmethod
    @abstractmethod
    def reflectance(
        F_i: SpectralFluxDensityQuantity,
        albedo: float,
        *,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> SpectralRadianceQuantity:
        r"""Bidirectional reflectance.

        The surface is illuminated by incident spectral flux density
        (irradiance), :math:`F_i`, at an angle of :math:`i`, and emitted toward
        an angle of :math:`e`, measured from the surface normal direction.
        :math:`\phi` is the source-target-observer (phase) angle.  Both the
        source and the emitted light are assumed to be collimated.


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
        I_r : `~astropy.units.Quantity`
            Spectral radiance / specific intensity received by the observer.

        """

    @staticmethod
    def _vectors_to_angles(
        n: np.ndarray,
        r: u.physical.length,
        ro: u.physical.length,
    ) -> tuple[u.Quantity, u.Quantity, u.Quantity]:
        n_hat = n / np.linalg.norm(n)
        r_hat = r / np.linalg.norm(r)
        ro_hat = ro / np.linalg.norm(ro)

        i = u.Quantity(np.arccos(np.dot(n_hat, r_hat)), "rad")
        e = u.Quantity(np.arccos(np.dot(n_hat, ro_hat)), "rad")
        phi = u.Quantity(np.arccos(np.dot(r_hat, ro_hat)), "rad")

        return i, e, phi

    @u.quantity_input
    def absorption_from_vectors(
        self,
        F_i: SpectralFluxDensityQuantity,
        epsilon: float,
        *,
        n: np.ndarray,
        r: u.physical.length,
        ro: Optional[u.physical.length] = None,
    ) -> SpectralFluxDensityQuantity:
        """Vector-based alternative to `absorption`.

        Input vectors do not need to be normalized.


        Parameters
        ----------
        F_i : `astropy.units.Quantity`
            Incident light (spectral flux density).

        epsilon : float
            Emissivity of the surface.

        n : `numpy.ndarray`
            Surface normal vector.

        r : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `~astropy.units.Quantity`, optional
            Radial vector from the surface to the observer (ignored).


        Returns
        -------
        F_a : `~astropy.units.Quantity`
            Absorbed spectral flux density.

        """

        i, _, _ = self._vectors_to_angles(n, r, [1, 0, 0])
        return self.absorption(F_i, epsilon, i)

    @u.quantity_input
    def emission_from_vectors(
        self,
        I_e: SpectralRadianceQuantity,
        epsilon: float,
        *,
        n: np.ndarray,
        ro: u.physical.length,
        r: Optional[u.physical.length] = None,
    ) -> SpectralRadianceQuantity:
        r"""Vector-based alternative to `emission`.

        Input vectors do not need to be normalized.


        Parameters
        ----------
        I_e : `astropy.units.Quantity`
            Emitted spectral radiance.

        epsilon : float
            Emissivity of the surface.

        n : `numpy.ndarray`
            Surface normal vector.

        ro : `numpy.ndarray`
            Radial vector from the surface to the observer.

        r : `~astropy.units.Quantity`, optional
            Radial vector from the surface to the light source (ignored).


        Returns
        -------
        F_e : `~astropy.units.Quantity`
            Spectral radiance / specific intensity received by the observer.

        """

        _, e, phi = self._vectors_to_angles(n, None, ro)
        return self.emission(I_e, epsilon, e, phi)

    @u.quantity_input
    def reflectance_from_vectors(
        self,
        F_i: SpectralFluxDensityQuantity,
        albedo: float,
        *,
        n: np.ndarray,
        r: u.physical.length,
        ro: u.physical.length,
    ) -> SpectralRadianceQuantity:
        """Vector-based alternative to `reflectance`.

        Input vectors do not need to be normalized.


        Parameters
        ----------
        F_i : `astropy.units.Quantity`
            Incident light (spectral flux density).

        albedo : float
            Surface albedo.

        n : `numpy.ndarray`
            Surface normal vector.

        r : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `numpy.ndarray`
            Radial vector from the surface to the observer.


        Returns
        -------
        F_r : `~astropy.units.Quantity`
            Spectral flux density / spectral irradiance received by the observer.

        """

        return self.reflectance(F_i, albedo, *self._vectors_to_angles(n, r, ro))
