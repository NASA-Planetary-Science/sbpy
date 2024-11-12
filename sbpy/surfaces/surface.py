# Licensed under a 3-clause BSD style license - see LICENSE.rst

from abc import ABC, abstractmethod

import numpy as np
from astropy import units as u

from ..data.phys import Phys
from ..data.decorators import dataclass_input
from ..units.typing import SpectralFluxDensityQuantity


class Surface(ABC):
    """Abstract base class for all small-body surfaces.


    Parameters
    ----------
    phys : Phys
        Surface physical parameters, e.g., albedo.

    """

    @dataclass_input
    def __init__(self, phys: Phys):
        self.phys = phys

    @staticmethod
    def _min_zero_cos(a: u.physical.angle) -> u.Quantity:
        """Use to ensure that cos(>=90 deg) equals 0."""

        # handle scalars separately
        if a.ndim == 0 and u.isclose(np.abs(a), 90 * u.deg):
            return u.Quantity(0)

        x = np.cos(a)
        x[u.isclose(np.abs(a), 90 * u.deg)] = 0

        return np.maximum(x, 0)

    @abstractmethod
    def absorptance(self, i: u.physical.angle) -> u.Quantity:
        r"""Absorption of incident light.

        The surface is illuminated by incident flux density, :math:`F_i`, at an
        angle of :math:`i`, measured from the surface normal direction.


        Parameters
        ----------
        i : `~astropy.units.Quantity`
            Angle from normal of incident light.


        Returns
        -------
        a : `~astropy.units.Quantity`
            Surface absorptance.

        """

    @abstractmethod
    def emittance(self, e: u.physical.angle, phi: u.physical.angle) -> u.Quantity:
        r"""Emission of light.

        The surface is observed at an angle of :math:`e`, measured from the
        surface normal direction, and at a solar phase angle of :math:`phi`.


        Parameters
        ----------
        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.


        Returns
        -------
        e : `~astropy.units.Quantity`
            Surface emittance.

        """

    @abstractmethod
    def reflectance(
        self, i: u.physical.angle, e: u.physical.angle, phi: u.physical.angle
    ) -> u.Quantity:
        r"""Reflectance.

        The surface is illuminated by incident flux density, :math:`F_i`, at an
        angle of :math:`i`, and emitted toward an angle of :math:`e`, measured
        from the surface normal direction.  :math:`\phi` is the
        source-target-observer (phase) angle.


        Parameters
        ----------
        i : `~astropy.units.Quantity`
            Angle from normal of incident light.

        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.


        Returns
        -------
        r : `~astropy.units.Quantity`
            Surface reflectance.

        """

    @abstractmethod
    def radiance(
        self,
        F_i: SpectralFluxDensityQuantity,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        """Observed radiance from a surface.


        Parameters
        ----------
        F_i : `astropy.units.Quantity`
            Incident light (spectral flux density).

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.

        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.


        Returns
        -------
        radiance : `~astropy.units.Quantity`
            Observed radiance.

        """

    @staticmethod
    def _vectors_to_angles(
        n: np.ndarray[3],
        rs: u.physical.length,
        ro: np.ndarray[3],
    ) -> tuple:
        n_hat = n / np.linalg.norm(n)
        rs_hat = rs / np.linalg.norm(rs)
        ro_hat = ro / np.linalg.norm(ro)

        i = u.Quantity(np.arccos(np.dot(n_hat, rs_hat)), "rad")
        e = u.Quantity(np.arccos(np.dot(n_hat, ro_hat)), "rad")
        phi = u.Quantity(np.arccos(np.dot(rs_hat, ro_hat)), "rad")

        return i, e, phi

    @u.quantity_input
    def radiance_from_vectors(
        self,
        F_i: SpectralFluxDensityQuantity,
        n: np.ndarray[3],
        rs: u.physical.length,
        ro: np.ndarray[3],
    ) -> u.Quantity:
        """Observed radiance from a surface with geometry defined by vectors.

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
        radiance : `~astropy.units.Quantity`
            Observed radiance.

        """

        return self.radiance(F_i, *self._vectors_to_angles(n, rs, ro))
