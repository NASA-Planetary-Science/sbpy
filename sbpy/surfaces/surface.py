# Licensed under a 3-clause BSD style license - see LICENSE.rst

from abc import ABC, abstractmethod
from typing import Union

import numpy as np
from astropy import units as u


def min_zero_cos(a: u.Quantity["angle"]) -> u.Quantity[""]:
    """Use to ensure that cos(>=90 deg) equals 0."""

    # handle scalars
    if a.ndim == 0:
        if a >= 90 * u.deg:
            return u.Quantity(0)
        return np.cos(a)

    # handle arrays
    x = u.Quantity(np.zeros(len(a)))
    i = a < 90 * u.deg
    if np.any(i):
        x[i] = np.cos(a[i])

    return x


class Surface(ABC):
    """Abstract base class for all small-body surfaces."""

    @abstractmethod
    def absorptance(
        self, epsilon: u.Quantity[""], i: u.Quantity["angle"]
    ) -> u.Quantity[""]:
        r"""Absorptance of directional, incident light.

        The surface is illuminated at an angle of :math:`i`, measured from the
        surface normal direction.


        Parameters
        ----------
        epsilon : `~astropy.units.Quantity`
            Surface emissivity.

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.


        Returns
        -------
        a : `~astropy.units.Quantity`
            Absorptance.

        """

    @abstractmethod
    def emittance(
        self,
        epsilon: u.Quantity[""],
        e: u.Quantity["angle"],
        phi: u.Quantity["angle"],
    ) -> u.Quantity[""]:
        r"""Emittance of directional light from a surface.

        The surface is observed at an angle of :math:`e`, measured from the
        surface normal direction.  Anisotropic emittance is characterized by the
        angle `phi`.


        Parameters
        ----------
        epsilon : `~astropy.units.Quantity`
            Surface emissivity.

        e : `~astropy.units.Quantity`
            Observed angle from normal.

        phi : `~astropy.units.Quantity`
            Angle to account for anisotropic emittance.


        Returns
        -------
        em : `~astropy.units.Quantity`
            Emittance.

        """

    @abstractmethod
    def reflectance(
        self,
        albedo: u.Quantity[""],
        i: u.Quantity["angle"],
        e: u.Quantity["angle"],
        phi: u.Quantity["angle"],
    ) -> u.Quantity[u.sr**-1]:
        r"""Bidirectional reflectance.

        The surface is illuminated at an angle of :math:`i`, and observed at an
        angle of :math:`e`, measured from the surface normal direction.
        :math:`\phi` is the source-target-observer (phase) angle.


        Parameters
        ----------
        albedo : `~astropy.units.Quantity`
            Surface albedo.

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.

        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.


        Returns
        -------
        bdr : `~astropy.units.Quantity`
            Bidirectional reflectance.

        """

    @staticmethod
    def _angle(a: u.Quantity, b: u.Quantity) -> u.Quantity["angle"]:
        a_hat = a / np.linalg.norm(a)
        b_hat = b / np.linalg.norm(b)
        return u.Quantity(np.arccos(np.dot(a_hat, b_hat)), "rad")

    # @u.quantity_input
    def absorptance_from_vectors(
        self,
        epsilon: u.Quantity[""],
        n: np.ndarray,
        r: u.Quantity["length"],
        ro: u.Quantity["length"] | None,
    ) -> u.Quantity[""]:
        """Vector-based alternative to `absorptance`.

        Input vectors do not need to be normalized.


        Parameters
        ----------
        epsilon : `~astropy.units.Quantity`
            Surface emissivity.

        n : `numpy.ndarray`
            Surface normal vector.

        r : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `~astropy.units.Quantity` or ``None``
            Radial vector from the surface to the observer.  Parameter is unused
            and may be ``None``.


        Returns
        -------
        a : `~astropy.units.Quantity`
            Absorptance.

        """

        i = self._angle(n, r)
        return self.absorptance(epsilon, i)

    # @u.quantity_input
    def emittance_from_vectors(
        self,
        epsilon: u.Quantity[""],
        n: np.ndarray,
        r: u.Quantity["length"],
        ro: u.Quantity["length"],
    ) -> u.Quantity[""]:
        r"""Vector-based alternative to `emittance`.

        Input vectors do not need to be normalized.


        Parameters
        ----------
        epsilon : `~astropy.units.Quantity`
            Surface emissivity.

        n : `numpy.ndarray`
            Surface normal vector.

        r : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `numpy.ndarray`
            Radial vector from the surface to the observer.


        Returns
        -------
        em : `~astropy.units.Quantity`
            Emittance.

        """

        e = self._angle(n, ro)
        phi = self._angle(r, ro)
        return self.emittance(epsilon, e, phi)

    # @u.quantity_input
    def reflectance_from_vectors(
        self,
        albedo: u.Quantity[""],
        n: np.ndarray,
        r: u.Quantity["length"],
        ro: u.Quantity["length"],
    ) -> u.Quantity[u.sr**-1]:
        """Vector-based alternative to `reflectance`.

        Input vectors do not need to be normalized.


        Parameters
        ----------
        albedo : `~astropy.units.Quantity`
            Surface albedo.

        n : `numpy.ndarray`
            Surface normal vector.

        r : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `numpy.ndarray`
            Radial vector from the surface to the observer.


        Returns
        -------
        bdr : `~astropy.units.Quantity`
            Bidirectional reflectance.

        """

        i = self._angle(n, r)
        e = self._angle(n, ro)
        phi = self._angle(r, ro)
        return self.reflectance(albedo, i, e, phi)
