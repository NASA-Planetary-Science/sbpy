# Licensed under a 3-clause BSD style license - see LICENSE.rst

from typing import Optional, Union
import numpy as np
import astropy.units as u

from .surface import Surface, min_zero_cos
from ..units.typing import SpectralFluxDensityQuantity, SpectralRadianceQuantity


class LambertianSurface(Surface):
    """Lambertian surface absorption, emission, and reflectance.

    The surface is illuminated at an angle of :math:`i`, and emitted toward an
    angle of :math:`e`, both measured from the surface normal direction.
    :math:`\phi` is the source-target-observer (phase) angle.  Both the source
    and the emitted light are assumed to be collimated.


    Examples
    --------

    Absorption of light:

    >>> import astropy.units as u
    >>> from sbpy.surfaces import LambertianSurface
    >>>
    >>> surface = LambertianSurface()
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>> epsilon = 0.9
    >>> F_i = 100 * u.W / u.m**2 / u.um  # incident light
    >>>
    >>> surface.absorption(F_i, epsilon, i=i)  # doctest: +FLOAT_CMP
    <Quantity 77.94228634 W / (um m2)>

    Thermal emission:

    >>> I_e = 100 * u.W / u.m**2 / u.um / u.sr
    >>> surface.emission(I_e, epsilon, e=e, phi=phi)  # doctest: +FLOAT_CMP
    <Quantity 45. W / (sr um m2)>

    Bidirectional reflectance:

    >>> albedo = 1 - epsilon
    >>> surface.reflectance(F_i, albedo, i=i, e=e, phi=phi)  # doctest: +FLOAT_CMP
    <Quantity 1.37832224 W / (sr um m2)>

    Using vector-based arguments:

    >>> n = [1, 0, 0]
    >>> r = [0.8660254, 0.5, 0] * u.au
    >>> ro = [0.5, -0.8660254, 0] * u.au
    >>> surface.reflectance_from_vectors(F_i, epsilon, n=n, r=r, ro=ro)  # doctest: +FLOAT_CMP
    <Quantity 0.45 W / (um m2)>

    """

    @staticmethod
    @u.quantity_input
    def absorption(
        F_i: SpectralFluxDensityQuantity,
        epsilon: float,
        *,
        i: u.physical.angle,
        **kwargs,
    ) -> u.Quantity:
        # use min_zero_cos(i) to ensure cos(>= 90 deg) = 0
        cos_i = min_zero_cos(i)
        return F_i * epsilon * cos_i

    @staticmethod
    @u.quantity_input
    def emission(
        I_e: SpectralFluxDensityQuantity,
        epsilon: float,
        *,
        e: u.physical.angle,
        **kwargs,
    ) -> u.Quantity:
        # use min_zero_cos(e) to ensure cos(>= 90 deg) = 0
        cos_e = min_zero_cos(e)
        return I_e * epsilon * cos_e

    @staticmethod
    @u.quantity_input
    def reflectance(
        F_i: SpectralFluxDensityQuantity,
        albedo: float,
        *,
        i: u.physical.angle,
        e: u.physical.angle,
        **kwargs,
    ) -> u.Quantity:
        # use min_zero_cos(theta) to ensure cos(>= 90 deg) = 0
        cos_i = min_zero_cos(i)
        cos_e = min_zero_cos(e)
        return F_i * albedo * cos_i * cos_e / np.pi / u.sr

    @u.quantity_input
    def emission_from_vectors(
        self,
        I_e: SpectralRadianceQuantity,
        epsilon: float,
        *,
        n: np.ndarray,
        ro: np.ndarray,
        **kwargs,
    ) -> SpectralRadianceQuantity:
        _, e, _ = self._vectors_to_angles(n, n, ro)
        return self.emission(I_e, epsilon, e=e)
