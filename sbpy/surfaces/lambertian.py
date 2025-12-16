# Licensed under a 3-clause BSD style license - see LICENSE.rst

from typing import Union
import numpy as np
import astropy.units as u

from .surface import Surface, min_zero_cos
from ..units.typing import SpectralFluxDensityQuantity, SpectralRadianceQuantity


class LambertianSurface(Surface):
    """Lambertian surface absorption, emission, and reflectance.

    The surface is illuminated at an angle of :math:`i`, and emitted toward an
    angle of :math:`e`, both measured from the surface normal direction.
    :math:`\\phi` is the source-target-observer (phase) angle.  Both the source
    and the emitted light are assumed to be collimated.


    Examples
    --------

    Absorption of light for a surface with an albedo of 0.1 and an emissivity of
    0.9:

    >>> import astropy.units as u
    >>> from sbpy.surfaces import LambertianSurface
    >>>
    >>> surface = LambertianSurface()
    >>> albedo = 0.1
    >>> epsilon = 1 - albedo
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>>
    >>> surface.absorption(epsilon, i)  # doctest: +FLOAT_CMP
    <Quantity 0.7794229>

    Emission:

    >>> surface.emission(epsilon, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.45>

    Bidirectional reflectance:

    >>> surface.reflectance(albedo, i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.01378322 1 / sr>

    Using vector-based arguments:

    >>> n = [1, 0, 0]
    >>> r = [0.8660254, 0.5, 0] * u.au
    >>> ro = [0.5, -0.8660254, 0] * u.au
    >>> surface.reflectance_from_vectors(albedo, n, r, ro)  # doctest: +FLOAT_CMP
    <Quantity 0.01378322 1 / sr>

    """

    @u.quantity_input
    def absorption(
        self,
        epsilon: float,
        i: u.physical.angle,
    ) -> u.Quantity[u.dimensionless_unscaled]:
        # use min_zero_cos(i) to ensure cos(>= 90 deg) = 0
        return epsilon * min_zero_cos(i)

    @u.quantity_input
    def emission(
        self,
        epsilon: float,
        e: u.physical.angle,
        phi: Union[u.physical.angle, None],
    ) -> u.Quantity[u.dimensionless_unscaled]:
        # use min_zero_cos(e) to ensure cos(>= 90 deg) = 0
        return epsilon * min_zero_cos(e)

    @u.quantity_input
    def reflectance(
        self,
        albedo: float,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity[u.sr**-1]:
        return albedo * min_zero_cos(i) * min_zero_cos(e) / np.pi / u.sr

    @u.quantity_input
    def emission_from_vectors(
        self,
        epsilon: float,
        n: np.ndarray,
        r: Union[u.physical.length, None],
        ro: u.physical.length,
    ) -> u.Quantity[u.dimensionless_unscaled]:
        e = self._angle(n, ro)
        return self.emission(epsilon, e, None)
