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

    Absorption of light:

    >>> import astropy.units as u
    >>> from sbpy.surfaces import LambertianSurface
    >>>
    >>> surface = LambertianSurface()
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>>
    >>> surface.absorption(i)  # doctest: +FLOAT_CMP
    <Quantity 0.8660254>

    Thermal emission:

    >>> surface.emission(e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.5>

    Bidirectional reflectance:

    >>> surface.reflectance(i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.13783222 1 / sr>

    Using vector-based arguments:

    >>> n = [1, 0, 0]
    >>> r = [0.8660254, 0.5, 0] * u.au
    >>> ro = [0.5, -0.8660254, 0] * u.au
    >>> surface.reflectance_from_vectors(n, r, ro)  # doctest: +FLOAT_CMP
    <Quantity 0.13783222 1 / sr>

    """

    @u.quantity_input
    def absorption(
        self,
        i: u.physical.angle,
    ) -> u.Quantity[u.dimensionless_unscaled]:
        # use min_zero_cos(i) to ensure cos(>= 90 deg) = 0
        return min_zero_cos(i)

    @u.quantity_input
    def emission(
        self,
        e: u.physical.angle,
        phi: Union[u.physical.angle, None],
    ) -> u.Quantity[u.dimensionless_unscaled]:
        # use min_zero_cos(e) to ensure cos(>= 90 deg) = 0
        return min_zero_cos(e)

    @u.quantity_input
    def reflectance(
        self,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity[u.sr**-1]:
        return self.absorption(i) * self.emission(e, phi) / np.pi / u.sr

    @u.quantity_input
    def emission_from_vectors(
        self,
        n: np.ndarray,
        r: Union[u.physical.length, None],
        ro: u.physical.length,
    ) -> u.Quantity[u.dimensionless_unscaled]:
        e = self._angle(n, ro)
        return self.emission(e, None)
