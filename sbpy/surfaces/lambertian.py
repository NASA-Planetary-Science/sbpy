# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u

from .surface import Surface, min_zero_cos
from ..units.typing import SpectralFluxDensityQuantity


class LambertianSurface(Surface):
    """Lambertian surface absorption, emission, and reflectance."""

    @staticmethod
    @u.quantity_input
    def absorption(
        F_i: SpectralFluxDensityQuantity,
        epsilon: float,
        i: u.physical.angle,
    ) -> u.Quantity:
        # use min_zero_cos(i) to ensure cos(>= 90 deg) = 0
        cos_i = min_zero_cos(i)
        return F_i * epsilon * cos_i

    @staticmethod
    @u.quantity_input
    def emission(
        X_e: SpectralFluxDensityQuantity,
        epsilon: float,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        # use min_zero_cos(e) to ensure cos(>= 90 deg) = 0
        cos_e = min_zero_cos(e)
        return X_e * epsilon * cos_e / np.pi / u.sr

    @staticmethod
    @u.quantity_input
    def reflectance(
        F_i: SpectralFluxDensityQuantity,
        albedo: float,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        # use min_zero_cos(theta) to ensure cos(>= 90 deg) = 0
        cos_i = min_zero_cos(i)
        cos_e = min_zero_cos(e)
        return F_i * albedo * cos_i * cos_e / np.pi / u.sr
