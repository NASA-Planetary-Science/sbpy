# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u

from .surface import Surface
from ..units.typing import SpectralFluxDensityQuantity


class LambertianSurface(Surface):
    """Abstract base class for Lambertian surfaces."""

    @u.quantity_input
    def absorption(
        self,
        F_i: SpectralFluxDensityQuantity,
        i: u.physical.angle,
    ) -> u.Quantity:
        # use self._min_zero_cos(i) to ensure cos(>= 90 deg) = 0
        cos_i = self._min_zero_cos(i)
        return F_i * (1 - self.phys["albedo"]) * cos_i

    @u.quantity_input
    def emission(
        self,
        X_e: SpectralFluxDensityQuantity,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        # use self._min_zero_cos(e) to ensure cos(>= 90 deg) = 0
        cos_e = self._min_zero_cos(e)
        return X_e * cos_e / np.pi / u.sr

    @u.quantity_input
    def reflectance(
        self,
        F_i: SpectralFluxDensityQuantity,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
    ) -> u.Quantity:
        # use self._min_zero_cos(theta) to ensure cos(>= 90 deg) = 0
        cos_i = self._min_zero_cos(i)
        cos_e = self._min_zero_cos(e)
        return F_i * self.phys["albedo"] * cos_i * cos_e / np.pi
