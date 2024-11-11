# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

import astropy.units as u

from .surface import Surface


class LambertianSurface(Surface):
    """Abstract base class for Lambertian surfaces."""

    @u.quantity_input
    def absorptance(self, i: u.physical.angle) -> u.Quantity:
        cos_i = np.maximum(self._cos(i), np.zeros_like(i.value))
        return (1 - self.phys["albedo"]) * cos_i

    @u.quantity_input
    def emittance(self, e: u.physical.angle, phi: u.physical.angle) -> u.Quantity:
        return np.maximum(self._cos(e), np.zeros_like(e.value))

    @u.quantity_input
    def reflectance(
        self, i: u.physical.angle, e: u.physical.angle, phi: u.physical.angle
    ) -> u.Quantity:
        cos_i = np.maximum(self._cos(i), np.zeros_like(i.value))
        return self.phys["albedo"] * cos_i * self.emittance(e, phi)
