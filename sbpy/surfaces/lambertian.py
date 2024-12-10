# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u

from .surface import Surface


class LambertianSurface(Surface):
    """Abstract base class for Lambertian surfaces."""

    @u.quantity_input
    def absorptance(self, i: u.physical.angle) -> u.Quantity:
        # use self._min_zero_cos(i) which ensures cos(>= 90 deg) = 0
        cos_i = self._min_zero_cos(i)
        return (1 - self.phys["albedo"]) * cos_i

    @u.quantity_input
    def emittance(self, e: u.physical.angle, phi: u.physical.angle) -> u.Quantity:
        # use self._min_zero_cos(e) which ensures cos(>= 90 deg) = 0
        return self._min_zero_cos(e)

    @u.quantity_input
    def reflectance(
        self, i: u.physical.angle, e: u.physical.angle, phi: u.physical.angle
    ) -> u.Quantity:
        # use self._min_zero_cos(e) which ensures cos(>= 90 deg) = 0
        cos_i = self._min_zero_cos(i)
        return self.phys["albedo"] * cos_i * self.emittance(e, phi)
