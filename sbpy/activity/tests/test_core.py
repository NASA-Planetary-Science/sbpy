# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..core import *

class TestCircularAperture:
    def test_coma_equivalent_radius(self):
        r = 1 * u.arcsec
        aper = CircularAperture(r)
        assert aper.coma_equivalent_radius() == 1 * u.arcsec

class TestAnnularAperture:
    def test_coma_equivalent_radius(self):
        shape = [1, 2] * u.arcsec
        aper = AnnularAperture(shape)
        assert aper.coma_equivalent_radius() == 1 * u.arcsec

class TestRectangularAperture:
    def test_coma_equivalent_radius(self):
        shape = (0.8, 2) * u.arcsec
        aper = RectangularAperture(shape)
        r = aper.coma_equivalent_radius()
        assert np.isclose(r.value, 0.66776816346357259)

class TestGaussianAperture:
    def test_coma_equivalent_radius(self):
        sigma = 1 * u.arcsec
        aper = GaussianAperture(sigma)
        assert aper.coma_equivalent_radius() == np.sqrt(np.pi / 2) * u.arcsec

