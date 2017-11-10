# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..core import *

class TestCircularAperture:
    def test_coma_equivalent_radius(self):
        r = 1 * u.arcsec
        aper = CircularAperture(r)
        assert aper.coma_equivalent_radius() == 1 * u.arcsec

class TestRectangularAperture:
    def test_coma_equivalent_radius(self):
        shape = (0.8, 2) * u.arcsec
        aper = RectangularAperture(shape)
        r = aper.coma_equivalent_radius()
        assert np.isclose(r.value, 0.66776816346357259)

