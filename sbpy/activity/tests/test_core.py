# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..core import *

def test_rho_as_angle():
    # arctan(100 km, 1 au) * 206264.806 = 0.13787950659645942
    rho = rho_as_angle(100 * u.km, {'delta': 1 * u.au})
    assert np.isclose(rho.to(u.arcsec).value, 0.13787950659645942)

def test_rho_as_distance():
    # 1 au * tan(1") = 725.2709438078363
    rho = rho_as_distance(1 * u.arcsec, {'delta': 1 * u.au})
    assert np.isclose(rho.to(u.km).value, 725.2709438078363)

def test_rho_roundtrip():
    a = 10 * u.arcsec
    eph = {'delta': 1 * u.au}
    b = rho_as_angle(rho_as_distance(a, eph), eph)
    assert np.isclose(a.value, b.to(u.arcsec).value)
    
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

