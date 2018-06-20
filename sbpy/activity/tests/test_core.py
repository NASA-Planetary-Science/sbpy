# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..core import *


def test_rho_as_angle():
    # arctan(100 km, 1 au) * 206264.806 = 0.13787950659645942
    rho = rho_as_angle(100 * u.km, {'delta': 1 * u.au})
    assert np.isclose(rho.to(u.arcsec).value, 0.13787950659645942)


def test_rho_as_length():
    # 1 au * tan(1") = 725.2709438078363
    rho = rho_as_length(1 * u.arcsec, {'delta': 1 * u.au})
    assert np.isclose(rho.to(u.km).value, 725.2709438078363)


def test_rho_roundtrip():
    a = 10 * u.arcsec
    eph = {'delta': 1 * u.au}
    b = rho_as_angle(rho_as_length(a, eph), eph)
    assert np.isclose(a.value, b.to(u.arcsec).value)


class TestCircularAperture:
    def test_str(self):
        assert str(CircularAperture(1 * u.arcsec)) == 'Circular aperture, radius 1.0 arcsec'

    def test_coma_equivalent_radius(self):
        r = 1 * u.arcsec
        aper = CircularAperture(r)
        assert aper.coma_equivalent_radius() == 1 * u.arcsec

    def test_as_length(self):
        r = 1 * u.arcsec
        aper = CircularAperture(r)
        eph = {'delta': 1 * u.au}
        assert aper.as_length(eph).dim == rho_as_length(r, eph)

    def test_as_angle(self):
        r = 100 * u.km
        aper = CircularAperture(r)
        eph = {'delta': 1 * u.au}
        assert aper.as_angle(eph).dim == rho_as_angle(r, eph)


class TestAnnularAperture:
    def test_str(self):
        assert str(AnnularAperture([1, 2] * u.arcsec)) == 'Annular aperture, radii 1.0–2.0 arcsec'

    def test_coma_equivalent_radius(self):
        shape = [1, 2] * u.arcsec
        aper = AnnularAperture(shape)
        assert aper.coma_equivalent_radius() == 1 * u.arcsec

    def test_as_length(self):
        shape = [1, 2] * u.arcsec
        aper = AnnularAperture(shape)
        eph = {'delta': 1 * u.au}
        assert all(aper.as_length(eph).dim == rho_as_length(shape, eph))

    def test_as_angle(self):
        shape = [100, 200] * u.km
        aper = AnnularAperture(shape)
        eph = {'delta': 1 * u.au}
        assert all(aper.as_angle(eph).dim == rho_as_angle(shape, eph))


class TestRectangularAperture:
    def test_str(self):
        assert str(RectangularAperture([1, 2] * u.arcsec)) == 'Rectangular aperture, dimensions 1.0×2.0 arcsec'

    def test_coma_equivalent_radius(self):
        shape = (0.8, 2) * u.arcsec
        aper = RectangularAperture(shape)
        r = aper.coma_equivalent_radius()
        assert np.isclose(r.value, 0.66776816346357259)

    def test_as_length(self):
        shape = [1, 2] * u.arcsec
        aper = RectangularAperture(shape)
        eph = {'delta': 1 * u.au}
        assert all(aper.as_length(eph).dim == rho_as_length(shape, eph))

    def test_as_angle(self):
        shape = [100, 200] * u.km
        aper = RectangularAperture(shape)
        eph = {'delta': 1 * u.au}
        assert all(aper.as_angle(eph).dim == rho_as_angle(shape, eph))


class TestGaussianAperture:
    def test_str(self):
        assert str(GaussianAperture(1 * u.arcsec)) == 'Gaussian aperture, 1-σ width 1.0 arcsec'

    def test_coma_equivalent_radius(self):
        sigma = 1 * u.arcsec
        aper = GaussianAperture(sigma)
        assert aper.coma_equivalent_radius() == np.sqrt(np.pi / 2) * u.arcsec

    def test_as_length(self):
        sig = 1 * u.arcsec
        aper = GaussianAperture(sig)
        eph = {'delta': 1 * u.au}
        assert aper.as_length(eph).dim == rho_as_length(sig, eph)

    def test_as_angle(self):
        sig = 100 * u.km
        aper = GaussianAperture(sig)
        eph = {'delta': 1 * u.au}
        assert aper.as_angle(eph).dim == rho_as_angle(sig, eph)

    def test_call(self):
        aper = GaussianAperture(1 * u.arcsec)
        assert np.isclose(aper(1 * u.arcsec), 0.6065306597126334)
