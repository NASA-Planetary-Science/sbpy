# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy.coordinates import Angle
from astropy import units as u

from ..surface import Surface


class TestingSurface(Surface):
    def absorptance(self, i):  # pragma: no cover
        pass

    def emittance(self, e, phi):  # pragma: no cover
        pass

    def reflectance(self, i: Angle, e: Angle, phi: Angle) -> u.Quantity:
        return np.cos(i) * np.cos(e) * np.sin(phi)

    def radiance(self, F_i, i, e, phi):
        return F_i * self.reflectance(i, e, phi) / u.sr


def test_min_zero_cos():
    a = Angle([-91, -90, 0, 30, 45, 60, 90, 91], "deg")
    result = Surface._min_zero_cos(a)
    expected = [0, 0, 1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0]
    assert np.allclose(result, expected)

    # test scalars
    for i in range(len(a)):
        assert np.isclose(Surface._min_zero_cos(a[i]), expected[i])


def test_radiance():
    surface = TestingSurface({})
    F_i = 1664 * u.W / (u.m**2 * u.um)
    result = surface.radiance(F_i, 30 * u.deg, 30 * u.deg, 60 * u.deg)
    assert u.isclose(result, F_i * np.sqrt(27) / 8 / u.sr)


def test_vectors_to_angles():
    n = [1, 0, 0]
    rs = [1, 1, 0] * u.au
    ro = [1, -1, 0] * u.au
    i, e, phi = Surface._vectors_to_angles(n, rs, ro)
    assert u.isclose(i, 45 * u.deg)
    assert u.isclose(e, 45 * u.deg)
    assert u.isclose(phi, 90 * u.deg)

    n = [1, 0, 0]
    rs = [1 / np.sqrt(3), 1, 0] * u.au
    ro = [np.sqrt(3), -1, 0] * u.au
    i, e, phi = Surface._vectors_to_angles(n, rs, ro)
    assert u.isclose(i, 60 * u.deg)
    assert u.isclose(e, 30 * u.deg)
    assert u.isclose(phi, 90 * u.deg)

    ro = [1 / np.sqrt(3), -1, 0] * u.au
    i, e, phi = Surface._vectors_to_angles(n, rs, ro)
    assert u.isclose(e, 60 * u.deg)
    assert u.isclose(phi, 120 * u.deg)


def test_radiance_from_vectors():
    F_i = 1664 * u.W / (u.m**2 * u.um)
    surface = TestingSurface({})

    n = [1, 0, 0]
    rs = [1, 1, 0] * u.au
    ro = [1, -1, 0] * u.au
    result = surface.radiance_from_vectors(F_i, n, rs, ro)
    assert u.isclose(result, F_i / 2 / u.sr)

    n = [1, 0, 0]
    rs = [1 / np.sqrt(3), 1, 0] * u.au
    ro = [np.sqrt(3), -1, 0] * u.au
    result = surface.radiance_from_vectors(F_i, n, rs, ro)
    assert u.isclose(result, F_i * np.sqrt(3) / 4 / u.sr)

    ro = [1 / np.sqrt(3), -1, 0] * u.au
    result = surface.radiance_from_vectors(F_i, n, rs, ro)
    assert u.isclose(result, F_i * np.sqrt(3) / 8 / u.sr)
