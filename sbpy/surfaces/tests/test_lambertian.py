# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import numpy as np
from astropy.coordinates import Angle
from astropy import units as u

from ..lambertian import LambertianSurface


class TestingSurface(LambertianSurface):
    def radiance():  # pragma: no cover
        pass


@pytest.fixture
def surface():
    return TestingSurface({"albedo": 0.1})


def test_emittance(surface: TestingSurface):
    e = Angle([0, 30, 45, 60, 90, 100], "deg")
    expected = [1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0]
    phi = np.random.rand(len(e)) * 360 * u.deg  # independent of phi
    result = surface.emittance(e, phi)
    assert u.allclose(result, expected)


def test_absorptance(surface: TestingSurface):
    i = Angle([0, 30, 45, 60, 90, 100], "deg")
    expected = 0.9 * np.array([1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0])
    result = surface.absorptance(i)
    assert u.allclose(result, expected)


def test_reflectance(surface: TestingSurface):
    i = Angle([0, 30, 45, 60, 90, 0, 90, 100] * u.deg)
    e = Angle([0, 60, 45, 30, 0, 90, 90, 0] * u.deg)
    phi = np.random.rand(len(i)) * 360 * u.deg  # independent of phi
    result = surface.reflectance(i, e, phi)
    expected = 0.1 * np.array([1, np.sqrt(3) / 4, 1 / 2, np.sqrt(3) / 4, 0, 0, 0, 0])
    assert u.allclose(result, expected)
