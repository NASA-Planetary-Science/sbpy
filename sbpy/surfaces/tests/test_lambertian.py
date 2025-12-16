# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy.coordinates import Angle
from astropy import units as u

from ..lambertian import LambertianSurface


def test_absorption():
    epsilon = 0.9
    i = Angle([0, 30, 45, 60, 90, 100], "deg")
    expected = 0.9 * np.array([1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0])

    surface = LambertianSurface()
    result = surface.absorption(epsilon, i)

    assert u.allclose(result, expected)


def test_absorption_from_vectors():
    epsilon = 0.9
    # equivalent to (i, e, phi) = (30, 60, 90) deg
    n = [1, 0, 0]
    r = [0.8660254, 0.5, 0] * u.au
    ro = [0.5, -0.8660254, 0] * u.au

    expected = 0.9 * np.sqrt(3) / 2

    surface = LambertianSurface()
    result = surface.absorption_from_vectors(epsilon, n, r, ro)

    assert u.allclose(result, expected)


def test_emission():
    epsilon = 0.9
    e = Angle([0, 30, 45, 60, 90, 100], "deg")
    expected = 0.9 * np.array([1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0])

    surface = LambertianSurface()
    result = surface.emission(epsilon, e, None)
    assert u.allclose(result, expected)


def test_emission_from_vectors():
    epsilon = 0.9
    # equivalent to (i, e, phi) = (30, 60, 90) deg
    n = [1, 0, 0]
    r = [0.8660254, 0.5, 0] * u.au
    ro = [0.5, -0.8660254, 0] * u.au

    expected = 0.45

    surface = LambertianSurface()
    result = surface.emission_from_vectors(epsilon, n, r, ro)
    assert u.allclose(result, expected)

    # incident vector is optional
    result = surface.emission_from_vectors(epsilon, n, None, ro)
    assert u.allclose(result, expected)


def test_reflectance():
    albedo = 0.1
    i = Angle([0, 30, 45, 60, 90, 0, 90, 100] * u.deg)
    e = Angle([0, 60, 45, 30, 0, 90, 90, 0] * u.deg)
    phi = np.random.rand(len(i)) * 360 * u.deg  # independent of phi
    expected = (
        np.array(
            [
                0.1,
                (0.1 * np.sqrt(3) / 2) / 2,
                (0.1 * np.sqrt(2) / 2) * np.sqrt(2) / 2,
                (0.1 / 2) * np.sqrt(3) / 2,
                0,
                0,
                0,
                0,
            ]
        )
        / np.pi
        / u.sr
    )

    surface = LambertianSurface()
    result = surface.reflectance(albedo, i, e, phi)

    assert u.allclose(result, expected)


def test_reflectance_from_vectors():
    albedo = 0.1
    # equivalent to (i, e, phi) = (30, 60, 90) deg
    n = [1, 0, 0]
    r = [0.8660254, 0.5, 0] * u.au
    ro = [0.5, -0.8660254, 0] * u.au

    expected = (0.1 * np.sqrt(3) / 2) / 2 / np.pi / u.sr

    surface = LambertianSurface()
    result = surface.reflectance_from_vectors(albedo, n, r, ro)

    assert u.allclose(result, expected)
