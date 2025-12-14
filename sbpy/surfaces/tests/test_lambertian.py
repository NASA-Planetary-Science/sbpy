# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy.coordinates import Angle
from astropy import units as u

from ..lambertian import LambertianSurface


def test_absorption():
    F_i = 1 * u.Jy
    epsilon = 0.9
    i = Angle([0, 30, 45, 60, 90, 100], "deg")
    expected = 0.9 * np.array([1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0]) * u.Jy

    surface = LambertianSurface()
    result = surface.absorption(F_i, epsilon, i=i)
    assert u.allclose(result, expected)


def test_emission():
    I_e = 1 * u.Jy / u.sr
    epsilon = 0.9
    e = Angle([0, 30, 45, 60, 90, 100], "deg")
    expected = (
        0.9 * np.array([1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0]) * u.Jy / u.sr
    )

    surface = LambertianSurface()
    result = surface.emission(I_e, epsilon, e=e)
    assert u.allclose(result, expected)

    phi = np.random.rand(len(e)) * 360 * u.deg  # independent of phi
    result = surface.emission(I_e, epsilon, e=e, phi=phi)
    assert u.allclose(result, expected)


def test_emission_from_vectors():
    I_e = 1 * u.Jy / u.sr
    epsilon = 0.9
    expected = (0.9 * 0.5) * u.Jy / u.sr

    # equivalent to (i, e, phi) = (30, 60, 90) deg
    n = [1, 0, 0]
    r = [0.8660254, 0.5, 0] * u.au
    ro = [0.5, -0.8660254, 0] * u.au

    surface = LambertianSurface()
    result = surface.emission_from_vectors(I_e, epsilon, n=n, r=r, ro=ro)
    assert u.allclose(result, expected)


def test_reflectance():
    F_i = 1 * u.Jy
    albedo = 0.1
    i = Angle([0, 30, 45, 60, 90, 0, 90, 100] * u.deg)
    e = Angle([0, 60, 45, 30, 0, 90, 90, 0] * u.deg)
    phi = np.random.rand(len(i)) * 360 * u.deg  # independent of phi

    surface = LambertianSurface()
    result = surface.reflectance(F_i, albedo, i=i, e=e, phi=phi)
    expected = (
        0.1
        * np.array([1, np.sqrt(3) / 4, 1 / 2, np.sqrt(3) / 4, 0, 0, 0, 0])
        / np.pi
        * u.Jy
        / u.sr
    )
    assert u.allclose(result, expected)
