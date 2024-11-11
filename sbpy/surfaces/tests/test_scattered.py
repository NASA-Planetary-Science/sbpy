# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy import units as u

from ...calib import Sun, solar_spectrum
from ..scattered import LambertianSurfaceScatteredSunlight


def test_scattered_light():

    surface = LambertianSurfaceScatteredSunlight({"albedo": 0.1})

    F_i = 1 * u.Unit("W/(m2 um)")
    n = np.array([1, 0, 0])
    rs = [1, 1, 0] * u.au
    ro = [1, -1, 0] * u.au

    # albedo * F_i / rh**2 * cos(45)**2
    # 0.1 * 1 / 2
    expected = 0.05 * u.W / (u.m**2 * u.um * u.sr)
    result = surface.radiance_from_vectors(F_i, n, rs, ro)
    assert u.isclose(result, expected)


def test_scattered_sunlight():

    surface = LambertianSurfaceScatteredSunlight({"albedo": 0.1})

    # fake an easy solar spectrum for testing
    wave = [0.5, 0.55, 0.6] * u.um
    spec = [0.5, 1.0, 1.5] * u.W / (u.m**2 * u.um)
    with solar_spectrum.set(Sun.from_array(wave, spec)):
        n = np.array([1, 0, 0])
        rs = [1, 1, 0] * u.au
        ro = [1, -1, 0] * u.au

        # albedo * F_i / rh**2 * cos(45)**2
        # 0.1 * 1 / 2 / 2
        expected = 0.025 * u.W / (u.m**2 * u.um * u.sr)
        result = surface.scattered_sunlight_from_vectors(0.55 * u.um, n, rs, ro)
        assert u.isclose(result, expected)

        # again to test branching to Sun.observe
        result = surface.scattered_sunlight_from_vectors(
            (0.549, 0.55, 0.551) * u.um, n, rs, ro
        )
        assert u.allclose(result[1], expected, rtol=0.01)
