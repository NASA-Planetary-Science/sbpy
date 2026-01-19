# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy import pi
from astropy import units as u

from ..sphere import Sphere
from ...calib import Sun, solar_spectrum
from ...surfaces.surface import Surface, min_zero_cos


class SimpleSurface(Surface):
    def absorptance(self, epsilon, i):
        return epsilon * min_zero_cos(i)

    def emittance(self, albedo, e, phi):
        pass  # pragma: no cover

    def reflectance(self, albedo, i, e, phi):
        return albedo * min_zero_cos(i) * min_zero_cos(e)


class TestSphere:
    def test_integrate_i_e_phi(self):
        s = Sphere(1 * u.km)

        def func(i, e, phase):
            return 1

        # test phase = 0 special case
        phase = 0 * u.deg
        area, _ = s.integrate_i_e_phi(func, phase)
        assert u.isclose(area, 4 * pi * u.km**2)

        # test phase > 0
        phase = 30 * u.deg
        area, _ = s.integrate_i_e_phi(func, phase)
        assert u.isclose(area, 4 * pi * u.km**2)

        # test hemisphere
        def func(i, e, phase):
            return 1.0 * (e.to_value(u.rad) < pi / 2)

        area, _ = s.integrate_i_e_phi(func, phase)
        assert u.isclose(area, 2 * pi * u.km**2)

        # test cos(i < 90)
        def func(i, e, phase):
            return min_zero_cos(i)

        area, _ = s.integrate_i_e_phi(func, 0 * u.deg)
        assert u.isclose(area, pi * u.km**2)

    def test_absorption_of_sunlight(self):
        surface = SimpleSurface()
        s = Sphere(1 * u.km)
        eph = dict(rh=2 * u.au, delta=2 * u.au, phase=0 * u.deg)
        wave = np.logspace(-0.5, 0.5) * u.um

        sun = Sun.from_array(wave, 1000 * u.W / u.m**2 / u.um)

        with solar_spectrum.set(sun):
            a, _ = s.absorption_of_sunlight(wave, 0.95, surface, eph, interpolate=True)

        expected = 950 / 4 * sun.fluxd.unit * pi * s.radius**2
        assert u.allclose(a, expected)

    def test_reflected_sunlight(self):
        surface = SimpleSurface()
        s = Sphere(1 * u.km)
        eph = dict(rh=2 * u.au, delta=2 * u.au, phase=0 * u.deg)
        wave = np.logspace(-0.5, 0.5) * u.um

        sun = Sun.from_array(wave, 1000 * u.W / u.m**2 / u.um)

        with solar_spectrum.set(sun):
            r, _ = s.reflected_sunlight(wave, 0.1, surface, eph, interpolate=True)

        # (%i2) integrate (2 * %pi * cos(x)**2 * sin(x), x, 0, %pi / 2);
        #                                      2 %pi
        # (%o2)                                -----
        #                                        3
        expected = 100 / 4 * sun.fluxd.unit * 2 / 3 * pi * s.radius**2 / (2 * u.au) ** 2
        assert u.allclose(r, expected)
