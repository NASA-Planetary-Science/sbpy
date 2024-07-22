# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
import pytest

import numpy as np
import astropy.units as u
from astropy.time import Time
import astropy.constants as const
from erfa import ErfaWarning

from ... import time  # for ephemeris time
from ..state import State
from ..models import (
    FreeExpansion,
    SolarGravity,
    SolarGravityAndRadiationPressure,
    SolverFailed,
)


scipy = pytest.importorskip("scipy")
scipy_version = [int(x) for x in scipy.__version__.split(".")]


def test_spice_prop2b():
    """Test case from SPICE NAIF toolkit prop2b, v2.2.0

    State at t0:
    R   (km):          0.00000   70710678.11865   70710678.11865
    V (km/s):          0.00000         -0.04464          0.04464

    State at tau/2:
    R   (km):         -0.00000  -70710678.11865  -70710678.11865
    V (km/s):          0.00000          0.04464         -0.04464

    """

    class EarthGravity(SolarGravity):
        _GM = 3.9860043543609598e5

    solver = EarthGravity()

    r1 = 1e8
    s = np.sqrt(solver.GM.to_value("km3/s2") / r1)
    half_period = np.pi * r1 / s

    r = [0, r1 / np.sqrt(2), r1 / np.sqrt(2)] * u.km
    v = [0, -s / np.sqrt(2), s / np.sqrt(2)] * u.km / u.s

    initial = State(r, v, Time("2023-01-01"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        t_f = initial.t + half_period * u.s
        final = solver.solve(initial, t_f)

    assert u.allclose(initial.r, [0, 70710678.11865, 70710678.11865] * u.km, rtol=1e-11)
    assert u.allclose(
        initial.v,
        [0, -0.04464, 0.04464] * u.km / u.s,
        rtol=1e-11,
        atol=0.00001 * u.km / u.s,
    )
    assert u.allclose(final.r, -r, rtol=1e-11)
    assert u.allclose(final.v, -v, rtol=1e-11, atol=0.00001 * u.km / u.s)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        t_f = initial.t + 5 * half_period * u.s
        final = solver.solve(initial, t_f)

    assert u.allclose(initial.r, [0, 70710678.11865, 70710678.11865] * u.km, rtol=1e-11)
    assert u.allclose(
        initial.v,
        [0, -0.04464, 0.04464] * u.km / u.s,
        rtol=1e-11,
        atol=0.00001 * u.km / u.s,
    )
    assert u.allclose(final.r, -r, rtol=1e-11)
    assert u.allclose(final.v, -v, rtol=1e-11)


class TestFreeExpansion:
    def test(self):
        r = [0, 1e6, 0] * u.km
        v = [0, -1, 1] * u.km / u.s

        solver = FreeExpansion()

        initial = State(r, v, Time("2023-01-01"))
        t_f = initial.t + 1e6 * u.s
        final = solver.solve(initial, t_f)

        assert u.allclose(final.r, [0, 0, 1e6] * u.km, atol=4e-4 * u.km)
        assert u.allclose(final.v, [0, -1, 1] * u.km / u.s)

        solver = FreeExpansion(method="Radau")

        initial = State(r, v, Time("2023-01-01"))
        final = solver.solve(initial, t_f)

        assert u.allclose(final.r, [0, 0, 1e6] * u.km, atol=4e-4 * u.km)
        assert u.allclose(final.v, [0, -1, 1] * u.km / u.s)

    def test_arbitrary_time(self):
        r = [0, 1e6, 0] * u.km
        v = [0, -1, 1] * u.km / u.s

        solver = FreeExpansion()

        initial = State(r, v, 0 * u.s)
        t_f = 1e6 * u.s
        final = solver.solve(initial, t_f)

        assert u.allclose(final.r, [0, 0, 1e6] * u.km, atol=2e-7 * u.km)
        assert u.allclose(final.v, [0, -1, 1] * u.km / u.s)
        assert final.t == t_f


class TestSolarGravity:
    @pytest.mark.parametrize("r1_au", ([0.3, 1, 3, 10, 30]))
    def test_circular_orbit(self, r1_au):
        solver = SolarGravity()

        r1 = r1_au * u.au.to("km")
        s = np.sqrt(solver.GM.to_value("km3/s2") / r1)
        half_period = np.pi * r1 / s

        r = [0, r1 / np.sqrt(2), r1 / np.sqrt(2)] * u.km
        v = [0, -s / np.sqrt(2), s / np.sqrt(2)] * u.km / u.s

        initial = State(r, v, Time("2023-01-01"))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ErfaWarning)
            t_f = initial.t + half_period * u.s
            final = solver.solve(initial, t_f)

        assert u.allclose(final.r, -initial.r, atol=10 * u.m)
        assert u.allclose(final.v, -initial.v, atol=1 * u.um / u.s)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ErfaWarning)
            t_f = initial.t + 2 * half_period * u.s
            final = solver.solve(initial, t_f)

        assert u.allclose(final.r, initial.r, atol=60 * u.m)
        assert u.allclose(final.v, initial.v, atol=1 * u.um / u.s)

        solver = SolarGravity(method="Radau")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ErfaWarning)
            t_f = initial.t + half_period * u.s
            final = solver.solve(initial, t_f)

        assert u.allclose(final.r, -initial.r, atol=1 * u.m)
        assert u.allclose(final.v, -initial.v, atol=1 * u.um / u.s)

    def test_GM(self):
        solver = SolarGravity()
        assert u.isclose(solver.GM, const.G * const.M_sun, rtol=1e-12)

    @pytest.mark.skipif(
        "scipy_version[0] < 2 and scipy_version[1] < 8", reason="requires scipy>=1.10"
    )
    def test_solverfailed(self):
        r = [0, 1, 0] * u.au
        v = [0, -1, 1] * u.km / u.s

        # force a solution failure with random derivatives
        class RandomGravity(SolarGravity):
            @classmethod
            def dx_dt(cls, t, rv, *args):
                return np.random.rand(6)

        solver = RandomGravity()

        initial = State(r, v, Time("2023-01-01"))
        t_f = initial.t + 1e6 * u.s
        with pytest.raises(SolverFailed):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                solver.solve(initial, t_f)


class TestSolarGravityAndRadiationPressure:
    def test_reduced_gravity(self):
        """Radiation pressure is essentially a reduced gravity problem."""
        solver = SolarGravityAndRadiationPressure()

        r1 = 1e8
        s = np.sqrt(solver.GM.to_value("km3/s2") / r1)

        initial = State([0, 0, r1] * u.km, [0, s, 0] * u.km / u.s, Time("2023-01-01"))
        t_f = initial.t + 1e6 * u.s
        beta = 0.1
        final1 = solver.solve(initial, t_f, beta)

        class ReducedGravity(SolarGravity):
            _GM = (1 - beta) * SolarGravity._GM

        solver = ReducedGravity()
        final2 = solver.solve(initial, t_f)

        assert u.allclose(final1.r, final2.r)
        assert u.allclose(final1.v, final2.v)

        solver = SolarGravityAndRadiationPressure(method="Radau")
        final2 = solver.solve(initial, t_f, beta)

        assert u.allclose(final1.r, final2.r)
        assert u.allclose(final1.v, final2.v)
