# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from astropy.coordinates.errors import ConvertError
from astropy.time import Time
from ..syndynes import Syndyne, Synchrone, Syndynes, Synchrones, SynGenerator
from ..state import State
from ..models import SolarGravity, SolarGravityAndRadiationPressure


pytest.importorskip("scipy")


@pytest.fixture
def example_syndynes():
    comet = State(
        [1, 0, 0] * u.au,
        [0, 30, 0] * u.km / u.s,
        Time("2023-12-07"),
        frame="heliocentriceclipticiau76",
    )
    betas = [0.1]
    ages = [1, 10, 100] * u.d
    observer = State(
        [0, 1, 0] * u.au, [30, 0, 0] * u.km / u.s, comet.t, frame=comet.frame
    )
    dust = SynGenerator(comet, betas, ages, observer=observer)
    return comet, betas, ages, dust, observer


class TestSyndyne:
    def test_init(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes

        syndyne = Syndyne(
            comet,
            betas[0],
            ages,
            dust.particles.r[:3],
            dust.particles.v[:3],
            dust.particles.t[:3],
            dust.initial_states[:3],
        )

        assert syndyne.source is comet
        assert syndyne.beta == betas[0]
        assert np.all(syndyne.ages == ages)
        assert np.all(dust.particles.r[:3] == syndyne.r)
        assert np.all(dust.particles.v[:3] == syndyne.v)
        assert np.all(dust.particles.t[:3] == syndyne.t)
        assert syndyne.coords is None
        assert syndyne.observer is None

        syndyne = Syndyne(
            comet,
            betas[0],
            ages,
            dust.particles.r[:3],
            dust.particles.v[:3],
            dust.particles.t[:3],
            dust.initial_states[:3],
            observer=observer,
        )
        assert np.all(syndyne.coords == observer.observe(dust.particles[:3]))
        assert syndyne.observer is observer

    def test_to_ephem(self, example_syndynes):
        comet, _, _, _, observer = example_syndynes

        ages = [1, 10] * u.day
        particles = State(
            [[0, 1, 2], [0, 10, 20]] * u.au,
            [[2, 1, 0], [20, 10, 0]] * u.km / u.s,
            comet.t,
            frame=comet.frame,
        )
        rh = np.sqrt([5, 500]) * u.au
        rhat = (particles.r.T / rh).T
        rdot = np.sqrt(np.sum((particles.v * rhat) ** 2, 1))

        d = particles - observer
        delta = np.sqrt(np.sum(d.r**2, 1))
        deltahat = (d.r.T / delta).T
        deltadot = np.sqrt(np.sum((d.v * deltahat) ** 2, 1))

        coords = observer.observe(particles)

        initial_states = State(
            [[0.1, 0.2, 0.3], [1.1, 1.2, 1.3]] * u.au,
            [[1.1, 1.2, 1.3], [10.1, 10.2, 10.3]] * u.km / u.s,
            comet.t,
            frame=comet.frame,
        )

        syndyne = Syndyne(
            comet,
            0.1,
            ages,
            particles.r,
            particles.v,
            particles.t,
            initial_states,
            observer=observer,
        )

        eph = syndyne.to_ephem()
        assert np.allclose(eph["beta_rad"], 0.1)
        assert u.allclose(eph["age"], ages, atol=1e-12 * u.s)
        assert np.allclose(eph["date"].et, comet.t.et, atol=1e-12)
        assert u.allclose(eph["r"], rh, atol=1e-12 * u.au)
        assert u.allclose(eph["rdot"], rdot, atol=1e-12 * u.km / u.s)
        assert u.allclose(eph["delta"], delta, atol=1e-12 * u.au)
        assert u.allclose(eph["deltadot"], deltadot, atol=1e-12 * u.km / u.s)
        assert u.allclose(eph["lon"], coords.lon)
        assert u.allclose(eph["lat"], coords.lat)
        assert u.allclose(eph["coords"].lon, coords.lon)
        assert u.allclose(eph["coords"].lat, coords.lat)
        assert u.allclose(eph["x"], particles.x)
        assert u.allclose(eph["y"], particles.y)
        assert u.allclose(eph["z"], particles.z)
        assert u.allclose(eph["vx"], particles.v_x)
        assert u.allclose(eph["vy"], particles.v_y)
        assert u.allclose(eph["vz"], particles.v_z)
        assert u.allclose(eph["x initial"], initial_states.x)
        assert u.allclose(eph["y initial"], initial_states.y)
        assert u.allclose(eph["z initial"], initial_states.z)
        assert u.allclose(eph["vx initial"], initial_states.v_x)
        assert u.allclose(eph["vy initial"], initial_states.v_y)
        assert u.allclose(eph["vz initial"], initial_states.v_z)
        assert np.allclose(eph["t initial"].et, initial_states.t.et)

        assert len(eph.meta["source"]["r"]) == 3
        assert u.allclose(eph.meta["source"]["r"], comet.r)
        assert u.allclose(eph.meta["source"]["v"], comet.v)
        assert np.allclose(eph.meta["source"]["t"].et, comet.t.et)
        assert eph.meta["source"]["frame"] == comet.frame

        assert u.allclose(eph.meta["observer"]["r"], observer.r)
        assert u.allclose(eph.meta["observer"]["v"], observer.v)
        assert np.allclose(eph.meta["observer"]["t"].et, observer.t.et)
        assert eph.meta["observer"]["frame"] == observer.frame

        # test conditional branches
        # no observer, relative time
        comet2 = State(comet.r, comet.v, 0 * u.s, frame=None)
        initial_states2 = State(initial_states.r, initial_states.v, -ages)
        syndyne = Syndyne(
            comet2,
            0.1,
            ages,
            particles.r,
            particles.v,
            comet2.t,
            initial_states2,
        )

        eph = syndyne.to_ephem()
        assert "date" not in eph
        assert np.all(eph["t_relative"] == 0 * u.s)
        assert eph.meta["observer"] is None
        assert np.all(
            [
                k not in eph
                for k in ("delta", "deltadot", "ra", "dec", "lon", "lat", "coords")
            ]
        )

    def test_init_1d(self, example_syndynes):
        comet, _, _, _, _ = example_syndynes

        ages = [1, 10] * u.day
        particles = State(
            [[0, 1, 2], [0, 10, 20]] * u.au,
            [[2, 1, 0], [20, 10, 0]] * u.km / u.s,
            comet.t,
            frame=comet.frame,
        )

        initial_states = State(
            [[0.1, 0.2, 0.3], [1.1, 1.2, 1.3]] * u.au,
            [[1.1, 1.2, 1.3], [10.1, 10.2, 10.3]] * u.km / u.s,
            comet.t,
            frame=comet.frame,
        )

        # 1D r and v
        with pytest.raises(
            ValueError, match="Syndyne only supports 2 dimensional r and v vectors"
        ):
            Syndyne(
                comet,
                0.1,
                ages,
                [0, 0, 1] * u.km,
                particles.v,
                comet.t,
                initial_states,
            )

        with pytest.raises(
            ValueError, match="Syndyne only supports 2 dimensional r and v vectors"
        ):
            Syndyne(
                comet,
                0.1,
                ages,
                particles.r,
                [0, 0, 1] * u.km / u.s,
                comet.t,
                initial_states,
            )


class TestSynchrone:
    def test_init(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes

        synchrone = Synchrone(
            comet,
            betas,
            ages[0],
            dust.particles.r[:1],
            dust.particles.v[:1],
            dust.particles.t[:1],
            dust.initial_states[:1],
        )

        assert synchrone.source is comet
        assert np.all(synchrone.betas == betas)
        assert synchrone.age == ages[0]
        assert np.all(dust.particles.r[:1] == synchrone.r)
        assert np.all(dust.particles.v[:1] == synchrone.v)
        assert np.all(dust.particles.t[:1] == synchrone.t)
        assert synchrone.coords is None
        assert synchrone.observer is None

        synchrone = Synchrone(
            comet,
            betas,
            ages[0],
            dust.particles.r[:1],
            dust.particles.v[:1],
            dust.particles.t[:1],
            dust.initial_states[:1],
            observer=observer,
        )
        assert np.all(synchrone.coords == observer.observe(dust.particles[:1]))
        assert synchrone.observer is observer

    def test_to_ephem(self, example_syndynes):
        comet, _, _, _, observer = example_syndynes

        betas = [0.1, 1]
        age = 10 * u.day
        particles = State(
            [[0, 1, 2], [0, 10, 20]] * u.au,
            [[2, 1, 0], [20, 10, 0]] * u.km / u.s,
            comet.t,
            frame=comet.frame,
        )

        initial_states = State(
            [[0.1, 0.2, 0.3], [1.1, 1.2, 1.3]] * u.au,
            [[1.1, 1.2, 1.3], [10.1, 10.2, 10.3]] * u.km / u.s,
            comet.t,
            frame=comet.frame,
        )

        synchrone = Synchrone(
            comet,
            betas,
            age,
            particles.r,
            particles.v,
            particles.t,
            initial_states,
            observer=observer,
        )

        eph = synchrone.to_ephem()

        # most everything is tested in TestSyndyne.test_to_ephem
        # we just need to make sure the betas and ages are right
        assert np.allclose(eph["beta_rad"], betas)
        assert u.allclose(eph["age"], age)


def random_state(length=1):
    return State(
        np.random.rand(length, 3).squeeze() * u.au,
        np.random.rand(length, 3).squeeze() * u.km / u.s,
        np.random.rand() * u.day,
    )


def test_syndynes():
    comet = random_state()
    betas = np.linspace(0, 1, 3)
    ages = np.arange(10) * u.day

    items = []
    for beta in betas:
        particles = random_state(len(ages))
        initial = random_state(len(ages))
        items.append(
            Syndyne(
                comet,
                beta,
                ages,
                particles.r,
                particles.v,
                particles.t,
                initial,
            )
        )

    syndynes = Syndynes(items)
    assert syndynes[2] is items[2]
    assert len(syndynes[:3]) == 3
    assert repr(syndynes) == "<Syndynes: betas=[0.  0.5 1. ]>"

    eph = syndynes.to_ephem()
    assert len(eph) == 30
    assert all(eph.meta["source"]["r"] == comet.r)

    syndynes2 = Syndynes(syndynes)
    assert np.all(syndynes2[0].r == syndynes[0].r)
    syndynes3 = Syndynes(tuple(items))
    assert np.all(syndynes3[0].r == syndynes[0].r)

    particles = random_state(len(betas))
    synchrone = Synchrones(
        [
            Synchrone(
                comet, betas, ages[0], particles.r, particles.v, particles.t, particles
            )
        ]
    )
    with pytest.raises(TypeError):
        Syndynes(items + [synchrone])

    # test emptpy syndynes objects
    assert len(Syndynes([])) == 0
    assert len(Syndynes([]).to_ephem()) == 0


def test_synchrones():
    comet = random_state()
    betas = np.linspace(0, 1, 3)
    ages = np.arange(10) * u.day

    items = []
    for age in ages:
        particles = random_state(len(ages))
        initial = random_state(len(ages))
        items.append(
            Synchrone(
                comet,
                betas,
                age,
                particles.r,
                particles.v,
                particles.t,
                initial,
            )
        )

    synchrones = Synchrones(items)
    assert synchrones[2] is items[2]
    assert len(synchrones[:3]) == 3
    assert repr(synchrones) == "<Synchrones: ages=[0. 1. 2. 3. 4. 5. 6. 7. 8. 9.] d>"


class TestSynGenerator:
    def test_init(self):
        comet = State([1, 0, 0] * u.au, [0, 30, 0] * u.km / u.s, Time("2023-12-07"))
        betas = [0.1]
        ages = [1, 10, 100] * u.d

        # no observer
        dust = SynGenerator(comet, betas, ages)
        assert dust.observer is None

        # observer and comet are both in the "Arbitrary" frame
        observer = State([0, 1, 0] * u.au, [30, 0, 0] * u.km / u.s, comet.t)
        dust = SynGenerator(comet, betas, ages, observer=observer)
        dust.syndyne(0)

        # cannot convert comet frame to observer frame:
        comet = State(comet.r, comet.v, comet.t, frame="heliocentriceclipticiau76")
        dust = SynGenerator(comet, betas, ages, observer=observer)
        with pytest.raises(
            ConvertError,
            match="Cannot transform from.*HeliocentricEclipticIAU76.*to.*ArbitraryFrame",
        ):
            dust.syndyne(0)

        # fix observer frame
        observer = State(observer.r, observer.v, observer.t, comet.frame)
        dust = SynGenerator(comet, betas, ages, observer=observer)
        dust.syndyne(0)

        # test source check
        with pytest.raises(ValueError, match="Only one source state vector allowed"):
            SynGenerator(State.from_states([comet, comet]), betas, ages)

    def test_at_epochs(self, example_syndynes):
        comet, betas, ages, _, observer = example_syndynes

        # produce a synchrone 6 days before the observational epoch
        date = Time("2023-12-01")
        dust = SynGenerator.at_epochs(comet, betas, date)
        assert dust.ages == 6 * u.day
        synchrone = dust.synchrone(0)
        assert synchrone.epoch.et == date.et

    def test_repr(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes
        assert (
            repr(dust)
            == """<SynGenerator:
 betas
    [0.1]
 ages
    [  1.  10. 100.] d>"""
        )

    def test_initialize_states(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes

        solver = SolarGravity()
        for i, age in enumerate(ages):
            initial = dust.initial_states[i]
            expected = solver.solve(comet, comet.t - age)
            assert u.allclose(initial.r, expected.r, atol=1 * u.cm, rtol=1e-11)
            assert u.allclose(initial.v, expected.v, atol=1 * u.um / u.s, rtol=1e-11)
            assert u.isclose(
                (initial.t - expected.t).to("s"), 0 * u.s, atol=1 * u.ns, rtol=1e-11
            )

    def test_solve(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes

        solver = SolarGravityAndRadiationPressure()
        for syndyne in dust.syndynes():
            for j, age in enumerate(syndyne.ages):
                initial = solver.solve(comet, comet.t - age, 0)
                expected = solver.solve(initial, comet.t, syndyne.beta)
                assert u.allclose(syndyne[j].r, expected.r, atol=1 * u.cm, rtol=1e-11)
                assert u.allclose(
                    syndyne[j].v, expected.v, atol=1 * u.um / u.s, rtol=1e-11
                )
                assert np.allclose(syndyne.t.et, expected.t.et, rtol=1e-11)

    def test_syndynes(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes

        for syndyne in dust.syndynes():
            assert syndyne.beta == betas[0]
            assert np.allclose(syndyne.t.et, comet.t.et, rtol=1e-11)
            assert np.allclose(syndyne.coords.obstime.et, comet.t.et, rtol=1e-11)

            # State.observe is already tested, so being generous here
            assert u.allclose(
                syndyne.coords.lon, [315, 315, 314] * u.deg, atol=0.1 * u.deg
            )
            assert u.allclose(syndyne.coords.lat, 0 * u.deg)

    def test_synchrones(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes

        expected_lon = [315, 315, 314] * u.deg
        for i, synchrone in enumerate(dust.synchrones()):
            assert synchrone.age == ages[i]
            assert np.allclose(synchrone.t.et, comet.t.et, rtol=1e-11)
            assert np.allclose(synchrone.coords.obstime.et, comet.t.et, rtol=1e-11)

            # State.observe is already tested, so a generous test here
            assert u.allclose(synchrone.coords.lon, expected_lon[i], atol=0.1 * u.deg)
            assert u.allclose(synchrone.coords.lat, 0 * u.deg)

    def test_orbit(self, example_syndynes):
        comet, betas, ages, dust, observer = example_syndynes

        dt = [-1, 0, 1] * u.d
        orbit, coords = dust.source_orbit(dt)
        assert u.allclose((orbit.t - comet.t).to("s"), dt)

        solver = SolarGravity()
        for i in range(len(dt)):
            expected = solver.solve(comet, comet.t + dt[i])
            assert u.allclose(orbit[i].r, expected.r, rtol=1e-11, atol=1 * u.cm)
            assert u.allclose(orbit[i].v, expected.v, rtol=1e-11, atol=1 * u.um / u.s)

        # test without observer
        dust.observer = None
        orbit2 = dust.source_orbit(dt)
        assert np.allclose(orbit2.r, orbit.r)
