# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from astropy.coordinates.errors import ConvertError
from astropy.time import Time
from astropy.wcs import WCS
from ..syndynes import (
    Syndyne,
    Syndynes,
    Synchrone,
    Synchrones,
    SourceOrbit,
    SynCollection,
    SynGenerator,
)
from ..state import State
from ..models import SolarGravity, SolarGravityAndRadiationPressure


pytest.importorskip("scipy")


def fake_state(source, betas, ages):
    """Fake dynamics: particles are displaced beta * 1e4 km from the source for
    a little more control on the tests."""
    r0 = 10000 * u.km
    v0 = 0.1 * u.km / u.s
    length = len(betas)
    return State(
        (betas * r0)[:, np.newaxis] + source.r,
        (betas * v0)[:, np.newaxis] + source.v,
        source.t - ages,
        frame=source.frame,
    )


def fake_syndynes():
    observer = State(
        [0, 1, 0] * u.au,
        [30, 0, 0] * u.km / u.s,
        Time("2025-12-29"),
        frame="icrs",
    )
    comet = State(
        [1, 0, 0] * u.au,
        [0, 30, 0] * u.km / u.s,
        observer.t,
        frame="heliocentriceclipticiau76",
    )

    betas = 10.0 ** -np.arange(4)
    ages = np.arange(10) * u.day

    items = []
    for beta in betas:
        particles = fake_state(comet, np.repeat(beta, 10), ages * 0)
        initial = fake_state(comet, np.zeros(10), ages)
        items.append(
            Syndyne(
                comet,
                beta,
                ages,
                particles.r,
                particles.v,
                particles.t,
                initial,
                observer=observer,
            )
        )

    return comet, observer, betas, ages, items


def fake_synchrones():
    observer = State(
        [0, 1, 0] * u.au,
        [30, 0, 0] * u.km / u.s,
        Time("2025-12-29"),
        frame="icrs",
    )
    comet = State(
        [1, 0, 0] * u.au,
        [0, 30, 0] * u.km / u.s,
        observer.t,
        frame="heliocentriceclipticiau76",
    )

    betas = 10.0 ** -np.arange(4)
    ages = np.arange(10) * u.day

    items = []
    for age in ages:
        particles = fake_state(comet, betas, age)
        initial = fake_state(comet, np.zeros(10), ages)
        items.append(
            Synchrone(
                comet,
                betas,
                age,
                particles.r,
                particles.v,
                particles.t,
                initial,
                observer=observer,
            )
        )

    return comet, observer, betas, ages, items


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
            [k not in eph for k in ("delta", "deltadot", "ra", "dec", "lon", "lat")]
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

    def test_getitem(self, example_syndynes):
        _, _, _, dust, _ = example_syndynes
        assert isinstance(dust.syndynes(), Syndynes)
        assert isinstance(dust.syndynes()[0], Syndyne)
        assert isinstance(dust.syndynes()[:1], Syndynes)


class TestSynCollection:
    def test_plot(self):
        """Test plot's label_format when plotting a mix of syndynes and synchrones"""

        plt = pytest.importorskip("matplotlib.pyplot")

        _, _, _, _, syndynes = fake_syndynes()
        _, _, _, _, synchrones = fake_synchrones()
        collection = SynCollection(syndynes + synchrones)

        _, ax = plt.subplots()
        collection.plot()

        for i, line in enumerate(ax.get_lines()):
            # when label="", mpl sets it to _child0, _child1, ...
            assert line.get_label().startswith("_child")


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

        # most everything is tested in TestSyndyne.test_to_ephem let's just
        # check the betas and ages
        assert np.allclose(eph["beta_rad"], betas)
        assert u.allclose(eph["age"], age)

    def test_getitem(self, example_syndynes):
        _, _, _, dust, _ = example_syndynes
        assert isinstance(dust.synchrones(), Synchrones)
        assert isinstance(dust.synchrones()[0], Synchrone)
        assert isinstance(dust.synchrones()[:1], Synchrones)


class TestSyndynes:
    def test_init(self):
        comet, _, betas, ages, items = fake_syndynes()

        # create from list
        syndynes = Syndynes(items)

        # create from syndynes object
        syndynes2 = Syndynes(syndynes)
        assert np.all(syndynes2[0].r == syndynes[0].r)

        # create from tuple
        syndynes3 = Syndynes(tuple(items))
        assert np.all(syndynes3[0].r == syndynes[0].r)

        # cannot mix syndynes and synchrones
        _, _, _, _, items2 = fake_synchrones()
        with pytest.raises(TypeError):
            Syndynes(items + items2)

        # test empty syndynes objects
        assert len(Syndynes([])) == 0

    def test_getitem(self):
        _, _, _, _, items = fake_syndynes()
        syndynes = Syndynes(items)

        # get a single Syndyne
        assert isinstance(syndynes[2], Syndyne)
        assert syndynes[2] is items[2]

        # get Syndynes from a slice
        assert len(syndynes[:3]) == 3
        assert isinstance(syndynes[:3], Syndynes)

        # get Syndynes from a tuple
        assert isinstance(syndynes[(0, 1)], Syndynes)

    def test_repr(self):
        _, _, _, _, items = fake_syndynes()
        syndynes = Syndynes(items)
        assert repr(syndynes) == "<Syndynes: betas=[1.    0.1   0.01  0.001]>"

    def test_to_ephem(self):
        comet, _, _, _, items = fake_syndynes()
        syndynes = Syndynes(items)

        eph = syndynes.to_ephem()
        assert len(eph) == 40
        assert all(eph.meta["source"]["r"] == comet.r)

        # [1] is 2nd index of 1st syndyne
        assert eph["x"][1] == items[0].x[1]

        # [13] is 4th index of 2nd syndyne
        assert eph["vz"][13] == items[1].v_z[3]

        assert len(Syndynes([]).to_ephem()) == 0

    def test_plot(self):
        plt = pytest.importorskip("matplotlib.pyplot")

        comet, observer, betas, ages, items = fake_syndynes()
        syndynes = Syndynes(items)

        _, ax = plt.subplots()
        syndynes.plot()

        # observer is sqrt(2) au away, particles are up to 1e-4 au = 14959 km =
        # 15"
        for i, line in enumerate(ax.get_lines()):
            assert line.get_label() == f"$\\beta={betas[i]}$"
            assert all(line.get_xdata() <= 15 * betas[i] * u.arcsec)
            assert all(line.get_ydata() <= 15 * betas[i] * u.arcsec)

        coords = observer.observe(comet)

        # 1 deg/pix
        wcs = WCS()
        wcs.wcs.ctype = "RA---TAN", "DEC--TAN"
        wcs.wcs.crval = coords.ra.deg, coords.dec.deg

        syndynes.plot(wcs=wcs)

        # 15" / (1 deg/pix) = 0.00417 pix
        for i, line in enumerate(ax.get_lines()[len(syndynes) :]):
            assert all(line.get_xdata() <= 4.17e-3 * betas[i])
            assert all(line.get_ydata() <= 4.17e-3 * betas[i])

        # cannot plot relative to source when observer is None
        state = fake_state(comet, [1], ages)
        syndyne = Syndyne(
            comet, [1], ages, state.r, state.v, state.t, comet, observer=None
        )
        with pytest.raises(ValueError):
            syndyne.plot()


class TestSynchrones:
    # most of the relevant tests are covered by TestSyndynes

    def test_init(self):
        _, _, _, _, items = fake_synchrones()
        synchrones = Synchrones(items)

        assert synchrones[2] is items[2]
        assert len(synchrones) == len(items)

    def test_getitem(self):
        _, _, _, _, items = fake_synchrones()
        synchrones = Synchrones(items)

        # index with integer returns one synchrone
        assert synchrones[2] is items[2]
        assert isinstance(synchrones[0], Synchrone)

        # index with slice returns synchrone collection
        assert isinstance(synchrones[:3], Synchrones)

        # index with tuple returns synchrone collection
        assert isinstance(synchrones[(0, 1)], Synchrones)

    def test_repr(self):
        _, _, _, _, items = fake_synchrones()
        synchrones = Synchrones(items)
        assert (
            repr(synchrones) == "<Synchrones: ages=[0. 1. 2. 3. 4. 5. 6. 7. 8. 9.] d>"
        )

    def test_plot(self):
        plt = pytest.importorskip("matplotlib.pyplot")

        comet, observer, betas, ages, items = fake_synchrones()
        synchrones = Synchrones(items)

        _, ax = plt.subplots()
        synchrones.plot()

        for i, line in enumerate(ax.get_lines()):
            assert line.get_label() == f"$\\Delta t={ages[i]}$"


class TestSourceOrbit:
    def test_dt(self):
        comet = State([1, 0, 0] * u.au, [0, 30, 0] * u.km / u.s, Time("2023-12-07"))
        dt = [-1, 0, 1] * u.day
        ages = -dt
        state = fake_state(comet, [0, 1, 2], ages)
        orbit = SourceOrbit(comet, dt, state.r, state.v, comet.t + dt)

        assert all(orbit.dt == dt)
        assert all(orbit.epoch == comet.t + dt)


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
        orbit = dust.source_orbit(dt)
        assert u.allclose((orbit.t - comet.t).to("s"), dt)

        solver = SolarGravity()
        for i in range(len(dt)):
            expected = solver.solve(comet, comet.t + dt[i])
            assert u.allclose(orbit[i].r, expected.r, rtol=1e-11, atol=1 * u.cm)
            assert u.allclose(orbit[i].v, expected.v, rtol=1e-11, atol=1 * u.um / u.s)
