# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import HeliocentricEclipticIAU76, SkyCoord

from ...data import Ephem
from ..state import ArbitraryFrame, State


class TestState:
    def test_init(self):
        t = Time("2022-08-02")
        state = State(
            [1, 1, 0] * u.au,
            [30, 0, 0] * u.km / u.s,
            t,
            frame="heliocentriceclipticiau76",
        )

        # test properties
        assert u.allclose(state.r, [1, 1, 0] * u.au)
        assert u.allclose(state.v, [30000, 0, 0] * u.m / u.s)
        assert np.isclose((state.t - t).jd, 0)
        assert not state.arbitrary_time

        # initialize with relative time
        state = State(
            [1, 1, 0] * u.au,
            [30, 0, 0] * u.km / u.s,
            0 * u.s,
        )
        assert state.arbitrary_time
        assert u.isclose(state.t, 0 * u.s)

        with pytest.raises(
            TypeError,
            match="`State` only supports time as a quantity with `ArbitraryFrame`",
        ):
            State(state.r, state.v, 0 * u.s, frame="icrs")

    def test_init_shape_mismatch(self):
        with pytest.raises(ValueError):
            State(
                [1, 1, 0] * u.au,
                [[0, 0, 30], [30, 0, 0]] * u.km / u.s,
                Time("2022-08-02"),
                frame="heliocentriceclipticiau76",
            )

        with pytest.raises(ValueError, match="only supports 1 and 2 dimensional"):
            State(
                [[[1], [1], [0]]] * u.au,
                [[[0], [0], [30]]] * u.km / u.s,
                Time("2022-08-02"),
            )

    def test_rv_shapes(self):
        good = [1, 2, 3]
        bad = [1, 2, 3, 4]
        t = Time("2023-12-03")

        State(good * u.km, good * u.km / u.s, t)

        with pytest.raises(ValueError):
            State(good * u.km, bad * u.km / u.s, t)

        with pytest.raises(ValueError):
            State(bad * u.km, good * u.km / u.s, t)

        bad = [[[1, 2, 3], [1, 2, 3]]]

        with pytest.raises(ValueError):
            State(good * u.km, bad * u.km / u.s, t)

        with pytest.raises(ValueError):
            State(bad * u.km, good * u.km / u.s, t)

    def test_frame_inputs(self):
        data = ([1, 0, 0] * u.au, [0, 30, 0] * u.km / u.s, Time("2024-01-21"))

        assert isinstance(State(*data, frame=None).frame, ArbitraryFrame)

        assert isinstance(
            State(*data, frame="heliocentriceclipticiau76").frame,
            HeliocentricEclipticIAU76,
        )

        with pytest.raises(ValueError, match="Invalid frame"):
            State(*data, frame="asdf")

        assert isinstance(
            State(*data, frame=HeliocentricEclipticIAU76).frame,
            HeliocentricEclipticIAU76,
        )

        assert isinstance(
            State(*data, frame=HeliocentricEclipticIAU76()).frame,
            HeliocentricEclipticIAU76,
        )

    def test_repr(self):
        state = State(
            [1, 1, 0] * u.au,
            [30, 0, 0] * u.km / u.s,
            Time("2022-08-02"),
            frame="heliocentriceclipticiau76",
        )
        assert (
            repr(state)
            == """<State (<HeliocentricEclipticIAU76 Frame (obstime=2022-08-02 00:00:00.000)>):
 r
  [1. 1. 0.] AU
 v
  [30.  0.  0.] km / s
 t
  2022-08-02 00:00:00.000>"""
        )

    def test_len(self):
        state = State(
            [1, 1, 0] * u.au,
            [30, 0, 0] * u.km / u.s,
            Time("2022-08-02"),
            frame="heliocentriceclipticiau76",
        )
        assert len(state) == 1

        state = State(
            [[1, 1, 0], [1, 1, 0]] * u.au,
            [[30, 0, 0], [30, 0, 0]] * u.km / u.s,
            Time(["2022-08-02", "2023-08-02"]),
            frame="heliocentriceclipticiau76",
        )
        assert len(state) == 2

    def test_getitem(self):
        state = State(
            [[1, 1, 0], [1, 1, 1]] * u.au,
            [[30, 0, 0], [0, 30, 0]] * u.km / u.s,
            Time(["2022-08-02", "2023-08-02"]),
            frame="heliocentriceclipticiau76",
        )
        assert u.allclose(state[0].r, [1, 1, 0] * u.au)
        assert u.allclose(state[0].v, [30, 0, 0] * u.km / u.s)
        assert np.isclose((state.t[0] - Time("2022-08-02")).jd, 0)

        assert u.allclose(state[1].r, [1, 1, 1] * u.au)
        assert u.allclose(state[1].v, [0, 30, 0] * u.km / u.s)
        assert np.isclose((state.t[1] - Time("2023-08-02")).jd, 0)

    def test_operators(self):
        t = Time("2022-08-02")
        a = State([1, 2, 3] * u.au, [10, 20, -30] * u.km / u.s, t)
        b = State([4, 5, 6] * u.au, [100, 200, -300] * u.km / u.s, Time("2023-08-02"))

        x = a + b
        assert u.allclose(x.r, [5, 7, 9] * u.au)
        assert u.allclose(x.v, [110, 220, -330] * u.km / u.s)
        assert np.isclose((x.t - a.t).jd, 0)

        x = -a
        assert u.allclose(x.r, [-1, -2, -3] * u.au)
        assert u.allclose(x.v, [-10, -20, 30] * u.km / u.s)

        x = a - b
        assert u.allclose(x.r, [-3, -3, -3] * u.au)
        assert u.allclose(x.v, [-90, -180, 270] * u.km / u.s)
        assert np.isclose((x.t - a.t).jd, 0)

    def test_abs(self):
        t = Time("2022-08-02")
        a = State([1, 2, 3] * u.au, [10, 20, -30] * u.km / u.s, t)
        b = State.from_states(
            [
                a,
                State([4, 5, 6] * u.au, [100, 200, -300] * u.km / u.s, t),
            ]
        )

        def length(*x):
            return np.sqrt(np.sum(np.array(x) ** 2))

        r, v = abs(a)
        assert u.isclose(r, length(1, 2, 3) * u.au, atol=1 * u.um, rtol=1e-10)
        assert u.isclose(
            v, length(10, 20, 30) * u.km / u.s, atol=1 * u.um / u.s, rtol=1e-10
        )

        r, v = abs(b)
        assert r.shape == (2,)
        assert v.shape == (2,)
        assert u.allclose(
            r, [length(1, 2, 3), length(4, 5, 6)] * u.au, atol=1 * u.um, rtol=1e-10
        )
        assert u.allclose(
            v,
            [length(10, 20, 30), length(100, 200, 300)] * u.km / u.s,
            atol=1 * u.um / u.s,
            rtol=1e-10,
        )

    def test_vector_properties(self):
        state = State(
            [1, 2, 4] * u.au,
            [30, -10, 5] * u.km / u.s,
            Time("2022-08-02"),
            frame="heliocentriceclipticiau76",
        )
        assert u.isclose(state.x, 1 * u.au)
        assert u.isclose(state.y, 2 * u.au)
        assert u.isclose(state.z, 4 * u.au)
        assert u.isclose(state.v_x, 30 * u.km / u.s)
        assert u.isclose(state.v_y, -10 * u.km / u.s)
        assert u.isclose(state.v_z, 5 * u.km / u.s)
        assert np.allclose(
            state.rv, [1.49597871e8, 2 * 1.49597871e8, 4 * 1.49597871e8, 30, -10, 5]
        )

        state = State(
            [[1, 2, 4], [7, 8, 9]] * u.au,
            [[30, -10, 5], [5, 4, 6]] * u.km / u.s,
            [Time("2022-08-02")] * 2,
            frame="heliocentriceclipticiau76",
        )
        assert u.allclose(state.x, [1, 7] * u.au)
        assert u.allclose(state.y, [2, 8] * u.au)
        assert u.allclose(state.z, [4, 9] * u.au)
        assert u.allclose(state.v_x, [30, 5] * u.km / u.s)
        assert u.allclose(state.v_y, [-10, 4] * u.km / u.s)
        assert u.allclose(state.v_z, [5, 6] * u.km / u.s)
        assert np.allclose(
            state.rv,
            [
                [1.49597871e8, 2 * 1.49597871e8, 4 * 1.49597871e8, 30, -10, 5],
                [7 * 1.49597871e8, 8 * 1.49597871e8, 9 * 1.49597871e8, 5, 4, 6],
            ],
        )

    def test_from_skycoord(self):
        """Test initialization from RA, Dec, distance.

        Coordinates are 2P/Encke at 2022-08-02 UTC

        Initialize using RA, Dec, and distance from Horizons.  These are
        light travel time corrected quantities.

        Compare to vectors from Horizons, explicitly adjusted for light travel
        time.
        RA, Dec is defined in the ICRF, so get barycentric coordinates from
        Horizons.

            from astropy.time import Time
            from sbpy.data import Ephem

            t = Time("2022-08-02")  # UTC
            eph = Ephem.from_horizons("2P", id_type="designation",
                                      epochs=t,
                                      closest_apparition=True,
                                      location="@0",
                                      extra_precision=True)
            [eph[k] for k in ("ra", "dec", "Delta", "RA*cos(Dec)_rate",
            "Dec_rate", "delta_rate")]

        [<MaskedQuantity [348.37706151] deg>,
         <MaskedQuantity [-1.8630357] deg>,
         <MaskedQuantity [3.90413272] AU>,
         <MaskedQuantity [6.308291] arcsec / h>,
         <MaskedQuantity [4.270122] arcsec / h>,
         <MaskedQuantity [-4.198173] km / s>]

            from astroquery.jplhorizons import Horizons
            import astropy.units as u
            from astropy.constants import c

            sun_vectors = Horizons("Sun", epochs=t.tdb.jd, location="@0").vectors(
                closest_apparition=True, refplane="ecliptic")

            ltt = eph["lighttime"]
            for i in range(3):
                t_delayed = t - ltt
                comet_vectors = Horizons("2P", id_type="designation", epochs=t_delayed.tdb.jd,
                    location="@0").vectors(closest_apparition=True, refplane="ecliptic")
                ltt = np.sqrt(np.sum([float(comet_vectors[x] - sun_vectors[x])**2 for x in "xyz"])) * u.au / c

            print([float(comet_vectors[x] - sun_vectors[x]) for x in "xyz"])
            print([float(comet_vectors[f"v{x}"] - sun_vectors[f"v{x}"]) for x in "xyz"])

        [3.831111550969744, -0.773239683630105, 0.19606224529093624]
        [-0.0017336347284832058, 0.003823198388051765, 0.0005453275902763597]

        """

        # barycentric coordinates
        coords = SkyCoord(
            ra=348.37706151 * u.deg,
            dec=-1.8630357 * u.deg,
            distance=3.90413272 * u.au,
            pm_ra_cosdec=6.308291 * u.arcsec / u.hr,
            pm_dec=4.270122 * u.arcsec / u.hr,
            radial_velocity=-4.198173 * u.km / u.s,
            obstime=Time("2022-08-02"),
        ).transform_to("heliocentriceclipticiau76")
        state = State.from_skycoord(coords)

        # heliocentric vectors, light travel time corrected
        r = [3.831111550969744, -0.773239683630105, 0.19606224529093624] * u.au
        v = (
            [-0.0017336347284832058, 0.003823198388051765, 0.0005453275902763597]
            * u.au
            / u.day
        )

        # if DE440s is used for the SkyCoord frame transformation, then the
        # agreement is better than 40 km
        assert u.allclose(state.r, r, atol=130 * u.km, rtol=1e-9)
        assert u.allclose(state.v, v, atol=2 * u.mm / u.s, rtol=1e-9)

    def test_to_skycoord(self):
        """Test conversion of coords to internal state and back to coords."""
        coords = SkyCoord(
            ra=348.37706 * u.deg,
            dec=-1.86304 * u.deg,
            distance=3.90413464 * u.au,
            pm_ra_cosdec=6.308283 * u.arcsec / u.hr,
            pm_dec=4.270114 * u.arcsec / u.hr,
            radial_velocity=-4.19816 * u.km / u.s,
            obstime=Time("2022-08-02"),
        )

        new_coords = State.from_skycoord(coords).to_skycoord()
        new_coords.representation_type = "spherical"

        assert u.isclose(new_coords.ra, coords.ra)
        assert u.isclose(new_coords.dec, coords.dec)
        assert u.isclose(new_coords.distance, coords.distance)
        assert u.isclose(new_coords.pm_ra * np.cos(new_coords.dec), coords.pm_ra_cosdec)
        assert u.isclose(new_coords.pm_dec, coords.pm_dec)
        assert u.isclose(new_coords.radial_velocity, coords.radial_velocity)
        assert np.isclose((coords.obstime - new_coords.obstime).jd, 0)

        # test skycoords representation using string
        new_coords = State(
            [1, 2, 3] * u.km,
            [4, 5, 6] * u.km / u.s,
            coords.obstime,
            frame="heliocentriceclipticiau76",
        ).to_skycoord()
        assert isinstance(new_coords.frame, HeliocentricEclipticIAU76)

    def test_observe(self):
        """Observe comet Encke from the Earth.

        Get heliocentric vectors for both (queries run 2024 Jan 22):

            from astropy.time import Time
            import astropy.constants as const
            from astroquery.jplhorizons import Horizons

            def print_rect(data):
                print([float(data[x]) for x in "xyz"])
                print([float(data[f"v{x}"]) for x in "xyz"])

            # Horizons.vectors uses TDB
            t = Time("2022-08-02", scale="tdb")
            q = Horizons(399,
                        epochs=t.tdb.jd,
                        location="@10")
            earth_data = q.vectors(refplane="ecliptic", aberrations="geometric", cache=False)
            print_rect(earth_data)

        [0.644307773374595, -0.7841972063903224, 3.391591147538755e-05]
        [0.01301350774835054, 0.01086343081039913, -2.045610960977488e-07]

            q = Horizons("2P",
                         id_type="designation",
                         epochs=t.tdb.jd,
                         location="@10")
            comet_data = q.vectors(refplane="ecliptic", closest_apparition=True, aberrations="geometric", cache=False)
            print_rect(comet_data)

        [3.831073725981583, -0.7731565414836856, 0.1960741286643133]
        [-0.001734044145851301, 0.003823283255413869, 0.0005453057034865457]

        Compare to Horizons's ICRS coordinates, adjusted for light travel time

            # Horizons.ephemerides uses UTC
            delta = np.sqrt(np.sum([(earth_data[x] - comet_data[x])**2 for x in "xyz"])) * u.au
            q = Horizons("2P",
                         id_type="designation",
                         epochs=t.utc.jd,
                         location="500")
            comet_eph = q.ephemerides(closest_apparition=True, extra_precision=True, cache=False)
            print(comet_eph["RA", "DEC", "delta"])

              RA          DEC          delta
             deg          deg            AU
        ------------- ----------- ----------------
        358.779201434 3.307642087 3.19284031826821

        Compare to spiceypy

            import spiceypy
            import sbpy.time

            heclip = [
                3.831073725981583 - 0.644307773374595,
                -0.7731565414836856 - -0.7841972063903224,
                0.1960741286643133 - -2.045610960977488e-07,
            ]
            t = Time("2022-08-02", scale="tdb")
            xmat = spiceypy.pxform("ECLIPJ2000", "J2000", t.et)
            icrf = spiceypy.mxv(xmat, heclip)
            delta, ra, dec = spiceypy.recrad(icrf)
            delta

        3.192811343804468

            np.degrees((ra, dec))

        array([358.78003306,   3.30890361])

        """

        frame = "heliocentriceclipticiau76"
        t = Time("2022-08-02", scale="tdb")
        r = [0.644307773374595, -0.7841972063903224, 3.391591147538755e-05]
        v = [0.01301350774835054, 0.01086343081039913, -2.045610960977488e-07]
        earth = State(r * u.au, v * u.au / u.day, t, frame=frame)

        r = [3.831073725981583, -0.7731565414836856, 0.1960741286643133]
        v = [-0.001734044145851301, 0.003823283255413869, 0.0005453057034865457]
        comet = State(r * u.au, v * u.au / u.day, t, frame=frame)

        # astropy and spiceypy agree within 3 arcsec, using the DE440s ephemeris
        # during the coordinate transformation does not improve the accuracy
        coords = earth.transform_to("icrs").observe(comet)
        coords_spiceypy = SkyCoord(
            358.78003306 * u.deg, 3.30890361 * u.deg, frame="icrs"
        )
        assert coords.separation(coords_spiceypy) < 2.3 * u.arcsec
        assert u.isclose(
            coords.distance, 3.192811343804468 * u.au, atol=1e5 * u.km, rtol=1e-8
        )

        # but astropy and Horizons disagree a little bit more
        coords_horizons = SkyCoord(358.779201434 * u.deg, 3.307642087 * u.deg)
        assert coords.separation(coords_horizons) < 5 * u.arcsec
        assert u.isclose(
            coords.distance, 3.19284031826821 * u.au, atol=1e5 * u.km, rtol=1e-8
        )

    def test_from_states(self):
        t = Time("2022-08-02")
        a = State(
            [1, 2, 3] * u.au,
            [4, 5, 6] * u.km / u.s,
            t,
            frame="icrs",
        )
        b = State(
            [1, 2, 3] * u.au,
            [4, 5, 6] * u.km / u.s,
            t,
            frame="heliocentriceclipticiau76",
        )
        c = State.from_states([a, b])
        b_icrs = b.transform_to("icrs")

        assert u.allclose(c[0].r, a.r)
        assert u.allclose(c[0].v, a.v)
        assert np.isclose((c[0].t - a.t).jd, 0)
        assert u.allclose(c[1].r, b_icrs.r)
        assert u.allclose(c[1].v, b_icrs.v)
        assert np.isclose((c[1].t - b.t).jd, 0)
        assert c.frame == a.frame

    def test_from_ephem(self):
        """Test Ephem to State.

        from astropy.time import Time
        from astroquery.jplhorizons import Horizons
        import spiceypy
        import sbpy.time

        t = Time("2023-12-06")
        q = Horizons("12P", id_type="designation", epochs=t.tdb.jd, location="@0")
        eph = q.vectors(refplane="earth", closest_apparition=True, aberrations="geometric")
        # au and au/dau
        rec = [float(eph[x]) for x in ("x", "y", "z", "vx", "vy", "vz")]

        # au, rad, and rad/day
        sph = spiceypy.xfmsta(rec, "rectangular", "latitudinal", " ")

        delta, ra, dec, deldot, dra, ddec = sph

        print(rec[:3])
        print(rec[3:])
        print(np.degrees((ra, dec)), "deg,", delta, "au")
        print(np.degrees((dra, ddec)), "deg/day,", deldot, "au/day")

        """

        r = [0.5849871061746752, -1.112950265888228, 1.957197574037426] * u.au
        v = (
            [-0.0008339472072142703, 0.01387010312203588, -0.006617657384488537]
            * u.au
            / u.day
        )
        t = Time("2023-12-06", scale="tdb")
        eph = Ephem.from_dict(
            {
                "x": r[0],
                "y": r[1],
                "z": r[2],
                "vx": v[0],
                "vy": v[1],
                "vz": v[2],
                "date": t,
            },
        )

        # initialize without specifying a frame
        state = State.from_ephem(eph)

        # and with the frame
        state = State.from_ephem(eph, frame="icrs")
        assert u.allclose(state.r, r)
        assert u.allclose(state.v, v)
        assert np.isclose((state.t - t).jd, 0)

        for k in ("x", "y", "z", "vx", "vy", "vz", "date"):
            incomplete = Ephem.from_table(eph.table.copy())
            del incomplete.table[k]
            with pytest.raises(ValueError):
                State.from_ephem(incomplete)

        # note, these are coordinates of 12P as seen by the solar system
        # barycenter
        eph = Ephem.from_dict(
            {
                "ra": -62.27275959 * u.deg,
                "dec": 57.28285294 * u.deg,
                "delta": 2.3262610671524557 * u.au,
                "RA*cos(Dec)_rate": 0.26043265 * u.deg / u.day,
                "Dec_rate": 0.17436216 * u.deg / u.day,
                "delta_rate": -0.012413330003023705 * u.au / u.day,
                "date": t,
            },
        )
        state = State.from_ephem(eph, "icrs")
        assert u.allclose(state.r, r)
        assert u.allclose(state.v, v)
        assert np.isclose((state.t - t).jd, 0)

        for k in ("ra", "dec", "delta", "RA*cos(Dec)_rate", "Dec_rate", "delta_rate"):
            incomplete = Ephem.from_table(eph.table.copy())
            del incomplete.table[k]
            with pytest.raises(ValueError):
                State.from_ephem(incomplete)
