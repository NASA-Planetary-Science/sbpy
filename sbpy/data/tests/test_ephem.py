# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import pytest
from urllib.parse import urlencode
from types import SimpleNamespace

import astropy.units as u
from astropy.time import Time

try:
    import astroquery
    from astroquery.utils.mocks import MockResponse
    from astroquery.query import BaseQuery
except ImportError:  # pragma: no cover
    astroquery = None


from ... import bib
from ..orbit import Orbit
from ..ephem import Ephem, EphemerisCLI

# retreived from Horizons on 23 Apr 2020
CERES = {
    "targetname": "1 Ceres",
    "H": u.Quantity(3.4, "mag"),
    "G": 0.12,
    "e": 0.07741102040801928,
    "q": u.Quantity(2.55375156, "au"),
    "incl": u.Quantity(10.58910839, "deg"),
    "Omega": u.Quantity(80.29081558, "deg"),
    "w": u.Quantity(73.7435117, "deg"),
    "n": u.Quantity(0.21401711, "deg / d"),
    "M": u.Quantity(154.70418799, "deg"),
    "nu": u.Quantity(158.18663933, "deg"),
    "a": u.Quantity(2.76802739, "AU"),
    "Q": u.Quantity(2.98230321, "AU"),
    "P": u.Quantity(1682.10848349, "d"),
    "epoch": Time(2458963.26397076, scale="tdb", format="jd"),
    "Tp": Time(2458240.40500675, scale="tdb", format="jd"),
}


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    return os.path.join(data_dir, filename)


def nonremote_request(self, request_type, url, **kwargs):
    """monkeypatch replacement astroquery request function"""

    params = urlencode(kwargs.get("params", {}))
    data = urlencode(kwargs.get("data", {}))

    key = (("URL", url), ("params", params), ("data", data))

    # URL to name translation, files should be saved to the ./data directory. To
    # add a new test URL-file pair, design the test as usual below and run the
    # test with the --remote-data and --capture=tee-sys parameters.  The error
    # should give enough info to formulate the key.
    #   curl "URL" -o file.txt
    #   curl "URL" --data "DATA" -o file.txt
    files = {
        (
            ("URL", "https://ssd.jpl.nasa.gov/api/horizons.api"),
            (
                "params",
                "format=text&EPHEM_TYPE=OBSERVER&QUANTITIES=%271%2C3%2C9%2C19%2C20%2C23%2C24%2C27%2C33%27&COMMAND=%22DES%3DC%2F1995+O1%3B+CAP%3C2460568%3B+NOFRAG%3B%22&SOLAR_ELONG=%220%2C180%22&LHA_CUTOFF=0&CSV_FORMAT=YES&CAL_FORMAT=BOTH&ANG_FORMAT=DEG&APPARENT=AIRLESS&REF_SYSTEM=ICRF&EXTRA_PREC=NO&CENTER=%27500%27&START_TIME=%222024-08-16+00%3A00%3A00.000%22&STOP_TIME=%222024-10-15+00%3A00%3A00.000%22&STEP_SIZE=%221d%22&SKIP_DAYLT=NO",
            ),
            ("data", ""),
        ): "TestEphemerisCLI-c1995o1-horizons.txt",
        (
            ("URL", "https://ssd.jpl.nasa.gov/api/horizons.api"),
            (
                "params",
                "format=text&EPHEM_TYPE=OBSERVER&QUANTITIES=%271%2C3%2C9%2C19%2C20%2C23%2C24%2C27%2C33%27&COMMAND=%22DES%3D2P%3B+CAP%3C2460568%3B+NOFRAG%3B%22&SOLAR_ELONG=%220%2C180%22&LHA_CUTOFF=0&CSV_FORMAT=YES&CAL_FORMAT=BOTH&ANG_FORMAT=DEG&APPARENT=AIRLESS&REF_SYSTEM=ICRF&EXTRA_PREC=NO&CENTER=%27500%27&START_TIME=%222024-08-16+00%3A00%3A00.000%22&STOP_TIME=%222024-10-15+00%3A00%3A00.000%22&STEP_SIZE=%221d%22&SKIP_DAYLT=NO",
            ),
            ("data", ""),
        ): "TestEphemerisCLI-2p-horizons.txt",
        (
            ("URL", "https://cgi.minorplanetcenter.net/cgi-bin/mpeph2.cgi"),
            ("params", ""),
            (
                "data",
                "ty=e&TextArea=2P&uto=0&igd=n&ibh=n&fp=y&adir=N&tit=&bu=&c=500&d=2024-08-16+000000&i=1&u=d&l=61&raty=a&s=t&m=h",
            ),
        ): "TestEphemerisCLI-2p-mpc.txt",
        (
            ("URL", "http://vo.imcce.fr/webservices/miriade/ephemcc_query.php"),
            (
                "params",
                "-name=2P&-type=Comet&-ep=2460538.5&-step=1.000000d&-nbd=61.0&-observer=500&-output=--jul&-tscale=UTC&-theory=INPOP&-teph=1&-tcoor=1&-rplane=1&-oscelem=ASTORB&-mime=votable",
            ),
            ("data", ""),
        ): "TestEphemerisCLI-2p-miriade.txt",
    }

    try:
        fn = files[key]
    except KeyError:  # pragma: no cover
        # printing this out helps with debugging
        print("url:", url + (("?" + params) if len(params) > 0 else ""))
        print("data:", data)
        print("files key:", key)
        raise KeyError("Request does not have corresponding mocked data")

    with open(data_path(fn), "rb") as inf:
        response = MockResponse(content=inf.read(), url=url)

    # for MPC module
    response.request = SimpleNamespace(body=data)

    return response


# use a pytest fixture to create a dummy 'requests.get' function,
# that mocks (monkeypatches) the actual 'requests.get' function:
@pytest.fixture
def patch_request(request):
    if astroquery is None:
        return None

    mp = request.getfixturevalue("monkeypatch")

    mp.setattr(BaseQuery, "_request", nonremote_request)
    return mp


class TestEphemFromOorb:
    def test_units(self):
        pytest.importorskip("pyoorb")

        orbit1 = Orbit.from_dict(CERES)
        eph1 = Ephem.from_oo(orbit1)

        orbit2 = Orbit.from_dict(
            {
                "targetname": orbit1["targetname"][0],
                "a": orbit1["a"].value[0] * u.au,
                "e": orbit1["e"][0],
                "i": orbit1["i"].value[0] * u.deg,
                "w": orbit1["w"].value[0] * u.deg,
                "Omega": orbit1["Omega"].value[0] * u.deg,
                "epoch": Time(orbit1["epoch"][0], format="jd"),
                "M": orbit1["M"].value[0] * u.deg,
                "H": orbit1["H"].value[0] * u.mag,
                "G": orbit1["G"][0],
            }
        )
        eph2 = Ephem.from_oo(orbit2)

        for k in [
            "ra",
            "dec",
            "RA*cos(Dec)_rate",
            "dec_rate",
            "alpha",
            "r",
            "delta",
            "V",
            "hlon",
            "hlat",
            "EL",
        ]:
            assert u.isclose(eph1[k], eph2[k])

    def test_basic(self):
        pytest.importorskip("pyoorb")

        orbit = Orbit.from_dict(CERES)
        oo_ephem = Ephem.from_oo(orbit, scope="basic")
        assert "dec_rate" not in oo_ephem.field_names

    def test_timescale(self):
        pytest.importorskip("pyoorb")

        orbit = Orbit.from_dict(CERES)
        epoch = Time.now()
        oo_ephem = Ephem.from_oo(orbit, epochs=epoch, scope="basic")
        assert oo_ephem["epoch"].scale == epoch.scale

    def test_bib(self):
        pytest.importorskip("pyoorb")

        with bib.Tracking():
            orbit = Orbit.from_dict(CERES)
            oo_ephem = Ephem.from_oo(orbit, scope="basic")
            assert "sbpy.data.ephem.Ephem.from_oo" in bib.show()
        bib.reset()


class TestEphemCLI:
    def test_format_epochs(self):
        """EphemerisCLI.parse_args will always have start and step.  Test
        number, stop, and neither."""

        epochs = {
            "start": Time("2024-07-23"),
            "step": u.Quantity(5, "day"),
        }
        epochs["stop"] = epochs["start"] + 10 * epochs["step"]

        # test stop
        args = SimpleNamespace(
            start=epochs["start"], step=epochs["step"], number=None, stop=epochs["stop"]
        )

        result = EphemerisCLI._format_epochs(args)
        for k in epochs.keys():
            assert result[k] == epochs[k]

        # test number
        args = SimpleNamespace(
            start=epochs["start"], step=epochs["step"], number=10, stop=None
        )

        result = EphemerisCLI._format_epochs(args)
        for k in epochs.keys():
            assert result[k] == epochs[k]

        # no number, no stop --> number = 60
        args = SimpleNamespace(
            start=epochs["start"], step=epochs["step"], number=None, stop=None
        )

        result = EphemerisCLI._format_epochs(args)
        epochs["stop"] = epochs["start"] + 60 * epochs["step"]
        for k in epochs.keys():
            assert result[k] == epochs[k]

    def test_radec_format(self, patch_request):
        """Test RA Dec formatting options."""

        pytest.importorskip("astroquery")

        cli = EphemerisCLI(
            ["horizons", "C/1995 O1", "--start=2024-08-16", "--radec=deg"]
        )

        row = str(cli.eph.table[0]).splitlines()[-1].split()
        assert row[2] == "339.94076"
        assert row[3] == "-85.76646"

        cli = EphemerisCLI(
            ["horizons", "C/1995 O1", "--start=2024-08-16", "--radec=hmsdms"]
        )
        row = str(cli.eph.table[0]).splitlines()[-1].split()
        assert row[2] == "22:39:45.78"
        assert row[3] == "-85:45:59.3"

    @pytest.mark.parametrize(
        "service, target",
        [("horizons", "2P/Encke"), ("mpc", "2P"), ("miriade", "2P")],
    )
    def test_comet(self, service, target, patch_request):
        """Test services with a comet designation"""

        pytest.importorskip("astroquery")

        cmd = f"{service} 2P --start=2024-08-16"
        if service == "miriade":
            cmd += " --type=comet"
        cli = EphemerisCLI(cmd.split())

        assert cli.eph["date"][0].iso == "2024-08-16 00:00:00.000"
        assert cli.eph.meta["target"] == target
