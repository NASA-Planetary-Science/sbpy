# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from types import SimpleNamespace

import astropy.units as u
from astropy.time import Time

from ..ephem.cli import EphemerisCLI
from .test_ephem import patch_request  # noqa: F401


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

    def test_radec_format(self, patch_request):  # noqa: F811
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

    def test_start_stop_order(self, patch_request):  # noqa: F811
        """Test start > stop."""

        pytest.importorskip("astroquery")

        with pytest.raises(ValueError):
            EphemerisCLI(
                [
                    "horizons",
                    "C/1995 O1",
                    "--start=2025-01-01",
                    "--stop=2024-01-01",
                    "--debug",
                ]
            )

    @pytest.mark.parametrize(
        "service, target",
        [("horizons", "1 Ceres (A801 AA)"), ("mpc", "1"), ("miriade", "Ceres")],
    )
    def test_asteroid(self, service, target, patch_request):  # noqa: F811
        """Test services with a comet designation"""

        pytest.importorskip("astroquery")

        cmd = f"{service} 1 --start=2024-08-16"
        if service == "horizons":
            cmd += " --id-type=smallbody"
        cli = EphemerisCLI(cmd.split())

        assert cli.eph["date"][0].iso == "2024-08-16 00:00:00.000"
        assert cli.eph.meta["target"] == target

    @pytest.mark.parametrize(
        "service, target",
        [("horizons", "2P/Encke"), ("mpc", "2P"), ("miriade", "2P")],
    )
    def test_comet(self, service, target, patch_request):  # noqa: F811
        """Test services with a comet designation"""

        pytest.importorskip("astroquery")

        cmd = f"{service} 2P --start=2024-08-16"
        if service == "miriade":
            cmd += " --type=comet"
        cli = EphemerisCLI(cmd.split())

        assert cli.eph["date"][0].iso == "2024-08-16 00:00:00.000"
        assert cli.eph.meta["target"] == target
