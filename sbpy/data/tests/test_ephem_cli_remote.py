# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from packaging.version import Version
from ..ephem.cli import EphemerisCLI

astroquery = pytest.importorskip("astroquery")
astroquery_version = Version(astroquery.__version__)


@pytest.mark.remote_data()
class TestEphemCLI:
    @pytest.mark.parametrize(
        "service",
        (
            "horizons",
            "mpc",
            pytest.param(
                "miriade",
                marks=pytest.mark.skipif(
                    astroquery_version < Version("0.4.12"),
                    reason="requires astroquery 0.4.12",
                ),
            ),
        ),
    )
    def test_services(self, service):
        """Spot check results from the services."""
        cli = EphemerisCLI(
            f"{service} ceres --start=2024-01-01 --step=3d --stop=2024-01-31".split()
        )
        assert len(cli.eph) == 11
        assert cli.eph["date"][0].iso == "2024-01-01 00:00:00.000"
        assert cli.eph["date"][-1].iso == "2024-01-31 00:00:00.000"

        cli = EphemerisCLI(
            f"{service} ceres --start=2024-01-01 --step=3d --number=10".split()
        )
        assert len(cli.eph) == 11
        assert cli.eph["date"][0].iso == "2024-01-01 00:00:00.000"
        assert cli.eph["date"][-1].iso == "2024-01-31 00:00:00.000"
