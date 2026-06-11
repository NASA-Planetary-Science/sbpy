# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from functools import partial

import pytest
import numpy as np
from astropy.io import ascii
from astropy.time import Time
import astropy.units as u

from ..obs import Obs, Ephem, MPC
from ..core import QueryError

pytest.importorskip("astroquery")


def _data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    return os.path.join(data_dir, filename)


def mpc_get_observations_from_cache(targetid, id_type=None, **kwargs):
    """Read and return cached table.

    To generate the local data:

    from astroquery.mpc import MPC

    objects = ["9530", "2024 YR4", "100P"]
    for i, obj in enumerate(objects):
        obs = MPC.get_observations(obj)
        obs = obs[:10 + i]  # only save 10 + i observations
        obs.write(f"TestObs-from_mpc-{obj}.ecsv")

    """

    fn = _data_path(f"TestObs-from_mpc-{targetid}.ecsv")
    if not os.path.exists(fn):
        raise ValueError

    return ascii.read(fn)


@pytest.fixture
def patch_astroquery(monkeypatch):
    """Replace astroquery MPC.get_observations"""

    monkeypatch.setattr(MPC, "get_observations", mpc_get_observations_from_cache)


@pytest.fixture
def patch_ephem(monkeypatch):
    """Replace Ephem.from_* methods."""

    def ephemeris_data(service, targetid, epochs=None, location=None, **kwargs):
        n = len(epochs)
        time_column = "Date" if service == "mpc" else "epoch"
        eph = Ephem.from_dict(
            {
                "service": [service] * n,
                "target": [targetid] * n,
                "location": [location] * n,
                time_column: Time(epochs),
                "saved_epochs": Time(epochs),
                "ra": np.arange(n) * u.deg,
                "dec": -np.arange(n) * u.deg,
                "rh": np.ones(n) * u.au,
                "delta": 2 * np.ones(n) * u.au,
                "phase": 30 * np.ones(n) * u.deg,
            }
        )
        return eph

    # returns dummy data
    monkeypatch.setattr(Ephem, "from_horizons", partial(ephemeris_data, "jplhorizons"))
    monkeypatch.setattr(Ephem, "from_mpc", partial(ephemeris_data, "mpc"))
    monkeypatch.setattr(Ephem, "from_miriade", partial(ephemeris_data, "miriade"))


class TestObs:
    def test_from_mpc_asteroid_numbered(self, patch_astroquery):
        obs = Obs.from_mpc("9530")
        assert len(obs) == 10

    def test_from_mpc_asteroid_designation(self, patch_astroquery):
        obs = Obs.from_mpc("2024 YR4")
        assert len(obs) == 11

    def test_from_mpc_comet(self, patch_astroquery):
        obs = Obs.from_mpc("100P")
        assert len(obs) == 12

    def test_from_mpc_error(self, patch_astroquery):
        with pytest.raises(QueryError):
            Obs.from_mpc("C/2222 X2")

    def test_supplement_from_horizons(self, patch_astroquery, patch_ephem):
        obs = Obs.from_mpc("9530")
        obs.supplement(
            "jplhorizons", id_field="number", modify_fieldnames="obs", location="X05"
        )
        assert len(obs) == 10
        assert not any(np.isclose(obs["ra"], obs["ra_obs"]))
        assert all([row["service"] == "jplhorizons" for row in obs])
        assert all([row["location"] == "X05" for row in obs])
        assert np.allclose(obs["saved_epochs"], obs["date"])
