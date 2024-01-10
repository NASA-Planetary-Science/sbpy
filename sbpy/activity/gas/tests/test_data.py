# Licensed under a 3-clause BSD style license - see LICENSE.rst

from unittest.mock import patch
import pytest
import numpy as np
import astropy.units as u
from .. import data


class TestOHFluorescenceSA88:
    def test_spline_interpolation(self):
        pytest.importorskip("scipy")
        model = data.OHFluorescenceSA88("0-0")
        LN = model(-0.5 * u.km / u.s)
        assert np.isclose(LN.value, 1.50307895e-15)

    @patch.dict("sys.modules", {"scipy": None})
    def test_linear_interpolation(self, monkeypatch):
        model = data.OHFluorescenceSA88("0-0")
        LN = model(-0.5 * u.km / u.s)
        assert np.isclose(LN.value, 1.515e-15)

    def test_tau(self):
        model = data.OHFluorescenceSA88("0-0")
        assert np.isclose(model.tau[0].value, 2.87e5)

    def test_inversion(self):
        model = data.OHFluorescenceSA88("0-0")
        assert np.isclose(model.inversion[0], -0.304)

    def test_rdot_error(self):
        model = data.OHFluorescenceSA88("0-0")
        with pytest.raises(ValueError):
            model(-61 * u.km / u.s)

    def test_rh_error(self):
        model = data.OHFluorescenceSA88("0-0")
        with pytest.raises(ValueError):
            model({"rdot": 1 * u.km / u.s, "rh": 0.4 * u.au})
