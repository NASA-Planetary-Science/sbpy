# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from .. import data


class TestOHFluorescenceSA88:
    def test_linear_interpolation(self, monkeypatch):
        monkeypatch.setattr(data, 'scipy', None)
        model = data.OHFluorescenceSA88('0-0')
        LN = model(-1 * u.km / u.s)
        assert np.isclose(LN.value, 1.54e-15)

    def test_tau(self):
        model = data.OHFluorescenceSA88('0-0')
        assert np.isclose(model.tau[0].value, 2.87e5)

    def test_inversion(self):
        model = data.OHFluorescenceSA88('0-0')
        assert np.isclose(model.inversion[0], -0.304)

    def test_rdot_error(self):
        model = data.OHFluorescenceSA88('0-0')
        with pytest.raises(ValueError):
            model(-61 * u.km / u.s)

    def test_rh_error(self):
        model = data.OHFluorescenceSA88('0-0')
        with pytest.raises(ValueError):
            model({'rdot': 1 * u.km / u.s, 'rh': 0.4 * u.au})
