# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
from astropy.modeling.blackbody import blackbody_nu, blackbody_lambda
import synphot
from ..sources import BlackbodySource, SinglePointSpectrumError
from ... import bib, units, utils


class TestBlackbodySource:
    @pytest.mark.parametrize('T', (
        300, 300 * u.K
    ))
    def test_init_temperature(self, T):
        BB = BlackbodySource(T)
        assert BB.T.value == 300

    def test_init_temperature_error(self):
        with pytest.raises(TypeError):
            BlackbodySource()

    def test_repr(self):
        BB = BlackbodySource(278)
        assert repr(BB) == '<BlackbodySource: T=278.0 K>'

    @pytest.mark.parametrize('B', (blackbody_nu, blackbody_lambda))
    def test_call(self, B):
        w = np.logspace(-0.5, 3) * u.um
        f = B(w, 300 * u.K) * np.pi * u.sr
        BB = BlackbodySource(300 * u.K)
        test = BB(w, unit=f.unit).value
        assert np.allclose(test, f.value)
