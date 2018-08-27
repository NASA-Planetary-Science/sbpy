import pytest
import numpy as np
import astropy.units as u
from .. import *

@pytest.mark.remote_data
def test_default_vega_string():
    with default_vega.set('Bohlin2014'):
        assert default_vega.get().description == sources.Bohlin2014['description']

@pytest.mark.remote_data
def test_vega_call_single():
    with default_vega.set('Bohlin2014'):
        vega = default_vega.get()
        f = vega(0.55 * u.um, unit='W/(m2 um)')
        assert np.isclose(f.value, 3.546923511485616e-08)

        f = vega(0.548 * u.um, unit='Jy')
        assert np.isclose(f.value, 3589.27)

