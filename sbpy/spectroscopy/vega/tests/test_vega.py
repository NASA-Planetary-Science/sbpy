import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
from .. import *


@remote_data
class TestVega:
    def test___repr__(self):
        with default_vega.set('Bohlin2014'):
            assert repr(Vega.from_default()
                        ) == '<Vega: Spectrum of Bohlin 2014>'

    def test_from_builtin(self):
        vega = Vega.from_builtin('Bohlin2014')
        assert vega.description == sources.Bohlin2014['description']

    def test_from_builtin_unknown(self):
        with pytest.raises(ValueError):
            sun = Sun.from_builtin('not a vega spectrum')

    def test_from_default(self):
        with default_vega.set('Bohlin2014'):
            vega = Vega.from_default()
            assert vega.description == sources.Bohlin2014['description']

    def test_call_single_wavelength(self):
        with default_vega.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(0.55 * u.um)
            assert np.isclose(f.value, 3.546923511485616e-08)  # W/(m2 μm)

    def test_call_single_frequency(self):
        with default_vega.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(3e14 * u.Hz)
            assert np.isclose(f.value, 2129.13636259)  # Jy


class Test_default_vega:
    @remote_data
    def test_set_string(self):
        with default_vega.set('Bohlin2014'):
            assert default_vega.get().description == source['description']

    def test_set_source(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        source = Vega.from_array(wave, fluxd, description='dummy source')
        with default_vega.set(source):
            assert default_vega.get().description == 'dummy source'
