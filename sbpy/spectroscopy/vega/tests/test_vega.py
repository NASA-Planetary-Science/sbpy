import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
from .. import *


class TestVega:
    @remote_data
    def test___repr__(self):
        with default_vega.set('Bohlin2014'):
            assert repr(Vega.from_default()
                        ) == '<Vega: Spectrum of Bohlin 2014>'

        vega = Vega.from_array([1, 2] * u.um, [1, 2] * u.Jy)
        assert repr(vega) == '<Vega>'

    @remote_data
    def test_from_builtin(self):
        vega = Vega.from_builtin('Bohlin2014')
        assert vega.description == sources.Bohlin2014['description']

    def test_from_builtin_unknown(self):
        with pytest.raises(ValueError):
            Vega.from_builtin('not a vega spectrum')

    @remote_data
    def test_from_default(self):
        with default_vega.set('Bohlin2014'):
            vega = Vega.from_default()
            assert vega.description == sources.Bohlin2014['description']

    @remote_data
    def test_call_single_wavelength(self):
        with default_vega.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(0.55 * u.um)
            assert np.isclose(f.value, 3.546923511485616e-08)  # W/(m2 μm)

    @remote_data
    def test_call_single_frequency(self):
        with default_vega.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(3e14 * u.Hz)
            assert np.isclose(f.value, 2129.13636259)  # Jy


class Test_default_vega:
    @remote_data
    def test_validate_str(self):
        assert isinstance(default_vega.validate('Bohlin2014'), Vega)

    def test_validate_Vega(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        vega = Vega.from_array(wave, fluxd, description='dummy source')
        assert isinstance(default_vega.validate(vega), Vega)

    def test_validate_error(self):
        with pytest.raises(TypeError):
            default_vega.validate(1)

    @remote_data
    def test_set_string(self):
        with default_vega.set('Bohlin2014'):
            assert default_vega.get(
            ).description == sources.Bohlin2014['description']

    def test_set_source(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        source = Vega.from_array(wave, fluxd, description='dummy source')
        with default_vega.set(source):
            assert default_vega.get().description == 'dummy source'
