# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import astropy.units as u
from .. import *


class Test_solar_spectrum:
    def test_validate_str(self):
        assert isinstance(solar_spectrum.validate('E490_2014'), Sun)

    def test_validate_Sun(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        sun = Sun.from_array(wave, fluxd, description='dummy source')
        assert isinstance(solar_spectrum.validate(sun), Sun)

    def test_validate_error(self):
        with pytest.raises(TypeError):
            solar_spectrum.validate(1)

    @pytest.mark.parametrize('name,source',
                             (('E490_2014', solar_sources.E490_2014),
                              ('E490_2014LR', solar_sources.E490_2014LR)))
    def test_set_string(self, name, source):
        with solar_spectrum.set(name):
            assert solar_spectrum.get().description == source['description']

    @pytest.mark.remote_data
    @pytest.mark.parametrize('name,source',
                             (('Kurucz1993', solar_sources.Kurucz1993),
                              ('Castelli1996', solar_sources.Castelli1996)))
    def test_set_string_remote(self, name, source):
        with solar_spectrum.set(name):
            assert solar_spectrum.get().description == source['description']

    def test_set_source(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        source = Sun.from_array(wave, fluxd, description='dummy source')
        with solar_spectrum.set(source):
            assert solar_spectrum.get().description == 'dummy source'


class Test_vega_spectrum:
    def test_validate_str(self):
        assert isinstance(vega_spectrum.validate('Bohlin2014'), Vega)

    def test_validate_Vega(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        vega = Vega.from_array(wave, fluxd, description='dummy source')
        assert isinstance(vega_spectrum.validate(vega), Vega)

    def test_validate_error(self):
        with pytest.raises(TypeError):
            vega_spectrum.validate(1)

    def test_set_string(self):
        with vega_spectrum.set('Bohlin2014'):
            assert vega_spectrum.get(
            ).description == vega_sources.Bohlin2014['description']

    def test_set_source(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        source = Vega.from_array(wave, fluxd, description='dummy source')
        with vega_spectrum.set(source):
            assert vega_spectrum.get().description == 'dummy source'
