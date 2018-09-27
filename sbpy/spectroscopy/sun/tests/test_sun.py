import pytest
import numpy as np
import astropy.units as u
from .. import *

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestSun:
    def test___repr__(self):
        with default_sun.set('E490_2014LR'):
            assert repr(Sun.from_default()
                        ) == '<Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>'

        sun = Sun.from_array([1, 2] * u.um, [1, 2] * u.Jy)
        assert repr(sun) == '<Sun>'

    def test_from_builtin(self):
        sun = Sun.from_builtin('E490_2014LR')
        assert sun.description == sources.E490_2014LR['description']

    def test_from_builtin_unknown(self):
        with pytest.raises(ValueError):
            sun = Sun.from_builtin('not a solar spectrum')

    def test_from_default(self):
        with default_sun.set('E490_2014LR'):
            sun = Sun.from_default()
            assert sun.description == sources.E490_2014LR['description']

    def test_call_single_wavelength(self):
        with default_sun.set('E490_2014'):
            sun = default_sun.get()
            f = sun(0.5555 * u.um)
            assert np.isclose(f.value, 1897)

    def test_call_single_frequency(self):
        with default_sun.set('E490_2014'):
            sun = default_sun.get()
            f = sun(3e14 * u.Hz)
            assert np.isclose(f.value, 2.49484251e+14)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_sun_wavelength_array(self):
        from scipy.integrate import trapz

        # compare Sun's rebinning with an integration over the spectrum
        sun = Sun.from_builtin('E490_2014')

        wave0 = sun.wave.to('um').value
        fluxd0 = sun.fluxd.to('W/(m2 um)').value

        wave = np.linspace(0.35, 0.55, 6)

        d = np.diff(wave)[0] / 2
        left_bins = wave - d
        right_bins = wave + d

        fluxd1 = np.zeros(len(wave))
        for i in range(len(wave)):
            j = (wave0 >= left_bins[i]) * (wave0 <= right_bins[i])
            fluxd1[i] = trapz(fluxd0[j] * wave0[j], wave0[j]) / trapz(
                wave0[j], wave0[j])

        fluxd2 = sun(wave * u.um).value

        assert np.allclose(fluxd1, fluxd2, 0.005)


class Test_default_sun:
    @pytest.mark.parametrize('name,source',
                             (('E490_2014', sources.E490_2014),
                              ('E490_2014LR', sources.E490_2014LR)))
    def test_set_string(self, name, source):
        with default_sun.set(name):
            assert default_sun.get().description == source['description']

    @pytest.mark.remote_data
    @pytest.mark.parametrize('name,source',
                             (('Kurucz1993', sources.Kurucz1993),
                              ('Castelli1996', sources.Castelli1996)))
    def test_set_string_remote(self, name, source):
        with default_sun.set(name):
            assert default_sun.get().description == source['description']

    def test_set_source(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        source = Sun.from_array(wave, fluxd, description='dummy source')
        with default_sun.set(source):
            assert default_sun.get().description == 'dummy source'
