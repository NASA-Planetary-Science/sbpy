import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
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
            Sun.from_builtin('not a solar spectrum')

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

    @remote_data
    def test_filt_units(self):
        """Colina et al. V=-26.75 mag, for zero-point flux density
           36.7e-10 ergs/s/cm2/Å.
        """
        sun = Sun.from_builtin('E490_2014')
        wave, fluxd = sun.filt('johnson_v', unit='erg/(s cm2 AA)')
        assert np.isclose(wave.value, 5502, rtol=0.001)
        assert np.isclose(fluxd.value, 183.94, rtol=0.0003)

    @remote_data
    def test_filt_vegamag(self):
        """Colina et al. V=-26.75 mag (Johnson system)

        Convert Johnson mag to Vega mag: V=-26.75 - 0.035 = -26.79?

        Not obvious we are using the same filter profile, but 0.02 mag
        agreement is good.

        Willmer 2018 using Haberreiter et al. 2017 solar spectrum and
        Bohlin 2014 Vega spectrum find -26.76 mag.

        """
        sun = Sun.from_builtin('E490_2014')
        wave, fluxd = sun.filt('johnson_v', unit='vegamag')
        assert np.isclose(fluxd.value, -26.79, rtol=0.001)

    @remote_data
    def test_filt_abmag(self):
        """Willmer 2018.

        Using Haberreiter et al. 2017 solar spectrum: -26.77.

        """
        sun = Sun.from_builtin('E490_2014')
        wave, fluxd = sun.filt('johnson_v', unit=u.ABmag)
        assert np.isclose(fluxd.value, -26.77, rtol=0.001)

    @remote_data
    def test_filt_stmag(self):
        """Willmer 2018.

        Using Haberreiter et al. 2017 solar spectrum: -26.76.

        """
        sun = Sun.from_builtin('E490_2014')
        wave, fluxd = sun.filt('johnson_v', unit=u.STmag)
        assert np.isclose(fluxd.value, -26.76, rtol=0.001)

    def test_meta(self):
        sun = Sun.from_builtin('E490_2014')
        assert sun.meta is None

    @remote_data
    def test_kurucz_nan_error(self):
        """Willmer 2018.

        Using Haberreiter et al. 2017 solar spectrum: -26.77.

        NaNs in Kurucz file should not affect this calulation.

        """
        sun = Sun.from_builtin('Kurucz1993')
        wave, fluxd = sun.filt('johnson_v', unit=u.ABmag)
        assert np.isclose(fluxd.value, -26.77, rtol=0.001)


class Test_default_sun:
    def test_validate_str(self):
        assert isinstance(default_sun.validate('E490_2014'), Sun)

    def test_validate_Sun(self):
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        sun = Sun.from_array(wave, fluxd, description='dummy source')
        assert isinstance(default_sun.validate(sun), Sun)

    def test_validate_error(self):
        with pytest.raises(TypeError):
            default_sun.validate(1)

    @pytest.mark.parametrize('name,source',
                             (('E490_2014', sources.E490_2014),
                              ('E490_2014LR', sources.E490_2014LR)))
    def test_set_string(self, name, source):
        with default_sun.set(name):
            assert default_sun.get().description == source['description']

    @remote_data
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
