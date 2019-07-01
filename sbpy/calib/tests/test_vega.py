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


class TestVega:
    def test___repr__(self):
        with vega_spectrum.set('Bohlin2014'):
            assert (repr(Vega.from_default()) ==
                    '<Vega: Dust-free template spectrum of Bohlin 2014>')

        vega = Vega.from_array([1, 2] * u.um, [1, 2] * u.Jy)
        assert repr(vega) == '<Vega>'

    def test_from_builtin(self):
        vega = Vega.from_builtin('Bohlin2014')
        assert vega.description == vega_sources.Bohlin2014['description']

    def test_from_builtin_unknown(self):
        with pytest.raises(UndefinedSourceError):
            Vega.from_builtin('not a vega spectrum')

    def test_from_default(self):
        with vega_spectrum.set('Bohlin2014'):
            vega = Vega.from_default()
            assert vega.description == vega_sources.Bohlin2014['description']

    def test_call_single_wavelength(self):
        with vega_spectrum.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(0.55 * u.um)
            assert np.isclose(f.value, 3.546923511485616e-08)  # W/(m2 μm)

    def test_call_single_frequency(self):
        with vega_spectrum.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(3e14 * u.Hz)
            assert np.isclose(f.value, 2129.13636259)  # Jy

    def test_show_builtin(self, capsys):
        Vega.show_builtin()
        captured = capsys.readouterr()
        for spec in vega_sources.available:
            assert spec in captured.out

    def test_observe_vega_fluxd(self):
        with vega_fluxd.set({'V': 3631 * u.Jy}):
            vega = Vega(None)
            fluxd = vega.observe('V', unit='Jy')
        assert np.isclose(fluxd.value, 3631)
