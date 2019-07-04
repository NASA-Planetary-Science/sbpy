import sys
import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
from ... import units as sbu
from ... import utils
from .. import core
from .. import *

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class Star(core.SpectralStandard):
    pass


class TestVega:
    def test___repr__(self):
        with vega_spectrum.set('Bohlin2014'):
            assert (repr(Vega.from_default()) ==
                    '<Vega: Dust-free template spectrum of Bohlin 2014>')

        vega = Vega.from_array([1, 2] * u.um, [1, 2] * u.Jy)
        assert repr(vega) == '<Vega>'

    @pytest.mark.parametrize('unit', ('W/(m2 um)', sbu.VEGA))
    def test_call_wavelength(self, unit):
        vega = Vega.from_default()
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(0.5 * w.value + 0.1, 'W/(m2 um)')
        s = Star.from_array(w, f)
        w = np.linspace(0.31, 0.99) * u.um

        test = (0.5 * w.value + 0.1)
        if unit == sbu.VEGA:
            test /= vega(w).value

        assert np.allclose(s(w, unit=unit).value, test, rtol=0.0005)

    def test_source_error(self, monkeypatch):
        monkeypatch.setattr(core, 'synphot', None)
        vega = Vega.from_default()
        with pytest.raises(UndefinedSourceError):
            vega.source

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

    def test_observe_vega_missing_lambda_pivot(self):
        with pytest.raises(u.UnitConversionError):
            with vega_fluxd.set({'filter1': 1 * u.Jy}):
                vega = Vega(None)
                fluxd = vega.observe(['filter1'], unit='W/(m2 um)')

    def test_observe_vega_list(self):
        with vega_fluxd.set({'filter1': 1 * u.Jy, 'filter2': 2 * u.Jy}):
            vega = Vega(None)
            fluxd = vega.observe(['filter1', 'filter2'], unit='Jy')

        assert np.allclose(fluxd.value, [1, 2])

    def test_color_index_wavelength(self):
        """Compare to Willmer 2018.

        ABmag:
        V - I = -0.013 - 0.414 = -0.427

        """
        w = [5476, 7993] * u.AA
        vega = Vega.from_default()
        lambda_eff, ci = vega.color_index(w, u.ABmag)
        assert np.allclose(lambda_eff.value, w.value)
        assert np.isclose(ci.value, -0.427, atol=0.02)

    def test_color_index_frequency(self):
        """Check that frequency is correctly coverted to wavelength."""
        w = [0.5476, 0.7993] * u.um
        f = w.to(u.Hz, u.spectral())
        vega = Vega.from_default()
        lambda_eff, ci = vega.color_index(f, u.ABmag)
        assert np.allclose(lambda_eff.value, w.value)

    def test_color_index_bandpass(self):
        """Compare to Willmer 2018."""
        bp = (utils.get_bandpass('johnson v'),
              utils.get_bandpass('cousins i'))
        vega = Vega.from_default()
        lambda_eff, ci = vega.color_index(bp, u.ABmag)
        # I bandpass seems to be very different:
        assert np.allclose(lambda_eff.to('AA').value,
                           [5476, 7993], atol=150)
        assert np.isclose(ci.value, -0.427, atol=0.02)

    def test_color_index_filter(self):
        vega = Vega.from_default()
        with vega_fluxd.set({
                'V': -0.013 * u.ABmag,
                'V(lambda eff)': 5476 * u.AA,
                'I': 0.414 * u.ABmag,
                'I(lambda eff)': 7993 * u.AA
        }):
            lambda_eff, ci = vega.color_index(('V', 'I'), u.ABmag)

        assert np.allclose(lambda_eff.value, [5476, 7993])
        assert np.isclose(ci.value, -0.427)
