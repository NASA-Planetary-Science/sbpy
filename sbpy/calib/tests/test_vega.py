import inspect
from unittest.mock import patch
import pytest
import numpy as np
import astropy.units as u
from ... import units as sbu
from ...photometry import bandpass
from .. import core
from .. import *


class Star(core.SpectralStandard):
    pass


class TestVega:
    def test___repr__(self):
        pytest.importorskip("synphot")
        with vega_spectrum.set('Bohlin2014'):
            assert (repr(Vega.from_default()) ==
                    '<Vega: Dust-free template spectrum of Bohlin 2014>')

        vega = Vega.from_array([1, 2] * u.um, [1, 2] * u.Jy)
        assert repr(vega) == '<Vega>'

    @pytest.mark.parametrize('fluxd0', (
        u.Quantity(3.44e-8, 'W/(m2 um)'), 1 * sbu.VEGA
    ))
    def test_call_wavelength(self, fluxd0):
        pytest.importorskip("synphot")
        vega = Vega.from_default()
        fluxd = vega(5557.5 * u.AA, unit=fluxd0.unit)
        assert np.isclose(fluxd.value, fluxd0.value)

    @patch.dict("sys.modules", {"synphot": None})
    def test_source_error(self):
        # Only test when synphot is not available.
        vega = Vega.from_default()
        with pytest.raises(UndefinedSourceError):
            vega.source

    def test_from_builtin(self):
        pytest.importorskip("synphot")
        vega = Vega.from_builtin('Bohlin2014')
        assert vega.description == vega_sources.VegaSpectra.Bohlin2014['description']

    def test_from_builtin_unknown(self):
        pytest.importorskip("synphot")
        with pytest.raises(UndefinedSourceError):
            Vega.from_builtin('not a vega spectrum')

    def test_from_default(self):
        pytest.importorskip("synphot")
        with vega_spectrum.set('Bohlin2014'):
            vega = Vega.from_default()
            assert vega.description == vega_sources.VegaSpectra.Bohlin2014['description']

    def test_call_single_wavelength(self):
        pytest.importorskip("synphot")
        with vega_spectrum.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(0.55 * u.um)
            assert np.isclose(f.value, 3.546923511485616e-08)  # W/(m2 μm)

    def test_call_single_frequency(self):
        pytest.importorskip("synphot")
        with vega_spectrum.set('Bohlin2014'):
            vega = Vega.from_default()
            f = vega(3e14 * u.Hz)
            assert np.isclose(f.value, 2129.13636259)  # Jy

    def test_show_builtin(self, capsys):
        pytest.importorskip("synphot")
        Vega.show_builtin()
        captured = capsys.readouterr()
        sources = inspect.getmembers(
            Vega._sources, lambda v: isinstance(v, dict))
        for k, v in sources:
            assert k in captured.out

    def test_observe_vega_fluxd(self):
        pytest.importorskip("synphot")
        with vega_fluxd.set({'V': 3631 * u.Jy}):
            vega = Vega(None)
            fluxd = vega.observe('V', unit='Jy')
        assert np.isclose(fluxd.value, 3631)

    def test_observe_vega_missing_lambda_pivot(self):
        pytest.importorskip("synphot")
        with pytest.raises(u.UnitConversionError):
            with vega_fluxd.set({'filter1': 1 * u.Jy}):
                vega = Vega(None)
                fluxd = vega.observe(['filter1'], unit='W/(m2 um)')

    def test_observe_vega_list(self):
        pytest.importorskip("synphot")
        with vega_fluxd.set({'filter1': 1 * u.Jy, 'filter2': 2 * u.Jy}):
            vega = Vega(None)
            fluxd = vega.observe(['filter1', 'filter2'], unit='Jy')

        assert np.allclose(fluxd.value, [1, 2])

    def test_color_index_wavelength(self):
        """Compare to Willmer 2018.

        ABmag:
        V - I = -0.013 - 0.414 = -0.427

        """
        pytest.importorskip("synphot")
        w = [5476, 7993] * u.AA
        vega = Vega.from_default()
        lambda_eff, ci = vega.color_index(w, u.ABmag)
        assert np.allclose(lambda_eff.value, w.value)
        assert np.isclose(ci.value, -0.427, atol=0.02)

    def test_color_index_frequency(self):
        """Check that frequency is correctly coverted to wavelength."""
        pytest.importorskip("synphot")
        w = [0.5476, 0.7993] * u.um
        f = w.to(u.Hz, u.spectral())
        vega = Vega.from_default()
        lambda_eff, ci = vega.color_index(f, u.ABmag)
        assert np.allclose(lambda_eff.value, w.value)

    def test_color_index_bandpass(self):
        pytest.importorskip("synphot")
        """Compare to Willmer 2018."""
        bp = (bandpass('johnson v'), bandpass('cousins i'))
        vega = Vega.from_default()
        lambda_eff, ci = vega.color_index(bp, u.ABmag)
        # I bandpass seems to be very different:
        assert np.allclose(lambda_eff.to('AA').value,
                           [5476, 7993], atol=150)
        assert np.isclose(ci.value, -0.427, atol=0.02)

    def test_color_index_filter(self):
        pytest.importorskip("synphot")
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
