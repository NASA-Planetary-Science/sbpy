# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u

from ... import units as sbu
from ... import bib
from ...photometry import bandpass
from ..core import SpectralStandard
from .. import *


class Star(SpectralStandard):
    pass


class TestSpectralStandard:
    def test_from_array(self):
        pytest.importorskip("synphot")
        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        s = Star.from_array(w, f)
        assert np.allclose(f.value, s(w, f.unit).value)

    # from_file() is tested by Sun and Vega

    def test_description(self):
        pytest.importorskip("synphot")
        s = Star(None, description="asdf")
        assert s.description == "asdf"

    def test_wave(self):
        synphot = pytest.importorskip("synphot")
        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=w, lookup_table=f
        )
        s = Star(source)
        assert np.allclose(s.wave.to(w.unit).value, w.value)

    def test_fluxd(self):
        synphot = pytest.importorskip("synphot")
        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=w, lookup_table=f
        )
        s = Star(source)
        assert np.allclose(s.fluxd.to(f.unit).value, f.value)

    def test_source(self):
        synphot = pytest.importorskip("synphot")
        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=w, lookup_table=f
        )
        s = Star(source)
        f = s.source(w)
        assert s.source == source

    # meta tested by Sun

    def test_call_frequency(self):
        pytest.importorskip("synphot")
        nu = u.Quantity(np.linspace(300, 1000), "THz")
        f = u.Quantity(0.5 * nu.value + 0.1, "Jy")
        s = Star.from_array(nu, f)
        nu = np.linspace(310, 999) * u.THz
        assert np.allclose(s(nu).value, 0.5 * nu.value + 0.1, rtol=0.002)

    def test_observe_units(self):
        pytest.importorskip("synphot")
        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)) * 0.35 / w.value, "W/(m2 um)")
        s = Star.from_array(w, f)
        w = [0.3, 0.35, 0.4] * u.um
        a = s.observe(w)
        b = s.observe(w, unit="W/(m2 Hz)")
        c = s.observe(w, unit=sbu.VEGAmag)
        d = s.observe(w, unit=u.ABmag)
        assert np.allclose(a.value, b.to(a.unit, u.spectral_density(w)).value)
        assert np.allclose(
            a.value, c.to(a.unit, sbu.spectral_density_vega(w)).value
        )
        assert c.unit == sbu.VEGAmag
        assert np.allclose(a.value, d.to(a.unit, u.spectral_density(w)).value)

    def test_observe_wavelength(self):
        pytest.importorskip("synphot")
        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)) * 0.35 / w.value, "W/(m2 um)")
        s = Star.from_array(w, f)
        w = [0.3, 0.35, 0.4] * u.um
        assert np.isclose(s.observe(w).value[1], 1.0, rtol=0.001)
        assert np.isclose(s.observe(w, interpolate=True).value[1], 1.0)

    def test_observe_frequency(self):
        pytest.importorskip("synphot")
        nu = u.Quantity(np.linspace(300, 1000), "THz")
        f = u.Quantity(np.ones(len(nu)) * nu.value / 350, "Jy")
        s = Star.from_array(nu, f)
        nu = [325, 350, 375] * u.THz
        assert np.isclose(s.observe(nu).value[1], 1.0, rtol=0.004)

    def test_observe_bandpass(self):
        synphot = pytest.importorskip("synphot")

        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        s = Star.from_array(w, f)

        bp = synphot.SpectralElement(
            synphot.Box1D, x_0=0.55 * u.um, width=0.1 * u.um
        )
        fluxd = s.observe(bp)
        assert np.allclose(fluxd.value, [1])

        bps = [
            synphot.SpectralElement(
                synphot.Box1D, x_0=0.55 * u.um, width=0.1 * u.um
            ),
            synphot.SpectralElement(
                synphot.Box1D, x_0=0.65 * u.um, width=0.1 * u.um
            ),
        ]
        fluxd = s.observe(bps, unit="W/(m2 um)")
        assert np.allclose(fluxd.value, [1, 1])

        lambda_eff, fluxd = s.observe_bandpass(bps, unit="W/(m2 um)")
        assert np.allclose(fluxd.value, [1, 1])

    def test_observe_singlepointspectrumerror(self):
        pytest.importorskip("synphot")

        from ...spectroscopy import sources

        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        s = Star.from_array(w, f)
        with pytest.raises(sources.SinglePointSpectrumError):
            s.observe(1 * u.um)
        with pytest.raises(sources.SinglePointSpectrumError):
            s.observe([1] * u.um)

    def test_observe_filter_name(self):
        pytest.importorskip("synphot")
        s = Star(None)
        s._fluxd_state = solar_fluxd
        with solar_fluxd.set({"B": 0.327 * u.ABmag}):
            fluxd = s.observe_filter_name("B", unit=u.ABmag)[-1]
            assert fluxd.value == 0.327

        with pytest.raises(FilterLookupError):
            fluxd = s.observe_filter_name("B", unit=u.ABmag)

    def test_observe_bad_wfb(self):
        synphot = pytest.importorskip("synphot")

        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=w, lookup_table=f
        )
        s = Star(source)
        with pytest.raises(TypeError):
            s.observe(np.arange(5))

    def test_bibcode(self):
        pytest.importorskip("synphot")

        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)) * 0.35 / w.value, "W/(m2 um)")
        s = Star.from_array(w, f, bibcode="asdf", description="fdsa")

        with bib.Tracking():
            s.source

        assert "asdf" in bib.show()
        bib.reset()

    def test_observe_vegamag(self):
        pytest.importorskip("synphot")

        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)), "W/(m2 um)")
        s = Star.from_array(w, f)
        V = bandpass("johnson v")
        mag = s.observe(V, unit=sbu.VEGAmag)
        # -18.60 is -2.5 * log10(3636e-11)
        assert np.isclose(mag.value, -18.60, atol=0.02)

    @pytest.mark.parametrize(
        "wfb, test, atol",
        (
            (("johnson v", "cousins i"), 0.0140 * sbu.VEGAmag, 0.004),
            ((600 * u.nm, 750 * u.nm), -0.2422 * u.ABmag, 0.001),
        ),
    )
    def test_color_index(self, wfb, test, atol):
        """Test color index.

        (1) Effective wavelengths and Vega zeropoints from Willmer
            (2018):

            -2.5 * log10(7993**3 / 5476**3) = -1.2318
            Vega = 21.1011 - 22.3469 = -1.2458
            -1.2318 - -1.2458 = 0.0140

        (2) -2.5 * log10((750**3 / 600**3)) = -0.7268
            -2.5 * log10((750**2 / 600**2)) = -0.4846
            -0.7268 - -0.4846 = -0.2422

        """

        pytest.importorskip("synphot")

        if isinstance(wfb[0], str):
            wfb = (bandpass(wfb[0]), bandpass(wfb[1]))

        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)) * w.value**-3, "W/(m2 um)")
        s = Star.from_array(w, f)
        eff_wave, ci = s.color_index(wfb, test.unit)
        assert np.isclose(ci.value, test.value, atol=atol)

    def test_color_index_typeerror(self):
        pytest.importorskip("synphot")

        w = u.Quantity(np.linspace(0.3, 1.0), "um")
        f = u.Quantity(np.ones(len(w)) * w.value**-3, "W/(m2 um)")
        s = Star.from_array(w, f)
        with pytest.raises(TypeError):
            s.color_index((None, None), u.ABmag)


class Test_solar_spectrum:
    def test_validate_str(self):
        pytest.importorskip("synphot")
        assert isinstance(solar_spectrum.validate("E490_2014"), Sun)

    def test_validate_Sun(self):
        pytest.importorskip("synphot")
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        sun = Sun.from_array(wave, fluxd, description="dummy source")
        assert isinstance(solar_spectrum.validate(sun), Sun)

    def test_validate_error(self):
        with pytest.raises(TypeError):
            solar_spectrum.validate(1)

    @pytest.mark.parametrize(
        "name,source",
        (
            ("E490_2014", solar_sources.SolarSpectra.E490_2014),
            ("E490_2014LR", solar_sources.SolarSpectra.E490_2014LR),
        ),
    )
    def test_set_string(self, name, source):
        pytest.importorskip("synphot")
        with solar_spectrum.set(name):
            assert solar_spectrum.get().description == source["description"]

    @pytest.mark.remote_data
    @pytest.mark.parametrize(
        "name,source",
        (
            ("Kurucz1993", solar_sources.SolarSpectra.Kurucz1993),
            ("Castelli1996", solar_sources.SolarSpectra.Castelli1996),
        ),
    )
    def test_set_string_remote(self, name, source):
        pytest.importorskip("synphot")
        with solar_spectrum.set(name):
            assert solar_spectrum.get().description == source["description"]

    def test_set_source(self):
        pytest.importorskip("synphot")
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        source = Sun.from_array(wave, fluxd, description="dummy source")
        with solar_spectrum.set(source):
            assert solar_spectrum.get().description == "dummy source"


class Test_vega_spectrum:
    def test_validate_str(self):
        pytest.importorskip("synphot")
        assert isinstance(vega_spectrum.validate("Bohlin2014"), Vega)

    def test_validate_Vega(self):
        pytest.importorskip("synphot")
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        vega = Vega.from_array(wave, fluxd, description="dummy source")
        assert isinstance(vega_spectrum.validate(vega), Vega)

    def test_validate_error(self):
        with pytest.raises(TypeError):
            vega_spectrum.validate(1)

    def test_set_string(self):
        pytest.importorskip("synphot")
        with vega_spectrum.set("Bohlin2014"):
            assert (
                vega_spectrum.get().description
                == vega_sources.VegaSpectra.Bohlin2014["description"]
            )

    def test_set_source(self):
        pytest.importorskip("synphot")
        wave = [1, 2] * u.um
        fluxd = [1, 2] * u.Jy
        source = Vega.from_array(wave, fluxd, description="dummy source")
        with vega_spectrum.set(source):
            assert vega_spectrum.get().description == "dummy source"


class TestSolarFluxd:
    def test_willmer2018(self):
        with solar_fluxd.set("Willmer2018"):
            filters = solar_fluxd.get()
            assert np.isclose(
                filters["PS1 r"].value, 10 ** (-0.4 * (21.1 - 26.66))
            )
            assert filters["PS1 r"].unit == "erg/(s cm2 AA)"
            assert filters["PS1 r(lambda eff)"].value == 0.6156
            assert filters["PS1 r(lambda eff)"].unit == u.um
            assert filters["PS1 r(lambda pivot)"].value == 0.6201
            assert filters["PS1 r(lambda pivot)"].unit == u.um


class TestVegaFluxd:
    def test_willmer2018(self):
        with vega_fluxd.set("Willmer2018"):
            filters = vega_fluxd.get()
            assert filters["PS1 r"].value == 2.53499e-09
            assert filters["PS1 r"].unit == "erg/(s cm2 AA)"
            assert filters["PS1 r(lambda eff)"].value == 0.6156
            assert filters["PS1 r(lambda eff)"].unit == u.um
            assert filters["PS1 r(lambda pivot)"].value == 0.6201
            assert filters["PS1 r(lambda pivot)"].unit == u.um
