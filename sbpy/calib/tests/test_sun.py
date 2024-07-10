# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
import pytest
import numpy as np
import astropy.units as u
from ...units import JMmag, VEGAmag
from ...photometry import bandpass
from .. import *


class TestSun:
    def test___repr__(self):
        pytest.importorskip("synphot")
        with solar_spectrum.set("E490_2014LR"):
            assert repr(Sun.from_default()) == (
                "<Sun: E490-00a (2014) low resolution reference "
                "solar spectrum (Table 4)>"
            )

        sun = Sun.from_array([1, 2] * u.um, [1, 2] * u.Jy)
        assert repr(sun) == "<Sun>"

    def test_from_builtin(self):
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("E490_2014LR")
        assert sun.description == solar_sources.SolarSpectra.E490_2014LR["description"]

    def test_from_builtin_unknown(self):
        with pytest.raises(UndefinedSourceError):
            Sun.from_builtin("not a solar spectrum")

    def test_from_default(self):
        pytest.importorskip("synphot")
        with solar_spectrum.set("E490_2014LR"):
            sun = Sun.from_default()
            assert (
                sun.description == solar_sources.SolarSpectra.E490_2014LR["description"]
            )

    def test_call_single_wavelength(self):
        pytest.importorskip("synphot")
        with solar_spectrum.set("E490_2014"):
            sun = solar_spectrum.get()
            f = sun(0.5555 * u.um)
            assert np.isclose(f.value, 1897)

    def test_call_single_frequency(self):
        pytest.importorskip("synphot")
        with solar_spectrum.set("E490_2014"):
            sun = solar_spectrum.get()
            f = sun(3e14 * u.Hz)
            assert np.isclose(f.value, 2.49484251e14)

    def test_sun_observe_wavelength_array(self):
        pytest.importorskip("scipy")
        pytest.importorskip("synphot")
        from scipy.integrate import trapezoid

        unit = "W/(m2 um)"

        # compare Sun's rebinning with an integration over the spectrum
        sun = Sun.from_builtin("E490_2014")

        wave0 = sun.wave.to("um").value
        fluxd0 = sun.fluxd.to(unit).value

        wave = np.linspace(0.35, 0.55, 6)

        d = np.diff(wave)[0] / 2
        left_bins = wave - d
        right_bins = wave + d

        fluxd1 = np.zeros(len(wave))
        for i in range(len(wave)):
            j = (wave0 >= left_bins[i]) * (wave0 <= right_bins[i])
            fluxd1[i] = trapezoid(fluxd0[j] * wave0[j], wave0[j]) / trapezoid(
                wave0[j], wave0[j]
            )

        fluxd2 = sun.observe(wave * u.um, unit=unit).value

        assert np.allclose(fluxd1, fluxd2, rtol=0.005)

    def test_filt_units(self):
        """Colina et al. V=-26.75 mag, for zero-point flux density
        36.7e-10 ergs/s/cm2/Å.
        """
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("E490_2014")
        V = bandpass("johnson v")
        weff, fluxd = sun.observe_bandpass(V, unit="erg/(s cm2 AA)")
        assert np.isclose(weff.value, 5502, rtol=0.001)
        assert np.isclose(fluxd.value, 183.94, rtol=0.0003)

    def test_filt_vegamag(self):
        """Colina et al. V=-26.75 mag (Johnson-Morgan system)

        Not obvious we are using the same filter profile, but 0.006 mag
        agreement is good.

        """
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("E490_2014")
        V = bandpass("johnson v")
        fluxd = sun.observe(V, unit=JMmag)
        assert np.isclose(fluxd.value, -26.75, atol=0.006)

    def test_filt_abmag(self):
        """Willmer 2018 V=-26.77.

        Willmer uses Haberreiter et al. 2017 solar spectrum in the
        optical.

        """
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("E490_2014")
        V = bandpass("johnson v")
        fluxd = sun.observe(V, unit=u.ABmag)
        assert np.isclose(fluxd.value, -26.77, atol=0.007)

    def test_filt_stmag(self):
        """Willmer 2018, V=-26.76

        Willmer uses Haberreiter et al. 2017 solar spectrum in the
        optical.

        """
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("E490_2014")
        V = bandpass("johnson v")
        fluxd = sun.observe(V, unit=u.STmag)
        assert np.isclose(fluxd.value, -26.76, atol=0.003)

    def test_filt_solar_fluxd(self):
        with solar_fluxd.set({"V": -26.76 * VEGAmag}):
            sun = Sun(None)
            fluxd = sun.observe("V", unit=VEGAmag)
        assert np.isclose(fluxd.value, -26.76)

    def test_meta(self):
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("E490_2014")
        assert sun.meta is None

    @pytest.mark.remote_data
    def test_kurucz_nan_error(self):
        """sbpy#113

        Using Haberreiter et al. 2017 solar spectrum: -26.77.

        NaNs in Kurucz file should not affect this calculation.

        """
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("Kurucz1993")
        V = bandpass("johnson v")
        fluxd = sun.observe(V, unit=u.ABmag)
        assert np.isclose(fluxd.value, -26.77, atol=0.005)

    @pytest.mark.remote_data
    def test_castelli96(self):
        """Verify Castelli1996 calibration.

        According to the FITS header:

        HISTORY   Created Thu 15:49:31 16-Nov-95
        COMMENT   solar model spectrum calculated by F. Castelli.
        COMMENT   Absolute flux normalized to a V flux of 184.2 ergs/s/cm^2/A
        COMMENT   For more details see Colina, Bohlin & Castelli 1996
                  CAL/SCS-008

        2022-06-05: sbpy calculates 184.5 ergs/s/cm^2/A; agreement within 0.2%
        """
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("Castelli1996")
        V = bandpass("johnson v")
        fluxd = sun.observe(V, unit="erg/(s cm2 AA)")
        assert np.isclose(fluxd.value, 184.2, rtol=0.002)

    @pytest.mark.remote_data
    def test_calspec(self):
        """Verify CALSPEC solar model calibration."""
        pytest.importorskip("synphot")
        sun = Sun.from_builtin("calspec")
        V = bandpass("johnson v")
        fluxd = sun.observe(V, unit=VEGAmag)
        assert np.isclose(fluxd.value, -26.75, rtol=0.002)

    def test_show_builtin(self, capsys):
        Sun.show_builtin()
        captured = capsys.readouterr()
        sources = inspect.getmembers(Sun._sources, lambda v: isinstance(v, dict))
        for k, v in sources:
            assert k in captured.out
