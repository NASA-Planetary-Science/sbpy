import pytest
import numpy as np
import astropy.units as u
import pytest
from ...units import hundred_nm
from ..core import SpectralGradient

try:
    import synphot
    from ...photometry import bandpass
except ModuleNotFoundError:
    # do nothing function so that the TestSpectralGradient class can compile
    synphot = None

    def bandpass(bp):
        pass


@pytest.mark.skipif(synphot is None, reason="requires synphot")
class TestSpectralGradient():
    def test_new(self):
        S = SpectralGradient(100 / u.um, wave=(525, 575) * u.nm)
        assert np.isclose(S.wave0.to(u.nm).value, 550)

    def test_new_wave_error(self):
        with pytest.raises(ValueError):
            SpectralGradient(100 / u.um, wave=525 * u.nm)

    @pytest.mark.parametrize('wfb, color, S0, atol', (
        ((bandpass('WFC3 F438W'), bandpass('WFC3 F606W')),
         0.09 * u.mag, 5.0 * u.percent / hundred_nm, 0.2),
        ((550, 650) * u.nm, 0.12 * u.mag, 11 * u.percent / hundred_nm, 0.5),
        ((550, 650) * u.nm, 0.11 * u.mag, 10 * u.percent / hundred_nm, 0.5),
        ((550, 650) * u.nm, -0.04 * u.mag,
         -3.2 * u.percent / hundred_nm, 1),
        ((bandpass('SDSS g'), bandpass('SDSS r')),
         0.21 * u.mag, 12 * u.percent / hundred_nm, 2),
        ((bandpass('SDSS g'), bandpass('SDSS r')),
         -0.15 * u.mag, -10 * u.percent / hundred_nm, 0.5),
         ((500, 600) * u.nm, 3, 100 * u.percent / hundred_nm, 0.01),
    ))
    def test_from_color(self, wfb, color, S0, atol):
        """Test from color to spectral gradient.

        Test 1: Li et al. 2013, ApJL, 779, L3.
        Test 2-4: Jewitt 2002, AJ, 123, 1039-1049 (assumes Δλ=100 nm
            for V-R).  Comets 2P/Encke from Luu & Jewitt 1990,
            49P/Arend-Rigaux from Millis et al. 1988, and 95P/Chiron
            from Luu 1993.
        Test 5, 6: Seccull et al. 2019, AJ, 157, 88.  sbpy computed
            13.6%/100 nm for the 0.21 mag color (atol increased to
            compensate for difference).

        """

        S = SpectralGradient.from_color(wfb, color)
        assert np.isclose(S.to(S0.unit).value, S0.value, atol=atol)

    @pytest.mark.parametrize('wfb, color0, S, atol', (
        ((bandpass('WFC3 F438W'), bandpass('WFC3 F606W')),
         0.09 * u.mag, 5.0 * u.percent / hundred_nm, 0.003),
        ((550, 650) * u.nm, 0.12 * u.mag, 11 * u.percent / hundred_nm, 0.005),
        ((550, 650) * u.nm, 0.11 * u.mag, 10 * u.percent / hundred_nm, 0.005),
        ((550, 650) * u.nm, -0.04 * u.mag,
         -3.2 * u.percent / hundred_nm, 0.01),
        ((bandpass('SDSS g'), bandpass('SDSS r')),
         0.21 * u.mag, 12 * u.percent / hundred_nm, 0.03),
        ((bandpass('SDSS g'), bandpass('SDSS r')),
         -0.15 * u.mag, -10 * u.percent / hundred_nm, 0.005),
    ))
    def test_to_color(self, wfb, color0, S, atol):
        """Inverse of from_color tests."""
        color = SpectralGradient(S, wave=wfb).to_color(wfb)
        assert np.isclose(color.to(color0.unit).value, color0.value, atol=atol)

    def test_renormalize(self):
        """Test wavelength renormalization.

        Let 650 be 10% brighter than 550, then 550 should be 0.1 / 1.1
        = 9.1% fainter than 650.

        """

        S1 = SpectralGradient(10 * u.percent / hundred_nm,
                              wave=(550, 650) * u.nm,
                              wave0=550 * u.nm)
        S2 = S1.renormalize(650 * u.nm)
        assert np.isclose(S2.to(u.percent / hundred_nm).value, 10 / 1.1)

    def test_renormalize_wave0_error(self):
        with pytest.raises(ValueError):
            SpectralGradient(10 * u.percent /
                             hundred_nm).renormalize(650 * u.nm)
