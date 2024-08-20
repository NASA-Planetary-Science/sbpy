# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
from unittest import mock
import pytest
import numpy as np
import astropy.units as u

try:
    import synphot
except ImportError:
    pass

from ..dust import Afrho, Efrho, phase_HalleyMarcus
from ..core import CircularAperture
from ...calib import solar_fluxd
from ...units import VEGAmag, JMmag
from ...data import Ephem
from ... import photometry


Wm2um = u.W / u.m**2 / u.um


def box(center, width):
    synphot = pytest.importorskip("synphot")
    return synphot.SpectralElement(synphot.Box1D, x_0=center * u.um, width=width * u.um)


@pytest.mark.parametrize(
    "phase, value", ((0, 1.0), (15, 5.8720e-01), (14.5, 0.5959274462322928))
)
def test_phase_HalleyMarcus(phase, value):
    assert np.isclose(phase_HalleyMarcus(phase * u.deg), value)


@pytest.mark.parametrize(
    "phase, value", ((0, 1.0), (15, 5.8720e-01), (14.5, (6.0490e-01 + 5.8720e-01) / 2))
)
@mock.patch.dict(sys.modules, {"scipy": None})
def test_phase_HalleyMarcus_linear_interp(phase, value):
    assert np.isclose(phase_HalleyMarcus(phase * u.deg), value)


class TestAfrho:
    def test_init(self):
        afrho = Afrho(1000 * u.cm)
        assert afrho.value == 1000
        assert afrho.unit == u.cm

    def test_scaler_ops(self):
        afrho = Afrho(1000 * u.cm)
        afrho = afrho / 2
        assert afrho == 500 * u.cm

    def test_quantity_ops(self):
        afrho = Afrho(1000 * u.cm)
        v = afrho * 2 * u.cm
        assert v == 2000 * u.cm**2
        assert not isinstance(v, Afrho)

    @pytest.mark.parametrize(
        "wfb, fluxd0, rho, rh, delta, S, afrho0, tol",
        (
            (
                0.55 * u.um,
                6.7641725e-14 * Wm2um,
                CircularAperture(1 * u.arcsec),
                1.5,
                1.0,
                None,
                1000,
                0.001,
            ),
            (666.21 * u.THz, 5.0717 * u.mJy, 1 * u.arcsec, 1.5, 1.0, None, 1000, 0.001),
            (
                "PS1 g",
                4.68e-15 * Wm2um,
                5000 * u.km,
                4.582,
                4.042,
                1707 * Wm2um,
                1682,
                0.01,
            ),
            (
                "WFC3 F606W",
                16.98 * VEGAmag,
                5000 * u.km,
                4.582,
                4.042,
                -26.92 * VEGAmag,
                1682,
                0.01,
            ),
            (
                "PS1 r",
                11.97 * u.ABmag,
                19.2 * u.arcsec,
                1.098,
                0.164,
                -26.93 * u.ABmag,
                34.9,
                0.01,
            ),
            (
                "WFC3 F606W",
                16.98 * VEGAmag,
                5000 * u.km,
                4.582,
                4.042,
                None,
                1682,
                0.02,
            ),
            (
                "WFC3 F438W",
                17.91 * VEGAmag,
                5000 * u.km,
                4.582,
                4.042,
                None,
                1550,
                0.01,
            ),
            ("Cousins I", 7.97 * JMmag, 10000 * u.km, 1.45, 0.49, None, 3188, 0.04),
            (
                "SDSS r",
                11.97 * u.ABmag,
                19.2 * u.arcsec,
                1.098,
                0.164,
                None,
                34.9,
                0.01,
            ),
            (
                "SDSS r",
                12.23 * u.STmag,
                19.2 * u.arcsec,
                1.098,
                0.164,
                None,
                34.9,
                0.01,
            ),
        ),
    )
    def test_from_fluxd(self, wfb, fluxd0, rho, rh, delta, S, afrho0, tol):
        """Flux density to afrho conversions.

        HST/WFC3 photometry of C/2013 A1 (Siding Spring) (Li et
        al. 2014).  Uncertainty is 5%.  Li et al. (2013) quotes the
        Sun in F606W as 1730 W/m2/um and -26.93.  Confirmed with
        J.-Y. Li that 1707 W/m2/um is a better value (via effective
        stimulus formula).  Testing with revised values: 1660 → 1682
        cm, and -26.93 → -26.92 mag.

        Woodward et al. photometry of C/2007 N3 (Lulin) in I-band.
        Their magnitude has been modified from 8.49 to 7.97 according
        to their phase correction (0.03 mag/deg, phase angle 17.77
        deg).  Uncertainty is 0.06 mag.

        Li et al. (2017) DCT photometry of 252P/LINEAR in r'.  An
        additional test is performed with this observation after
        converting ABmag to STmag:

            synphot.units.convert_flux(6182, 11.97 * u.ABmag, u.STmag)

        """

        pytest.importorskip("synphot")

        eph = Ephem.from_dict(dict(rh=rh * u.au, delta=delta * u.au))

        cal = {}
        if isinstance(wfb, str):
            if S is None:
                wfb = photometry.bandpass(wfb)
            else:
                cal = {wfb: S}

        with solar_fluxd.set(cal):
            # test from_fluxd
            afrho = Afrho.from_fluxd(wfb, fluxd0, rho, eph).to("cm")
            assert np.isclose(afrho.value, afrho0, rtol=tol)

            # test to_fluxd
            fluxd = Afrho(afrho0 * u.cm).to_fluxd(wfb, rho, eph, unit=fluxd0.unit)
            k = "atol" if isinstance(fluxd, u.Magnitude) else "rtol"
            assert np.isclose(fluxd.value, fluxd0.value, **{k: tol})

    def test_source_fluxd_S_error(self):
        eph = Ephem.from_dict(dict(rh=1 * u.au, delta=1 * u.au))
        with pytest.raises(ValueError):
            with solar_fluxd.set({"J": 5 * u.um}):
                assert Afrho.from_fluxd("J", 1 * u.Jy, 1 * u.arcsec, eph)

    def test_phasecor(self):
        a0frho = Afrho(100 * u.cm)
        filt = "J"
        eph = {"rh": 1 * u.au, "delta": 1 * u.au, "phase": 100 * u.deg}
        aper = 10 * u.arcsec

        with solar_fluxd.set({filt: 1000 * Wm2um}):
            # fluxd at 0 deg phase
            fluxd0 = a0frho.to_fluxd(filt, aper, eph, phasecor=False)

            # fluxd at 100 deg phase
            fluxd = a0frho.to_fluxd(filt, aper, eph, phasecor=True)
            assert np.isclose(fluxd / fluxd0, phase_HalleyMarcus(eph["phase"]))

            # convert back to 0 deg phase
            afrho = Afrho.from_fluxd(filt, fluxd, aper, eph, phasecor=True)
            assert np.isclose(a0frho.value, afrho.value)

    def test_to_phase(self):
        afrho = Afrho(10 * u.cm).to_phase(15 * u.deg, 0 * u.deg).to("cm")
        assert np.isclose(afrho.value, 5.8720)

    def test_from_fluxd_PR125(self):
        """Regression test for PR#125: User requested Phi was ignored."""
        afrho = Afrho(100 * u.cm)
        filt = "J"
        eph = {"rh": 1 * u.au, "delta": 1 * u.au, "phase": 100 * u.deg}
        aper = 10 * u.arcsec

        def Phi(phase):
            # nonsense phase function:
            return 1 + u.Quantity(phase, "deg").value

        with solar_fluxd.set({filt: 1000 * Wm2um}):
            f0 = afrho.to_fluxd(filt, aper, eph, phasecor=False, Phi=Phi)
            f1 = afrho.to_fluxd(filt, aper, eph, phasecor=True, Phi=Phi)
            assert np.isclose(f0.value * 101, f1.value)

            a0 = Afrho.from_fluxd(filt, f0, aper, eph, phasecor=False, Phi=Phi)
            a1 = Afrho.from_fluxd(filt, f0, aper, eph, phasecor=True, Phi=Phi)
            assert np.isclose(a0.value / 101, a1.value)


class TestEfrho:
    def test_init(self):
        efrho = Efrho(1000 * u.cm)
        assert efrho.value == 1000
        assert efrho.unit == u.cm

    def test_scaler_ops(self):
        efrho = Efrho(1000 * u.cm)
        efrho = efrho / 2
        assert efrho == 500 * u.cm

    def test_quantity_ops(self):
        efrho = Efrho(1000 * u.cm)
        v = efrho * 2 * u.cm
        assert v == 2000 * u.cm**2
        assert not isinstance(v, Efrho)

    @pytest.mark.parametrize(
        "efrho0,wfb,fluxd0,rho,rh,delta,unit,B,T,tol",
        (
            (
                1000,
                10 * u.um,
                3.824064e-15 * Wm2um,
                1 * u.arcsec,
                1.5,
                1.0,
                None,
                None,
                None,
                0.001,
            ),
            (
                33.0,
                25.624 * u.THz,
                0.0060961897 * u.Jy,
                1 * u.arcsec,
                1.5,
                1.0,
                None,
                None,
                None,
                0.001,
            ),
            (
                33.0,
                11.7 * u.um,
                6.0961897 * u.mJy,
                1 * u.arcsec,
                1.5,
                1.0,
                u.mJy,
                None,
                None,
                0.001,
            ),
            (
                33.0,
                box(11.7, 0.1),
                6.0961897 * u.mJy,
                1 * u.arcsec,
                1.5,
                1.0,
                u.mJy,
                None,
                None,
                0.001,
            ),
            (
                616.1,
                box(11.7, 0.1),
                5 * VEGAmag,
                1 * u.arcsec,
                1.0,
                1.0,
                VEGAmag,
                None,
                None,
                0.001,
            ),
            (
                78750,
                11.7 * u.um,
                5 * u.ABmag,
                1 * u.arcsec,
                1.0,
                1.0,
                u.ABmag,
                None,
                None,
                0.001,
            ),
            (
                3.596e7,
                11.7 * u.um,
                5 * u.STmag,
                1 * u.arcsec,
                1.0,
                1.0,
                u.STmag,
                None,
                None,
                0.001,
            ),
            (
                31.1,
                23.7 * u.um,
                26 * u.mJy,
                18300 * u.km,
                2.52,
                2.02,
                u.mJy,
                None,
                192.3,
                0.003,
            ),
            (
                744,
                11.7 * u.um,
                6.54 * u.Jy,
                415 * u.km,
                1.11,
                0.154,
                u.Jy,
                None,
                289,
                0.003,
            ),
            (
                814,
                11.6 * u.um,
                5.94 * u.Jy,
                326 * u.km,
                1.06,
                0.15,
                u.Jy,
                None,
                290,
                0.002,
            ),
            (
                2518,
                11.7 * u.um,
                30.1 * u.Jy,
                348 * u.km,
                1.09,
                0.129,
                u.Jy,
                None,
                298,
                0.001,
            ),
            (
                2720,
                10.0 * u.um,
                13.6 * u.Jy,
                2048 * u.km,
                0.38,
                1.13,
                u.Jy,
                None,
                490,
                0.003,
            ),
            (
                52400,
                7.8 * u.um,
                1.20 * u.Jy,
                3260 * u.km,
                2.8,
                3.0,
                u.Jy,
                None,
                245,
                0.001,
            ),
            (
                53700,
                10.3 * u.um,
                2910 * u.Jy,
                7178 * u.km,
                0.34,
                0.99,
                u.Jy,
                None,
                676,
                0.001,
            ),
            (
                127e3,
                10.3 * u.um,
                1727 * u.Jy,
                12620 * u.km,
                0.59,
                1.54,
                u.Jy,
                None,
                459,
                0.002,
            ),
            (
                1.28e6,
                10.3 * u.um,
                29.4e3 * u.Jy,
                13400 * u.km,
                0.92,
                1.37,
                u.Jy,
                None,
                494.5,
                0.002,
            ),
        ),
    )
    def test_fluxd(self, efrho0, wfb, fluxd0, rho, rh, delta, unit, B, T, tol):
        """Last set of tests are compared to εfρ results of Kelley et al. 2013:

        Comet          rh (au) λ (μm)  ρ (km)  εfρ (cm)  σ (cm) Reference
        -------------- ------- ------ ------- --------- ------- ------------------
        2P/Encke          2.52   23.7  18,300      31.1         Reach et al 2007
        73P-B/S-W 3       1.11   11.7     415       744       8 Harker et al. 2011
        103P/Hartley 2    1.06   11.6     326       814      26 Meech et al. 2011
        73P-C/S-W 3       1.09   11.7     348     2,518      25 Harker et al. 2011
        2P/Encke          0.38   10.0   2,048     2,720      60 Gerhz et al. 1989
        C/1995 O1         2.8    7.8    3,260    52,400   5,200 Wooden et al. 1999
        C/1996 B2         0.34   10.3   7,180    53,700   1,500 Mason et al. 1998
        1P/Halley         0.59   10.3  12,620   127,000         Gehrz et al. 1992
        C/1995 O1         0.92   10.3  13,400 1,280,000 120,000 Mason et al. 2001

        Referred to pre-publication notes for the Kelley et al. table:
          * Improved precision of some aperture sizes.
          * Corrected aperture sizes:
            * 73P-B: 430 → 415
            * Encke: 2000 → 2048
            * C/1996 B2: 14400 → 7180
          * Added a significant figure to Encke εfρ: 31 → 31.1.

        """

        pytest.importorskip("synphot")

        eph = dict(rh=rh * u.au, delta=delta * u.au)
        efrho = Efrho.from_fluxd(wfb, fluxd0, rho, eph, Tscale=1.1, B=B, T=T).to("cm")
        assert np.isclose(efrho0, efrho.value, rtol=tol)

        fluxd = Efrho(efrho0 * u.cm).to_fluxd(
            wfb, rho, eph, unit=unit, Tscale=1.1, B=B, T=T
        )
        k = "atol" if isinstance(fluxd, u.Magnitude) else "rtol"
        assert np.isclose(fluxd0.value, fluxd.value, **{k: tol})

    def test_source_fluxd_B_error(self):
        eph = Ephem.from_dict(dict(rh=1 * u.au, delta=1 * u.au))
        with pytest.raises(ValueError):
            assert Efrho.from_fluxd(1 * u.um, 1 * u.Jy, 1 * u.arcsec, eph, B=5 * u.m)
