# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from distutils.version import LooseVersion
import astropy
from astropy import units as u
from astropy.modeling import Parameter
from ..core import *
from ...data import Ephem
from ...units import hundred_nm
from ...utils import get_bandpass

req_ver = LooseVersion('3.0.2')


def test_ref2mag():
    assert np.isclose(ref2mag(0.1/np.pi, 460), 3.3208379018205285)
    if LooseVersion(astropy.__version__) >= req_ver:
        assert u.isclose(ref2mag(0.1/np.pi, 460*u.km),
                         3.3208379018205285*u.mag)


def test_mag2ref():
    assert np.isclose(mag2ref(3.32, 460), 0.03185556322261901)
    if LooseVersion(astropy.__version__) >= req_ver:
        assert u.isclose(mag2ref(3.32, 460*u.km), 0.03185556322261901/u.sr)


class TestSpectralGradient():
    def test_new(self):
        S = SpectralGradient(100 / u.um, wave=(525, 575) * u.nm)
        assert np.isclose(S.wave0.to(u.nm).value, 550)

    def test_new_wave_error(self):
        with pytest.raises(ValueError):
            SpectralGradient(100 / u.um, wave=525 * u.nm)

    @pytest.mark.parametrize('wfb, color, S0, atol', (
        ((get_bandpass('F438W'), get_bandpass('F606W')),
         0.09 * u.mag, 5.0 * u.percent / hundred_nm, 0.2),
        ((550, 650) * u.nm, 0.12 * u.mag, 11 * u.percent / hundred_nm, 0.5),
        ((550, 650) * u.nm, 0.11 * u.mag, 10 * u.percent / hundred_nm, 0.5),
        ((550, 650) * u.nm, -0.04 * u.mag,
         -3.2 * u.percent / hundred_nm, 1),
        ((get_bandpass('SDSS g'), get_bandpass('SDSS r')),
         0.21 * u.mag, 12 * u.percent / hundred_nm, 2),
        ((get_bandpass('SDSS g'), get_bandpass('SDSS r')),
         -0.15 * u.mag, -10 * u.percent / hundred_nm, 0.5),
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
        ((get_bandpass('F438W'), get_bandpass('F606W')),
         0.09 * u.mag, 5.0 * u.percent / hundred_nm, 0.003),
        ((550, 650) * u.nm, 0.12 * u.mag, 11 * u.percent / hundred_nm, 0.005),
        ((550, 650) * u.nm, 0.11 * u.mag, 10 * u.percent / hundred_nm, 0.005),
        ((550, 650) * u.nm, -0.04 * u.mag,
         -3.2 * u.percent / hundred_nm, 0.01),
        ((get_bandpass('SDSS g'), get_bandpass('SDSS r')),
         0.21 * u.mag, 12 * u.percent / hundred_nm, 0.03),
        ((get_bandpass('SDSS g'), get_bandpass('SDSS r')),
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


def test_spline():
    nodes = np.array([10, 20, 30, 40, 50])
    values = np.array([1.2, 3.5, 8.9, 15.3, 30])
    dy = np.array([0.5, 2])
    y = spline(nodes, values, dy)
    x_test = np.linspace(5, 55, 10)
    y_test = np.array([-1.3, 1.45926293, 2.73994342, 4.12536376,
                       7.223823, 10.4137113, 13.8713955, 19.75606139,
                       28.90522609, 40.])
    assert np.isclose(y(x_test), y_test).all()


class TestLinear():
    def test_init(self):
        linphase = LinearPhaseFunc(5, 2.29, radius=300)
        assert np.isclose(linphase.H.value, 5)
        assert np.isclose(linphase.S.value, 2.29)

    def test_mag(self):
        linphase = LinearPhaseFunc(5, 2.29, radius=300)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array(
            [5.0, 5.799360797413403, 6.598721594826806,
             7.398082392240209, 8.197443189653612, 8.996803987067015,
             9.796164784480418, 10.59552558189382, 11.394886379307223,
             12.194247176720626])
        assert np.isclose(linphase.mag(pha_test), mag_test).all()

    def test_ref(self):
        linphase = LinearPhaseFunc(5, 2.29, radius=300)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array(
            [1.59389035e-02, 7.63333149e-03, 3.65569373e-03,
             1.75075544e-03, 8.38457717e-04, 4.01547427e-04, 1.92305864e-04,
             9.20975780e-05, 4.41066314e-05, 2.11231932e-05])
        ref_norm_test = np.array(
            [1., 0.47891196, 0.22935666,
             0.10984165, 0.05260448, 0.02519291, 0.01206519, 0.00577816,
             0.00276723, 0.00132526])
        assert np.isclose(linphase.ref(pha_test), ref_test).all()
        assert np.isclose(linphase.ref(
            pha_test, normalized=0), ref_norm_test).all()

    def test_props(self):
        linphase = LinearPhaseFunc(5, 2.29, radius=300)
        assert np.isclose(linphase.geoalb, 0.05007354222252798)
        assert np.isclose(linphase.bondalb, 0.018404727835791654)
        assert np.isclose(linphase.phaseint, 0.3675539420399024)


class TestHG:
    def test_init(self):
        ceres = HG(3.34, 0.12, radius=480*u.km, M_sun=-26.74)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(ceres.radius, 480*u.km)
        assert np.isclose(ceres.M_sun, -26.74)
        assert isinstance(ceres.H, Parameter)
        assert np.isclose(ceres.H.value, 3.34)
        assert isinstance(ceres.G, Parameter)
        assert np.isclose(ceres.G.value, 0.12)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [3.34, 4.37957695, 4.97935011,
             5.55797517, 6.20811269, 7.03560828, 8.2295693, 10.19445176,
             14.29255427, np.inf])
        assert np.isclose(HG.evaluate(pha_test, 3.34, 0.12), phi_test).all()

    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi*0.999, 10)
        deriv_test = np.array(
            [[1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             [0., -1.33897222, -2.00745805, -2.43171634, -2.6155218,
              -2.46808816, -1.71295215, 0.06121683, 1.22716693,
              1.23379114]])
        assert np.isclose(np.array(HG.fit_deriv(
            pha_test, 3.34, 0.12)), deriv_test).all()

    def test__check_unit(self):
        ceres = HG(3.34, 0.12, radius=480*u.km, M_sun=-26.74)
        assert ceres._unit == 'mag'

    def test_props(self):
        ceres = HG(3.34, 0.12, radius=480)
        assert np.isclose(ceres.geoalb, 0.09023361346774741)
        assert np.isclose(ceres.bondalb, 0.03287666899906162)
        assert np.isclose(ceres.phaseint, 0.36435057552929395)
        ceres.radius = 480*u.km
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(ceres.geoalb,
                             0.09023361*u.dimensionless_unscaled)
            assert u.isclose(ceres.bondalb, 0.03287667 *
                             u.dimensionless_unscaled)
        assert np.isclose(ceres.phaseint, 0.36435057552929323)

    def test_fit(self):
        pass

    def test_distance_module(self):
        pass

    def test_mag(self):
        ceres = HG(3.34, 0.12, radius=480)
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10),
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem(eph_dict)
        mag1_test = np.array(
            [6.773181346311468, 7.746327508868813,
             8.297741273359549, 8.813984838536113, 9.366879943505342,
             10.024055427421063, 10.886692329621765, 12.143261499943726,
             14.18326309145893, 18.48388800989832])
        assert np.isclose(ceres.mag(eph_test), mag1_test).all()
        assert np.isclose(ceres.mag(eph_dict), mag1_test).all()
        pha_test = np.linspace(0, np.pi*0.9, 10)
        mag2_test = np.array(
            [3.34, 4.313146162557345,
             4.864559927048081, 5.380803492224645, 5.9336985971938745,
             6.590874081109595, 7.453510983310297, 8.710080153632259,
             10.750081745147462, 15.050706663586855])
        assert np.isclose(ceres.mag(pha_test), mag2_test).all()

    def test_ref(self):
        ceres = HG(3.34, 0.12, radius=480)
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10),
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem(eph_dict)
        ref1_test = np.array(
            [0.02872225123287084, 0.0117208743870569,
             0.007053334911678305, 0.004384267859347806,
             0.0026347477224558857, 0.0014383641304522461,
             0.0006498514365555728, 0.0002042614521939071,
             3.1202240400267656e-05, 5.942043286373853e-07])
        assert np.isclose(ceres.ref(eph_test), ref1_test).all()
        assert np.isclose(ceres.ref(eph_dict), ref1_test).all()
        pha_test = np.linspace(0, np.pi*0.9, 10)
        ref2_test = np.array(
            [2.87222512e-02, 1.17208744e-02, 7.05333491e-03,
             4.38426786e-03, 2.63474772e-03, 1.43836413e-03,
             6.49851437e-04, 2.04261452e-04, 3.12022404e-05,
             5.94204329e-07])
        ref2_norm_test = np.array(
            [1.00000000e+00, 4.08076452e-01,
             2.45570407e-01, 1.52643601e-01, 9.17319364e-02, 5.00783911e-02,
             2.26253657e-02, 7.11161011e-03, 1.08634383e-03, 2.06879441e-05])
        assert np.isclose(ceres.ref(pha_test), ref2_test).all()
        assert np.isclose(ceres.ref(pha_test, normalized=0),
                          ref2_norm_test).all()


class TestHG1G2:
    def test_init(self):
        themis = HG1G2(7.063, 0.62, 0.14, radius=100.*u.km, M_sun=-26.74)
        assert themis._unit == 'mag'
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.radius, 100*u.km)
        assert np.isclose(themis.M_sun, -26.74)
        assert isinstance(themis.H, Parameter)
        assert np.isclose(themis.H.value, 7.063)
        assert isinstance(themis.G1, Parameter)
        assert np.isclose(themis.G1.value, 0.62)
        assert isinstance(themis.G2, Parameter)
        assert np.isclose(themis.G2.value, 0.14)

    def test__G1_G2(self):
        themis = HG1G2(7.063, 0.62, 0.14, radius=100*u.km, M_sun=-26.74)
        assert np.isclose(themis._G1, 0.62)
        assert np.isclose(themis._G2, 0.14)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [7.063, 8.07436233, 8.68048572,
             9.29834638, 9.96574599, 10.72080704, 11.52317465, 12.15094612,
             18.65369516, 18.65389398])
        assert np.isclose(HG1G2.evaluate(
            pha_test, 7.063, 0.62, 0.14), phi_test).all()

    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi*0.999, 10)
        deriv_test = np.array(
            [[1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00],
             [2.43068538e-09, -1.21097348e+00, -1.20062175e+00,
              -1.14122248e+00, -1.10961252e+00, -1.16263148e+00,
              -1.41529341e+00, -1.70796793e+00, 0.00000000e+00,
              0.00000000e+00],
             [-6.92019713e-10, -2.07410290e+00, -2.43821941e+00,
              -2.70127334e+00, -2.84126028e+00, -2.60646207e+00,
              -1.48753067e+00, -1.91400644e-01, -7.75525861e+00,
              -7.75525861e+00]])
        assert np.isclose(np.array(HG1G2.fit_deriv(
            pha_test, 7.063, 0.62, 0.14)), deriv_test).all()

    def test_props(self):
        themis = HG1G2(7.063, 0.62, 0.14, radius=100)
        assert np.isclose(themis.geoalb, 0.06739859193616704)
        assert np.isclose(themis.bondalb, 0.02521731797010077)
        assert np.isclose(themis.phaseint, 0.374152)
        themis.radius = 100*u.km
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.geoalb, 0.06739859 *
                             u.dimensionless_unscaled)
            assert u.isclose(themis.bondalb, 0.02521732 *
                             u.dimensionless_unscaled)
        assert np.isclose(themis.phaseint, 0.374152)

    def test_fit(self):
        pass

    def test_distance_module(self):
        pass

    def test_mag(self):
        themis = HG1G2(7.063, 0.62, 0.14, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array([7.063, 8.07436233, 8.68048572,
                             9.29834638, 9.96574599, 10.72080704,
                             11.52317465, 12.15094612,
                             18.65369516, 18.65389398])
        assert np.isclose(themis.mag(pha_test), mag_test).all()

    def test_ref(self):
        themis = HG1G2(7.063, 0.62, 0.14, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array(
            [2.14536381e-02, 8.45193252e-03, 4.83622683e-03,
             2.73755213e-03, 1.48048003e-03, 7.38546998e-04,
             3.52720817e-04, 1.97843827e-04, 4.95704528e-07,
             4.95613763e-07])
        assert np.isclose(themis.ref(pha_test), ref_test).all()


class TestHG12:
    def test_init(self):
        themis = HG12(7.121, 0.68, radius=100*u.km, M_sun=-26.74)
        assert themis._unit == 'mag'
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.radius, 100*u.km)
        assert np.isclose(themis.M_sun, -26.74)
        assert isinstance(themis.H, Parameter)
        assert np.isclose(themis.H.value, 7.121)
        assert isinstance(themis.G12, Parameter)
        assert np.isclose(themis.G12.value, 0.68)

    def test__G1_G2(self):
        themis = HG12(7.121, 0.68, radius=100*u.km, M_sun=-26.74)
        assert np.isclose(themis._G1, 0.669592)
        assert np.isclose(themis._G2, 0.1407)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [7.121, 8.07252953, 8.67890827,
             9.2993879, 9.96817595, 10.72086969, 11.51208664, 12.12722017,
             18.70628001, 18.70647883])
        assert np.isclose(HG12.evaluate(pha_test, 7.121, 0.68),
                          phi_test).all()

    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [[1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00],
             [8.27405304e-10, -3.90147030e-01, -3.86830716e-01,
              -3.68566572e-01, -3.58865952e-01, -3.75398096e-01,
              -4.52786343e-01, -5.39541077e-01, 0.00000000e+00,
              0.00000000e+00]])
        assert np.isclose(np.array(HG12.fit_deriv(
            pha_test, 7.121, 0.68)), phi_test).all()

    def test_props(self):
        themis = HG12(7.121, 0.68, radius=100)
        assert np.isclose(themis.geoalb, 0.06389263856216909)
        assert np.isclose(themis.bondalb, 0.02522850358089249)
        assert np.isclose(themis.phaseint, 0.3948577512)
        themis.radius = 100*u.km
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.geoalb, 0.06389264 *
                             u.dimensionless_unscaled)
            assert u.isclose(themis.bondalb, 0.0252285 *
                             u.dimensionless_unscaled)
        assert np.isclose(themis.phaseint, 0.3948577512)

    def test_fit(self):
        pass

    def test_distance_module(self):
        pass

    def test_mag(self):
        themis = HG12(7.121, 0.68, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array([7.121, 8.07252953, 8.67890827,
                             9.2993879, 9.96817595, 10.72086969,
                             11.51208664, 12.12722017,
                             18.70628001, 18.70647883])
        assert np.isclose(themis.mag(pha_test), mag_test).all()

    def test_ref(self):
        themis = HG12(7.121, 0.68, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array(
            [2.03376585e-02, 8.46621202e-03, 4.84325842e-03,
             2.73492734e-03, 1.47717032e-03, 7.38504380e-04,
             3.56341412e-04, 2.02214774e-04, 4.72268466e-07,
             4.72181992e-07])
        assert np.isclose(themis.ref(pha_test), ref_test).all()
