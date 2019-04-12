# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from distutils.version import LooseVersion
import astropy
from astropy import units as u
from astropy.modeling import Parameter
from ..core import *
from ...data import Ephem, Phys
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
        ((get_bandpass('WFC3 F438W'), get_bandpass('WFC3 F606W')),
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
        ((get_bandpass('WFC3 F438W'), get_bandpass('WFC3 F606W')),
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
        eph = linphase.mag(pha_test)
        assert np.isclose(eph['mag'], mag_test).all()
        assert np.isclose(eph['alpha'], pha_test).all()
        assert set(eph.column_names) == {'alpha', 'mag'}

        linphase = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg, radius=300)
        pha_test = np.linspace(0, 180, 10) * u.deg
        eph = linphase.mag(pha_test)
        mag_test = np.array([5., 5.8, 6.6, 7.4, 8.2, 9., 9.8, 10.6, 11.4, 12.2])*u.mag
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(eph['mag'], mag_test).all()
            assert u.isclose(eph['alpha'], pha_test).all()
        assert set(eph.column_names) == {'alpha', 'mag'}

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
        eph = linphase.ref(pha_test)
        assert np.isclose(eph['ref'], ref_test).all()
        assert np.isclose(eph['alpha'], pha_test).all()
        assert set(eph.column_names) == {'alpha', 'ref'}
        eph_norm = linphase.ref(pha_test, normalized=0)
        assert np.isclose(eph_norm['ref'], ref_norm_test).all()
        assert np.isclose(eph_norm['alpha'], pha_test).all()
        assert set(eph_norm.column_names) == {'alpha', 'ref'}

        linphase = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg, radius=300)
        pha_test = np.linspace(0, 180, 10) * u.deg
        eph = linphase.ref(pha_test)
        ref_test = np.array(
          [1.59389035e-02, 7.62883887e-03, 3.65139185e-03, 1.74766602e-03,
           8.36485548e-04, 4.00367155e-04, 1.91627768e-04, 9.17188165e-05,
           4.38993856e-05, 2.10115670e-05])/u.sr
        ref_norm_test = np.array(
          [1., 0.47863009, 0.22908677, 0.10964782, 0.05248075,
           0.02511886, 0.01202264, 0.0057544, 0.00275423,
           0.00131826]) * u.dimensionless_unscaled
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(eph['ref'], ref_test).all()
            assert u.isclose(eph['alpha'], pha_test).all()
        assert set(eph.column_names) == {'alpha', 'ref'}
        eph_norm = linphase.ref(pha_test, normalized=0*u.deg)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(eph_norm['ref'], ref_norm_test).all()
            assert u.isclose(eph_norm['alpha'], pha_test).all()
        assert set(eph_norm.column_names) == {'alpha', 'ref'}

    def test_props(self):
        linphase = LinearPhaseFunc(5, 2.29, radius=300)
        assert np.isclose(linphase.geoalb, 0.05007354222252798)
        assert np.isclose(linphase.bondalb, 0.018404727835791654)
        assert np.isclose(linphase.phaseint, 0.3675539420399024)

    def test__distance_module(self):
        r = [0.5, 1, 1.2, 2]
        delta = [0.3, 1, 1, 2]
        m = LinearPhaseFunc(5*u.mag, 0.04*u.mag/u.deg)
        module = np.array([4.1195437, -0., -0.39590623, -3.01029996])
        assert np.isclose(m._distance_module({'r': r, 'delta': delta}), module).all()
        m = HG(5*u.mag, 0.3*u.dimensionless_unscaled)
        module = np.array([4.1195437, -0., -0.39590623, -3.01029996])
        assert np.isclose(m._distance_module({'r': r, 'delta': delta}), module).all()


class TestHG:
    def test_init(self):
        # initialize with numbers
        ceres = HG(3.34, 0.12, radius=480*u.km, M_sun=-26.74)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(ceres.radius, 480*u.km)
        assert np.isclose(ceres.M_sun, -26.74)
        assert isinstance(ceres.H, Parameter)
        assert np.isclose(ceres.H.value, 3.34)
        assert isinstance(ceres.G, Parameter)
        assert np.isclose(ceres.G.value, 0.12)
        # initialize with Quantity
        ceres = HG(3.34 * u.mag, 0.12 * u.dimensionless_unscaled)
        assert np.isclose(ceres.H.value, 3.34)
        assert ceres.H.unit == u.mag
        assert np.isclose(ceres.G.value, 0.12)
        assert ceres.G.unit == u.dimensionless_unscaled
        # test initialization from `sbpy.data.DataClass`
        phys = Phys.from_sbdb(['Ceres', 'Pallas', '12893', '3552'])
        m = HG(data = phys)
        for i, tgt in enumerate(m.meta['targetname']):
            index = phys['targetname']==tgt
            assert np.isclose(m.param_sets[:,i],
                phys[index]['H','G'].as_array().view
                (dtype=((float,float)))).all()
            assert np.isclose(m.radius[i].value,
                phys[index]['diameter'].value/2)
        # test initialization failure with `sbpy.data.DataClass` when H or G
        # is not present
        phys = Phys.from_sbdb(['12893','3552'])
        with pytest.raises(KeyError):
            m = HG(data=phys)
        # test initialization failure with `sbpy.data.DataClass` when no valid
        # parameter set is found
        phys = Phys.from_sbdb(['Ceres','12893','3552'])
        phys[0]['H'] = np.nan
        with pytest.raises(ValueError):
            m = HG(data=phys)

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
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789,
           31.57894737, 37.89473684, 44.21052632, 50.52631579, 56.84210526,
           63.15789474, 69.47368421, 75.78947368, 82.10526316, 88.42105263,
           94.73684211, 101.05263158, 107.36842105, 113.68421053,
           120.]) * u.deg
        data = np.array(
          [3.14451639, 4.06262914, 4.1154297, 4.54870242, 4.42265052,
           4.71990531, 5.1628504, 5.16098737, 5.20971821, 5.3032115,
           5.52976173, 5.64255607, 5.84536878, 6.13724017, 6.33675472,
           6.63099954, 7.2461781, 7.32734464, 8.00147425, 8.40595306]) * u.mag
        m0 = HG()
        m = m0.fit(pha, data)
        assert isinstance(m, HG)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 3.436677) & (m.H.unit == u.mag)
        assert isinstance(m.G, Parameter) & np.isclose(m.G.value, 0.1857588) & (m.G.unit == u.dimensionless_unscaled)

    def test_from_data(self):
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789, 31.57894737,
           37.89473684, 44.21052632, 50.52631579, 56.84210526, 63.15789474,
           69.47368421, 75.78947368, 82.10526316, 88.42105263, 94.73684211,
           101.05263158, 107.36842105, 113.68421053, 120.]) * u.deg
        data = np.array(
          [3.14451639, 4.06262914, 4.1154297, 4.54870242, 4.42265052,
           4.71990531, 5.1628504, 5.16098737, 5.20971821, 5.3032115,
           5.52976173, 5.64255607, 5.84536878, 6.13724017, 6.33675472,
           6.63099954, 7.2461781, 7.32734464, 8.00147425, 8.40595306]) * u.mag
        m = HG.from_data(pha, data)
        assert isinstance(m, HG)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 3.436677) & (m.H.unit == u.mag)
        assert isinstance(m.G, Parameter) & np.isclose(m.G.value, 0.1857588) & (m.G.unit == u.dimensionless_unscaled)

    def test_mag(self):
        ceres = HG(3.34, 0.12, radius=480)
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10),
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem(eph_dict)
        mag1_test = np.array(
            [6.773181346311468, 7.746327508868813, 8.297741273359549,
             8.813984838536113, 9.366879943505342, 10.024055427421063,
             10.886692329621765, 12.143261499943726, 14.18326309145893,
             18.48388800989832])
        eph1 = ceres.mag(eph_test)
        eph2 = ceres.mag(eph_dict)
        assert set(eph1.column_names) == {'alpha', 'delta', 'mag', 'r'}
        assert np.isclose(eph1['mag'], mag1_test).all()
        assert set(eph2.column_names) == {'alpha', 'delta', 'mag', 'r'}
        assert np.isclose(eph2['mag'], mag1_test).all()
        pha_test = np.linspace(0, np.pi*0.9, 10)
        mag2_test = np.array(
            [3.34, 4.313146162557345, 4.864559927048081, 5.380803492224645,
             5.9336985971938745, 6.590874081109595, 7.453510983310297,
             8.710080153632259, 10.750081745147462, 15.050706663586855])
        eph3 = ceres.mag(pha_test)
        assert np.isclose(eph3['mag'], mag2_test).all()
        assert set(eph3.column_names) == {'alpha', 'mag'}

        ceres = HG(3.34 * u.mag, 0.12 * u.dimensionless_unscaled, radius=480)
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10)*u.rad,
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem(eph_dict)
        mag1_test = np.array(
            [6.773181346311468, 7.746327508868813, 8.297741273359549,
             8.813984838536113, 9.366879943505342, 10.024055427421063,
             10.886692329621765, 12.143261499943726, 14.18326309145893,
             18.48388800989832]) * u.mag
        eph1 = ceres.mag(eph_test)
        eph2 = ceres.mag(eph_dict)
        assert set(eph1.column_names) == {'alpha', 'delta', 'mag', 'r'}
        assert set(eph2.column_names) == {'alpha', 'delta', 'mag', 'r'}
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(eph1['mag'], mag1_test).all()
            assert u.isclose(eph2['mag'], mag1_test).all()
        pha_test = np.linspace(0, np.pi*0.9, 10)
        mag2_test = np.array(
            [3.34, 4.313146162557345, 4.864559927048081, 5.380803492224645,
             5.9336985971938745, 6.590874081109595, 7.453510983310297,
             8.710080153632259, 10.750081745147462,
             15.050706663586855]) * u.mag
        eph3 = ceres.mag(pha_test)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(eph3['mag'], mag2_test).all()
        assert set(eph3.column_names) == {'alpha', 'mag'}

    def test_ref(self):
        ceres = HG(3.34, 0.12, radius=480)
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10),
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem(eph_dict)
        ref1_test = np.array(
            [0.02872225123287084, 0.0117208743870569, 0.007053334911678305,
             0.004384267859347806, 0.0026347477224558857,
             0.0014383641304522461, 0.0006498514365555728,
             0.0002042614521939071, 3.1202240400267656e-05,
             5.942043286373853e-07])
        eph1 = ceres.ref(eph_test)
        eph2 = ceres.ref(eph_dict)
        assert set(eph1.column_names) == {'alpha', 'delta', 'ref', 'r'}
        assert np.isclose(eph1['ref'], ref1_test).all()
        assert set(eph2.column_names) == {'alpha', 'delta', 'ref', 'r'}
        assert np.isclose(eph2['ref'], ref1_test).all()
        pha_test = np.linspace(0, np.pi*0.9, 10)
        ref2_test = np.array(
            [2.87222512e-02, 1.17208744e-02, 7.05333491e-03, 4.38426786e-03,
             2.63474772e-03, 1.43836413e-03, 6.49851437e-04, 2.04261452e-04,
             3.12022404e-05, 5.94204329e-07])
        ref2_norm_test = np.array(
            [1.00000000e+00, 4.08076452e-01, 2.45570407e-01, 1.52643601e-01,
             9.17319364e-02, 5.00783911e-02, 2.26253657e-02, 7.11161011e-03,
             1.08634383e-03, 2.06879441e-05])
        eph3 = ceres.ref(pha_test)
        eph4 = ceres.ref(pha_test, normalized=0*u.deg)
        assert set(eph3.column_names) == {'alpha', 'ref'}
        assert np.isclose(eph3['ref'], ref2_test).all()
        assert set(eph4.column_names) == {'alpha', 'ref'}
        assert np.isclose(eph4['ref'], ref2_norm_test).all()

        ceres = HG(3.34, 0.12, radius=480)
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10)*u.rad,
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem(eph_dict)
        ref1_test = np.array(
            [0.02872225123287084, 0.0117208743870569, 0.007053334911678305,
             0.004384267859347806, 0.0026347477224558857,
             0.0014383641304522461, 0.0006498514365555728,
             0.0002042614521939071, 3.1202240400267656e-05,
             5.942043286373853e-07])
        eph1 = ceres.ref(eph_test)
        eph2 = ceres.ref(eph_dict)
        assert set(eph1.column_names) == {'alpha', 'delta', 'ref', 'r'}
        assert set(eph2.column_names) == {'alpha', 'delta', 'ref', 'r'}
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(eph1['ref'], ref1_test*u.dimensionless_unscaled).all()
            assert u.isclose(eph2['ref'], ref1_test*u.dimensionless_unscaled).all()
        pha_test = np.linspace(0, np.pi*0.9, 10)
        ref2_test = np.array(
            [2.87222512e-02, 1.17208744e-02, 7.05333491e-03, 4.38426786e-03,
             2.63474772e-03, 1.43836413e-03, 6.49851437e-04, 2.04261452e-04,
             3.12022404e-05, 5.94204329e-07])*u.dimensionless_unscaled
        ref2_norm_test = np.array(
            [1.00000000e+00, 4.08076452e-01, 2.45570407e-01, 1.52643601e-01,
             9.17319364e-02, 5.00783911e-02, 2.26253657e-02, 7.11161011e-03,
             1.08634383e-03, 2.06879441e-05])
        eph3 = ceres.ref(pha_test)
        eph4 = ceres.ref(pha_test, normalized=0*u.deg)
        assert set(eph3.column_names) == {'alpha', 'ref'}
        assert set(eph4.column_names) == {'alpha', 'ref'}
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(eph3['ref'], ref2_test*u.dimensionless_unscaled).all()
            assert u.isclose(eph4['ref'], ref2_norm_test*u.dimensionless_unscaled).all()

    def test_g_validate(self):
        with pytest.warns(RuntimeWarning):
            m = HG(0, 1.2)


class TestHG1G2:
    def test_init(self):
        # initialization with numbers
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
        # initialization with Quantity
        themis = HG1G2(7.062 * u.mag, 0.62 * u.dimensionless_unscaled,
            0.14 * u.dimensionless_unscaled)
        assert np.isclose(themis.H.value, 7.062)
        assert themis.H.unit == u.mag
        assert np.isclose(themis.G1.value, 0.62)
        assert themis.G1.unit == u.dimensionless_unscaled
        assert np.isclose(themis.G2.value, 0.14)
        assert themis.G2.unit == u.dimensionless_unscaled
        # initialization with Phys, will cause exception because G1, G2 are
        # not generally unavailable.
        phys = Phys.from_sbdb(['Ceres'])
        with pytest.raises(KeyError):
            m = HG1G2(data=phys)

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
            [[1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00],
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
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789, 31.57894737,
           37.89473684, 44.21052632, 50.52631579, 56.84210526, 63.15789474,
           69.47368421, 75.78947368, 82.10526316, 88.42105263, 94.73684211,
           101.05263158, 107.36842105, 113.68421053, 120.]) * u.deg
        data = np.array(
          [7.14679706, 7.32220201, 7.85637226, 7.98824651, 8.2029765,
           8.27574759, 8.49437766, 9.05650671, 8.79649221, 9.33071561,
           9.24703668, 9.49069761, 9.57246629, 10.12429626, 10.14465944,
           10.51021594, 10.63215313, 11.15570421, 11.44890748,
           11.43888611]) * u.mag
        m0 = HG1G2()
        m = m0.fit(pha, data)
        assert isinstance(m, HG1G2)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 7.1167) & (m.H.unit == u.mag)
        assert isinstance(m.G1, Parameter) & np.isclose(m.G1.value, 0.63922) & (m.G1.unit == u.dimensionless_unscaled)
        assert isinstance(m.G2, Parameter) & np.isclose(m.G2.value, 0.17262569) & (m.G2.unit == u.dimensionless_unscaled)

    def test_from_data(self):
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789, 31.57894737,
           37.89473684, 44.21052632, 50.52631579, 56.84210526, 63.15789474,
           69.47368421, 75.78947368, 82.10526316, 88.42105263, 94.73684211,
           101.05263158, 107.36842105, 113.68421053, 120.]) * u.deg
        data = np.array(
          [7.14679706, 7.32220201, 7.85637226, 7.98824651, 8.2029765,
           8.27574759, 8.49437766, 9.05650671, 8.79649221, 9.33071561,
           9.24703668, 9.49069761, 9.57246629, 10.12429626, 10.14465944,
           10.51021594, 10.63215313, 11.15570421, 11.44890748,
           11.43888611]) * u.mag
        m = HG1G2.from_data(pha, data)
        assert isinstance(m, HG1G2)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 7.1167) & (m.H.unit == u.mag)
        assert isinstance(m.G1, Parameter) & np.isclose(m.G1.value, 0.63922) & (m.G1.unit == u.dimensionless_unscaled)
        assert isinstance(m.G2, Parameter) & np.isclose(m.G2.value, 0.17262569) & (m.G2.unit == u.dimensionless_unscaled)

    def test_mag(self):
        themis = HG1G2(7.063, 0.62, 0.14, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array([7.063, 8.07436233, 8.68048572,
                             9.29834638, 9.96574599, 10.72080704,
                             11.52317465, 12.15094612,
                             18.65369516, 18.65389398])
        assert np.isclose(themis.mag(pha_test)['mag'], mag_test).all()

        themis = HG1G2(7.063*u.mag, 0.62*u.dimensionless_unscaled, 0.14*u.dimensionless_unscaled, radius=100)
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        mag_test = np.array([7.063, 8.07436233, 8.68048572,
                             9.29834638, 9.96574599, 10.72080704,
                             11.52317465, 12.15094612,
                             18.65369516, 18.65389398])*u.mag
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.mag(pha_test)['mag'], mag_test).all()

    def test_ref(self):
        themis = HG1G2(7.063, 0.62, 0.14, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array(
            [2.14536381e-02, 8.45193252e-03, 4.83622683e-03,
             2.73755213e-03, 1.48048003e-03, 7.38546998e-04,
             3.52720817e-04, 1.97843827e-04, 4.95704528e-07,
             4.95613763e-07])
        assert np.isclose(themis.ref(pha_test)['ref'], ref_test).all()

        themis = HG1G2(7.063*u.mag, 0.62*u.dimensionless_unscaled, 0.14*u.dimensionless_unscaled, radius=100)
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        ref_test = np.array(
            [2.14536381e-02, 8.45193252e-03, 4.83622683e-03,
             2.73755213e-03, 1.48048003e-03, 7.38546998e-04,
             3.52720817e-04, 1.97843827e-04, 4.95704528e-07,
             4.95613763e-07])/u.sr
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.ref(pha_test)['ref'], ref_test).all()

    def test_g1g2_validator(self):
        with pytest.warns(RuntimeWarning):
            m = HG1G2(0, -0.2, 0.5)
            m = HG1G2(0, 0.5, -0.2)
            m = HG1G2(0, 0.6, 0.6)


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
             [2.74006218e-09, 1.10544425e-01, 3.31104623e-01,
              5.38630727e-01, 6.48853898e-01, 4.61007626e-01,
             -4.18315316e-01, -1.40406486e+00, 4.72646358e+00,
              4.72646358e+00]])
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
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789, 31.57894737,
           37.89473684, 44.21052632, 50.52631579, 56.84210526, 63.15789474,
           69.47368421, 75.78947368, 82.10526316, 88.42105263, 94.73684211,
           101.05263158, 107.36842105, 113.68421053, 120.]) * u.deg
        data = np.array(
          [6.95036472, 7.71609702, 8.04175457, 7.88226545, 8.28192813,
           8.50954834, 8.36880691, 8.73216696, 8.90742914, 9.05696656,
           9.20869753, 9.52578025, 9.8427691, 9.91588852, 10.3636637,
           10.26459992, 10.79316978, 10.79202241, 11.36950747, 11.61018708]) * u.mag
        m0 = HG12()
        m = m0.fit(pha, data)
        assert isinstance(m, HG12)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 7.13939) & (m.H.unit == u.mag)
        assert isinstance(m.G12, Parameter) & np.isclose(m.G12.value, 0.44872) & (m.G12.unit == u.dimensionless_unscaled)

    def test_from_data(self):
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789, 31.57894737,
           37.89473684, 44.21052632, 50.52631579, 56.84210526, 63.15789474,
           69.47368421, 75.78947368, 82.10526316, 88.42105263, 94.73684211,
           101.05263158, 107.36842105, 113.68421053, 120.]) * u.deg
        data = np.array(
          [6.95036472, 7.71609702, 8.04175457, 7.88226545, 8.28192813,
           8.50954834, 8.36880691, 8.73216696, 8.90742914, 9.05696656,
           9.20869753, 9.52578025, 9.8427691, 9.91588852, 10.3636637,
           10.26459992, 10.79316978, 10.79202241, 11.36950747,
           11.61018708]) * u.mag
        m = HG12.from_data(pha, data)
        assert isinstance(m, HG12)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 7.13939) & (m.H.unit == u.mag)
        assert isinstance(m.G12, Parameter) & np.isclose(m.G12.value, 0.44872) & (m.G12.unit == u.dimensionless_unscaled)

    def test_mag(self):
        themis = HG12(7.121, 0.68, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array(
          [7.121, 8.07252953, 8.67890827, 9.2993879, 9.96817595, 10.72086969,
          11.51208664, 12.12722017, 18.70628001, 18.70647883])
        assert np.isclose(themis.mag(pha_test)['mag'], mag_test).all()

        themis = HG12(7.121*u.mag, 0.68*u.dimensionless_unscaled, radius=100)
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        mag_test = mag_test * u.mag
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.mag(pha_test)['mag'], mag_test).all()

    def test_ref(self):
        themis = HG12(7.121, 0.68, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array(
            [2.03376585e-02, 8.46621202e-03, 4.84325842e-03,
             2.73492734e-03, 1.47717032e-03, 7.38504380e-04,
             3.56341412e-04, 2.02214774e-04, 4.72268466e-07,
             4.72181992e-07])
        assert np.isclose(themis.ref(pha_test)['ref'], ref_test).all()

        themis = HG12(7.121*u.mag, 0.68*u.dimensionless_unscaled, radius=100)
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        ref_test = ref_test / u.sr
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.ref(pha_test)['ref'], ref_test).all()

    def test_g_validator(self):
        with pytest.warns(RuntimeWarning):
            m = HG12(0, -0.71)
            m = HG12(0, 1.31)


class TestHG12_Pen16:
    def test_init(self):
        themis = HG12_Pen16(7.121, 0.68, radius=100*u.km, M_sun=-26.74)
        assert themis._unit == 'mag'
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.radius, 100*u.km)
        assert np.isclose(themis.M_sun, -26.74)
        assert isinstance(themis.H, Parameter)
        assert np.isclose(themis.H.value, 7.121)
        assert isinstance(themis.G12, Parameter)
        assert np.isclose(themis.G12.value, 0.68)

    def test__G1_G2(self):
        themis = HG12_Pen16(7.121, 0.68, radius=100*u.km, M_sun=-26.74)
        assert np.isclose(themis._G1, 0.5731968132)
        assert np.isclose(themis._G2, 0.17124272)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [7.121, 8.07252953, 8.67890827, 9.2993879, 9.96817595,
             10.72086969, 11.51208664, 12.12722017, 18.70628001, 18.70647883])
        assert np.isclose(HG12_Pen16.evaluate(pha_test, 7.121, 0.68),
                          phi_test).all()

    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [[1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
              1.00000000e+00],
             [2.41923634e-09, 8.87922520e-02, 2.87810959e-01,
              4.70632181e-01, 5.65925756e-01, 4.02776061e-01,
              -4.11888156e-01, -1.43862153e+00, 3.39292564e+00,
              3.39292564e+00]])
        assert np.isclose(np.array(HG12_Pen16.fit_deriv(
            pha_test, 7.121, 0.68)), phi_test).all()

    def test_props(self):
        themis = HG12_Pen16(7.121, 0.68, radius=100)
        assert np.isclose(themis.geoalb, 0.06389263856216909)
        assert np.isclose(themis.bondalb, 0.024306474259348763)
        assert np.isclose(themis.phaseint, 0.38042683486452)
        themis.radius = 100*u.km
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.geoalb, 0.06389264 *
                             u.dimensionless_unscaled)
            assert u.isclose(themis.bondalb, 0.02430647 *
                             u.dimensionless_unscaled)
        assert np.isclose(themis.phaseint, 0.38042683486452)

    def test_fit(self):
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789, 31.57894737,
           37.89473684, 44.21052632, 50.52631579, 56.84210526, 63.15789474,
           69.47368421, 75.78947368, 82.10526316, 88.42105263, 94.73684211,
           101.05263158, 107.36842105, 113.68421053, 120.]) * u.deg
        data = np.array(
          [7.15663893, 7.4389134, 8.00006177, 7.9044872, 8.16865497,
           8.51010016, 8.63386712, 8.65893367, 8.84895152, 9.24495642,
           9.16195702, 9.54770054, 9.60599559, 10.06129054, 10.22544773,
           10.49122575, 10.78544483, 11.12145723, 11.18055954, 11.40468613]) * u.mag
        m0 = HG12_Pen16()
        m = m0.fit(pha, data)
        assert isinstance(m, HG12_Pen16)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 7.091456) & (m.H.unit == u.mag)
        assert isinstance(m.G12, Parameter) & np.isclose(m.G12.value, 0.631243) & (m.G12.unit == u.dimensionless_unscaled)

    def test_from_data(self):
        pha = np.array(
          [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789, 31.57894737,
           37.89473684, 44.21052632, 50.52631579, 56.84210526, 63.15789474,
           69.47368421, 75.78947368, 82.10526316, 88.42105263, 94.73684211,
           101.05263158, 107.36842105, 113.68421053, 120.]) * u.deg
        data = np.array(
          [7.15663893, 7.4389134, 8.00006177, 7.9044872, 8.16865497,
           8.51010016, 8.63386712, 8.65893367, 8.84895152, 9.24495642,
           9.16195702, 9.54770054, 9.60599559, 10.06129054, 10.22544773,
           10.49122575, 10.78544483, 11.12145723, 11.18055954, 11.40468613]) * u.mag
        m = HG12_Pen16.from_data(pha, data)
        assert isinstance(m, HG12_Pen16)
        assert isinstance(m.H, Parameter) & np.isclose(m.H.value, 7.091456) & (m.H.unit == u.mag)
        assert isinstance(m.G12, Parameter) & np.isclose(m.G12.value, 0.631243) & (m.G12.unit == u.dimensionless_unscaled)

    def test_mag(self):
        themis = HG12_Pen16(7.121, 0.68, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array(
          [7.121, 8.07252953, 8.67890827, 9.2993879, 9.96817595, 10.72086969,
           11.51208664, 12.12722017, 18.70628001, 18.70647883])
        assert np.isclose(themis.mag(pha_test)['mag'], mag_test).all()

        themis = HG12_Pen16(7.121*u.mag, 0.68*u.dimensionless_unscaled, radius=100)
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        mag_test = mag_test*u.mag
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.mag(pha_test)['mag'], mag_test).all()

    def test_ref(self):
        themis = HG12_Pen16(7.121, 0.68, radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array(
            [2.03376585e-02, 8.46621202e-03, 4.84325842e-03, 2.73492734e-03,
             1.47717032e-03, 7.38504380e-04, 3.56341412e-04, 2.02214774e-04,
             4.72268466e-07, 4.72181992e-07])
        assert np.isclose(themis.ref(pha_test)['ref'], ref_test).all()

        themis = HG12_Pen16(7.121*u.mag, 0.68*u.dimensionless_unscaled, radius=100)
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        ref_test = ref_test / u.sr
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.ref(pha_test)['ref'], ref_test).all()

    def test_g_validator(self):
        with pytest.warns(RuntimeWarning):
            m = HG12(0, -0.71)
            m = HG12(0, 1.31)
