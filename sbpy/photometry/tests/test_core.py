# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from distutils.version import LooseVersion
import astropy
from astropy import units as u
from astropy.modeling import Parameter
from ...calib import solar_fluxd
from ..core import *
from ...data import Ephem, Phys

req_ver = LooseVersion('3.0.2')


class TestDiskIntegratedPhaseFunc():
    def test__unit(self):
        class TestClass(DiskIntegratedPhaseFunc):
            p = Parameter(default=1.)
            @staticmethod
            def evaluate(a, p):
                return a * p

        temp = TestClass()
        with pytest.raises(ValueError):
            temp._check_unit()

    def test_ref_phasefunc(self):
        solar_fluxd.set({'V': -26.77 * u.mag})
        class ExpPhase(DiskIntegratedPhaseFunc):
            _unit = 'ref'
            p = Parameter(default=0.1 / u.sr)
            nu = Parameter(default=0.1 / u.rad)
            @staticmethod
            def evaluate(a, p, nu):
                return p * np.exp(-nu * a)

        exp_phase = ExpPhase()
        pha_test = np.linspace(0, 120, 10) * u.deg
        ref_test = [0.1, 0.09769976, 0.09545244, 0.0932568, 0.09111168,
                    0.08901589, 0.08696831, 0.08496784, 0.08301337,
                    0.08110387] / u.sr
        ref = exp_phase(pha_test)
        assert np.allclose(ref.value, ref_test.value)
        assert ref.unit == ref_test.unit
        ref_n_test = [1.01760649, 0.99419913, 0.97133019, 0.94898729,
                    0.92715833, 0.90583148, 0.88499521, 0.86463822,
                    0.84474949, 0.82531824]
        ref_n = exp_phase.to_ref(pha_test, normalized=10 * u.deg)
        assert np.allclose(ref_n, ref_n_test)

        with pytest.raises(ValueError):
            mag = exp_phase.to_mag(1 * u.rad)
        with pytest.raises(ValueError):
            mag = exp_phase.to_mag(1 * u.rad, unit=u.mag)
        exp_phase.radius = 100 * u.km
        with pytest.raises(ValueError):
            mag = exp_phase.to_mag(1 * u.rad, unit=u.mag)
        exp_phase.wfb = 'V'
        mag = exp_phase.to_mag(pha_test, unit=u.mag)
        mag_test = [5.36175238, 5.38701861, 5.41228484, 5.43755106,
                5.46281729, 5.48808352, 5.51334975, 5.53861598, 5.56388221,
                5.58914844] * u.mag
        assert np.allclose(mag.value, mag_test.value)
        assert mag.unit == mag_test.unit

        eph_dict = {'alpha': pha_test,
                    'r': np.repeat(0.8 * u.au, 10),
                    'delta': np.repeat(0.5 * u.au, 10)}
        eph_test = Ephem.from_dict(eph_dict)
        ref1 = exp_phase.to_ref(eph_test)
        assert np.allclose(ref1.value, ref_test.value)


class TestLinear():
    def test_init(self):
        linphase = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg,
                radius=300 * u.km, wfb='V')
        assert np.isclose(linphase.H.value, 5)
        assert linphase.H.unit == u.mag
        assert np.isclose(linphase.S.value, 0.04)
        assert linphase.S.unit == u.mag/u.deg
        assert linphase.radius == 300 * u.km
        assert linphase.wfb == 'V'

    def test_to_mag(self):
        linphase = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg,
                radius=300 * u.km)
        pha_test = np.linspace(0, np.pi, 10) * u.rad
        mag_test = [5., 5.8, 6.6, 7.4, 8.2, 9., 9.8, 10.6, 11.4, 12.2] * u.mag
        eph = linphase.to_mag(pha_test, append_results=True)
        assert np.allclose(eph['mag'].value, mag_test.value)
        assert eph['mag'].unit == mag_test.unit
        assert np.allclose(eph['alpha'].value, pha_test.value)
        assert eph['alpha'].unit == pha_test.unit
        assert set(eph.field_names) == {'alpha', 'mag'}
        eph = linphase.to_mag(eph, append_results=True)
        assert set(eph.field_names) == {'alpha', 'mag', 'mag1'}

    def test_to_ref(self):
        linphase = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg,
                radius=300 * u.km, wfb='V')
        pha_test = np.linspace(0, 180, 10) * u.deg
        eph = linphase.to_ref(pha_test, append_results=True)
        ref_test = [1.55045242e-02, 7.42093183e-03, 3.55188129e-03,
            1.70003727e-03, 8.13688994e-04, 3.89456039e-04, 1.86405380e-04,
            8.92192241e-05, 4.27030055e-05, 2.04389434e-05] / u.sr
        ref_norm_test = np.array(
            [1., 0.47863009, 0.22908677, 0.10964782, 0.05248075,
             0.02511886, 0.01202264, 0.0057544, 0.00275423,
             0.00131826]) * u.dimensionless_unscaled
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(eph['ref'], ref_test)
            assert u.allclose(eph['alpha'], pha_test)
        assert set(eph.field_names) == {'alpha', 'ref'}
        eph_norm = linphase.to_ref(pha_test, normalized=0 * u.deg,
            append_results=True)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(eph_norm['ref'], ref_norm_test)
            assert u.allclose(eph_norm['alpha'], pha_test)
        assert set(eph_norm.field_names) == {'alpha', 'ref'}
        eph = linphase.to_ref(eph, append_results=True)
        assert set(eph.field_names) == {'alpha', 'ref', 'ref1'}
        # test exception
        linphase = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg,
                radius=300 * u.km)
        with pytest.raises(ValueError):
            linphase.to_ref(10 * u.deg)

    def test_props(self):
        linphase = LinearPhaseFunc(5 * u.mag, 2.29 * u.mag/u.rad,
                radius=300 * u.km, wfb='V')
        assert np.isclose(linphase.geomalb, 0.0487089)
        assert np.isclose(linphase.bondalb, 0.01790315)
        assert np.isclose(linphase.phaseint, 0.36755394203990327)

    def test__distance_module(self):
        r = [0.5, 1, 1.2, 2]
        delta = [0.3, 1, 1, 2]
        m = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag / u.deg)
        module_test = [0.0225, 1., 1.44, 16.]
        module = m._distance_module(Ephem.from_dict({'r': r, 'delta': delta}))
        assert np.allclose(module, module_test)

    def test_fit(self):
        pha = np.linspace(0, 60, 100) * u.deg
        mag = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg)(pha) + \
                (np.random.rand(100)*0.2-0.1) * u.mag
        from astropy.modeling.fitting import LevMarLSQFitter
        fitter = LevMarLSQFitter()
        m0 = LinearPhaseFunc(3 * u.mag, 0.02 * u.mag/u.deg)
        m = fitter(m0, pha, mag)
        assert isinstance(m, LinearPhaseFunc)

    def test_fit_deriv(self):
        assert np.allclose(LinearPhaseFunc.fit_deriv(1, 1, 2), [1, 1])


class TestHG:
    def test_init(self):
        # initialize with numbers
        ceres = HG(3.34 * u.mag, 0.12, radius=480 * u.km, wfb='V')
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(ceres.radius, 480 * u.km)
        assert isinstance(ceres.H, Parameter)
        assert np.isclose(ceres.H.value, 3.34)
        assert isinstance(ceres.G, Parameter)
        assert np.isclose(ceres.G.value, 0.12)
        assert ceres.wfb == 'V'

    @pytest.mark.remote_data
    def test_from_phys(self):
        # test initialization from `sbpy.data.DataClass`
        phys = Phys.from_sbdb('Ceres')
        m = HG.from_phys(phys)
        assert np.all(m.meta['targetname'] == phys['targetname'])
        assert np.isclose(m.H.value, phys['H'])
        assert np.isclose(m.G.value, phys['G'])
        assert np.isclose(m.radius.value, phys['diameter'].value/2)
        # test the case when target name is unknown
        phys.table.remove_column('targetname')
        m = HG.from_phys(phys)
        # test initialization failure with `sbpy.data.DataClass` when H or G
        # is not present
        phys = Phys.from_sbdb('12893')
        with pytest.raises(KeyError):
            m = HG.from_phys(phys)

    def test_to_phys(self):
        m = HG(3.34 * u.mag, 0.12)
        p = m.to_phys()
        assert isinstance(p, Phys)
        assert set(p.field_names.keys()) == {'H', 'G'}
        assert np.isclose(p['H'].value, 3.34)
        assert p['H'].unit == u.mag
        assert np.isclose(p['G'], 0.12)
        m = HG(3.34 * u.mag, 0.12, radius=480 * u.km, wfb='V',
                meta={'targetname': '1 Ceres'})
        p = m.to_phys()
        assert isinstance(p, Phys)
        assert p['targetname'] == '1 Ceres'
        assert np.isclose(p['H'].value, 3.34)
        assert p['H'].unit == u.mag
        assert np.isclose(p['G'], 0.12)
        assert np.isclose(p['diameter'].value, 960)
        assert p['diameter'].unit == u.km
        assert np.isclose(p['pv'], 0.0877745)
        assert np.isclose(p['A'], 0.03198069)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [3.34, 4.37957695, 4.97935011,
             5.55797517, 6.20811269, 7.03560828, 8.2295693, 10.19445176,
             14.29255427, np.inf])
        assert np.allclose(HG.evaluate(pha_test, 3.34, 0.12), phi_test)

    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi*0.999, 10)
        deriv_test = np.array(
            [[1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             [0., -1.33897222, -2.00745805, -2.43171634, -2.6155218,
              -2.46808816, -1.71295215, 0.06121683, 1.22716693,
              1.23379114]])
        assert np.allclose(np.array(HG.fit_deriv(pha_test, 3.34, 0.12)),
            deriv_test)
        assert np.allclose(HG.fit_deriv(1, 3.4, 0.2), [1.0, -2.031224359464])

    def test__check_unit(self):
        ceres = HG(3.34 * u.mag, 0.12)
        assert ceres._unit == 'mag'

    def test_props(self):
        ceres = HG(3.34 * u.mag, 0.12, radius=480 * u.km, wfb='V')
        assert np.isclose(ceres.geomalb, 0.0877745)
        assert np.isclose(ceres.bondalb, 0.03198069)
        assert np.isclose(ceres.phaseint, 0.3643505755292945)

    def test_from_obs(self):
        pha = [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789,
             31.57894737, 37.89473684, 44.21052632, 50.52631579, 56.84210526,
             63.15789474, 69.47368421, 75.78947368, 82.10526316, 88.42105263,
             94.73684211, 101.05263158, 107.36842105, 113.68421053,
             120.] * u.deg
        data = [3.14451639, 4.06262914, 4.1154297, 4.54870242, 4.42265052,
             4.71990531, 5.1628504, 5.16098737, 5.20971821, 5.3032115,
             5.52976173, 5.64255607, 5.84536878, 6.13724017, 6.33675472,
             6.63099954, 7.2461781, 7.32734464, 8.00147425,
             8.40595306] * u.mag
        from astropy.modeling.fitting import LevMarLSQFitter
        fitter = LevMarLSQFitter()
        # test fit with one column
        m = HG.from_obs({'alpha': pha, 'mag': data}, fitter)
        assert isinstance(m, HG)
        assert isinstance(m.H, Parameter) & np.isclose(
            m.H.value, 3.436677) & (m.H.unit == u.mag)
        assert isinstance(m.G, Parameter) & np.isclose(
            m.G.value, 0.1857588) & (m.G.unit == u.dimensionless_unscaled)
        # test fit with one column and `init` parameters
        m = HG.from_obs({'alpha': pha, 'mag': data}, fitter, init=[3, 0.1])
        assert isinstance(m, HG)
        assert isinstance(m.H, Parameter) & np.isclose(
            m.H.value, 3.4366849) & (m.H.unit == u.mag)
        assert isinstance(m.G, Parameter) & np.isclose(
            m.G.value, 0.18576319) & (m.G.unit == u.dimensionless_unscaled)
        # test fit with more than one column
        m = HG.from_obs({'alpha': pha, 'mag': data, 'mag1': data,
            'mag2': data}, fitter, fields=['mag', 'mag1', 'mag2'])
        assert isinstance(m, HG)
        assert isinstance(m.H, Parameter) & np.allclose(
            m.H.value, [3.436677]*3) & (m.H.unit == u.mag)
        assert isinstance(m.G, Parameter) & np.allclose(
            m.G.value, [0.1857588]*3) & (m.G.unit == u.dimensionless_unscaled)
        assert 'fields' in m.meta
        assert m.meta['fields'] == ['mag', 'mag1', 'mag2']
        # test fit with more than one column with `init` parameters
        m = HG.from_obs({'alpha': pha, 'mag': data, 'mag1': data,
            'mag2': data}, fitter, fields=['mag', 'mag1', 'mag2'],
            init=[[3., 3., 3.], [0.1, 0.1, 0.1]])
        assert isinstance(m, HG)
        assert isinstance(m.H, Parameter) & np.allclose(
            m.H.value, [3.4366849]*3) & (m.H.unit == u.mag)
        assert isinstance(m.G, Parameter) & np.allclose(
            m.G.value, [0.18576319]*3) & (m.G.unit == u.dimensionless_unscaled)
        assert 'fields' in m.meta
        assert m.meta['fields'] == ['mag', 'mag1', 'mag2']

    def test_to_mag(self):
        ceres = HG(3.34 * u.mag, 0.12, radius=480 * u.km, wfb='V')
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10) * u.rad,
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem.from_dict(eph_dict)
        mag1_test = np.array(
            [6.773181346311468, 7.746327508868813, 8.297741273359549,
             8.813984838536113, 9.366879943505342, 10.024055427421063,
             10.886692329621765, 12.143261499943726, 14.18326309145893,
             18.48388800989832]) * u.mag
        eph1 = ceres.to_mag(eph_test, append_results=True)
        assert set(eph1.field_names) == {'alpha', 'delta', 'mag', 'r'}
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(eph1['mag'], mag1_test)
        pha_test = np.linspace(0, np.pi*0.9, 10) * u.rad
        mag2_test = np.array(
            [3.34, 4.313146162557345, 4.864559927048081, 5.380803492224645,
             5.9336985971938745, 6.590874081109595, 7.453510983310297,
             8.710080153632259, 10.750081745147462,
             15.050706663586855]) * u.mag
        eph3 = ceres.to_mag(pha_test, append_results=True)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(eph3['mag'], mag2_test)
        assert set(eph3.field_names) == {'alpha', 'mag'}
        mag4 = ceres.to_mag(pha_test)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(mag4, mag2_test)

    def test_to_ref(self):
        ceres = HG(3.34 * u.mag, 0.12, radius=480 * u.km, wfb='V')
        eph_dict = {'alpha': np.linspace(0, np.pi*0.9, 10)*u.rad,
                    'r': np.repeat(2.7*u.au, 10),
                    'delta': np.repeat(1.8*u.au, 10)}
        eph_test = Ephem.from_dict(eph_dict)
        ref1_test = [2.79394901e-02, 1.14014480e-02, 6.86111195e-03,
            4.26478439e-03, 2.56294353e-03, 1.39916471e-03, 6.32141181e-04,
            1.98694761e-04, 3.03518927e-05, 5.78010611e-07] / u.sr
        eph1 = ceres.to_ref(eph_test, append_results=True)
        assert set(eph1.field_names) == {'alpha', 'delta', 'ref', 'r'}
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(eph1['ref'], ref1_test)
        pha_test = np.linspace(0, np.pi*0.9, 10) * u.rad
        ref2_norm_test = np.array(
            [1.00000000e+00, 4.08076452e-01, 2.45570407e-01, 1.52643601e-01,
             9.17319364e-02, 5.00783911e-02, 2.26253657e-02, 7.11161011e-03,
             1.08634383e-03, 2.06879441e-05])
        eph3 = ceres.to_ref(pha_test, append_results=True)
        eph4 = ceres.to_ref(pha_test, normalized=0*u.deg, append_results=True)
        assert set(eph3.field_names) == {'alpha', 'ref'}
        assert set(eph4.field_names) == {'alpha', 'ref'}
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(eph3['ref'], ref1_test)
            assert u.allclose(eph4['ref'], ref2_norm_test)
        ref5 = ceres.to_ref(pha_test)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(ref5, ref1_test)

    def test_g_validate(self):
        with pytest.warns(NonmonotonicPhaseFunctionWarning):
            m = HG(0, 1.2)

    def test_hgphi_exception(self):
        with pytest.raises(ValueError):
            tmp = HG._hgphi(0, 0)
        with pytest.raises(ValueError):
            tmp = HG._hgphi(0, 3)


class TestHG1G2:
    def test_init(self):
        # initialization with numbers
        themis = HG1G2(7.063 * u.mag, 0.62, 0.14, radius=100 * u.km, wfb='V')
        assert themis._unit == 'mag'
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.radius, 100 * u.km)
        assert isinstance(themis.H, Parameter)
        assert np.isclose(themis.H.value, 7.063)
        assert isinstance(themis.G1, Parameter)
        assert np.isclose(themis.G1.value, 0.62)
        assert isinstance(themis.G2, Parameter)
        assert np.isclose(themis.G2.value, 0.14)
        assert themis.wfb == 'V'

    @pytest.mark.remote_data
    def test_from_phys(self):
        # initialization with Phys, will cause exception because G1, G2 are
        # not generally unavailable.
        phys = Phys.from_sbdb(['Ceres'])
        with pytest.raises(KeyError):
            m = HG1G2.from_phys(phys)

    def test__G1_G2(self):
        themis = HG1G2(7.063 * u.mag, 0.62, 0.14)
        assert np.isclose(themis._G1, 0.62)
        assert np.isclose(themis._G2, 0.14)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [7.063, 8.07436233, 8.68048572,
             9.29834638, 9.96574599, 10.72080704, 11.52317465, 12.15094612,
             18.65369516, 18.65389398])
        assert np.allclose(HG1G2.evaluate(
            pha_test, 7.063, 0.62, 0.14), phi_test)

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
        assert np.allclose(np.array(HG1G2.fit_deriv(
            pha_test, 7.063, 0.62, 0.14)), deriv_test)
        assert np.allclose(HG1G2.fit_deriv(0.2, 3.4, 0.62, 0.14),
                [1.0, -1.1563984700303085, -1.6666940099848913])

    def test_props(self):
        themis = HG1G2(7.063 * u.mag, 0.62, 0.14, radius=100 * u.km, wfb='V')
        assert np.isclose(themis.geomalb, 0.06556179)
        assert np.isclose(themis.bondalb, 0.02453008)
        assert np.isclose(themis.phaseint, 0.374152)

    def test_from_obs(self):
        pha = [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789,
               31.57894737, 37.89473684, 44.21052632, 50.52631579, 56.84210526,
               63.15789474, 69.47368421, 75.78947368, 82.10526316, 88.42105263,
               94.73684211, 101.05263158, 107.36842105, 113.68421053,
               120.] * u.deg
        data = [7.14679706, 7.32220201, 7.85637226, 7.98824651, 8.2029765,
                8.27574759, 8.49437766, 9.05650671, 8.79649221, 9.33071561,
                9.24703668, 9.49069761, 9.57246629, 10.12429626, 10.14465944,
                10.51021594, 10.63215313, 11.15570421, 11.44890748,
                11.43888611] * u.mag
        from astropy.modeling.fitting import LevMarLSQFitter
        fitter = LevMarLSQFitter()
        m = HG1G2.from_obs({'alpha': pha, 'mag': data}, fitter)
        assert isinstance(m, HG1G2)
        assert isinstance(m.H, Parameter) & np.isclose(
            m.H.value, 7.1167) & (m.H.unit == u.mag)
        assert isinstance(m.G1, Parameter) & np.isclose(
            m.G1.value, 0.63922) & (m.G1.unit == u.dimensionless_unscaled)
        assert isinstance(m.G2, Parameter) & np.isclose(
            m.G2.value, 0.17262569) & (m.G2.unit == u.dimensionless_unscaled)

    def test_to_mag(self):
        themis = HG1G2(7.063 * u.mag, 0.62, 0.14, radius=100 * u.km, wfb='V')
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        mag_test = [7.063, 8.07436233, 8.68048572, 9.29834638, 9.96574599,
            10.72080704, 11.52317465, 12.15094612, 18.65369516,
            18.65389398] * u.mag
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(themis.to_mag(pha_test), mag_test)

    def test_to_ref(self):
        themis = HG1G2(7.063 * u.mag, 0.62, 0.14, radius=100 * u.km, wfb='V')
        pha_test = np.linspace(0, np.pi, 10)*u.rad
        ref_test = [2.08689669e-02, 8.22159390e-03, 4.70442623e-03,
                2.66294623e-03, 1.44013284e-03, 7.18419542e-04, 3.43108196e-04,
                1.92452033e-04, 4.82195204e-07, 4.82106912e-07] / u.sr
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(themis.to_ref(pha_test), ref_test)

    def test_g1g2_validator(self):
        with pytest.warns(NonmonotonicPhaseFunctionWarning):
            m = HG1G2(0, -0.2, 0.5)
            m = HG1G2(0, 0.5, -0.2)
            m = HG1G2(0, 0.6, 0.6)


class TestHG12:
    def test_init(self):
        themis = HG12(7.121 * u.mag, 0.68, radius=100 * u.km, wfb='V')
        assert themis._unit == 'mag'
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.radius, 100 * u.km)
        assert isinstance(themis.H, Parameter)
        assert np.isclose(themis.H.value, 7.121)
        assert isinstance(themis.G12, Parameter)
        assert np.isclose(themis.G12.value, 0.68)
        assert themis.wfb == 'V'

    def test__G1_G2(self):
        themis = HG12(7.121 * u.mag, 0.68)
        assert np.isclose(themis._G1, 0.669592)
        assert np.isclose(themis._G2, 0.1407)
        themis = HG12(7.121 * u.mag, 0.1)
        assert np.isclose(themis._G1, 0.13691)
        assert np.isclose(themis._G2, 0.53088)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [7.121, 8.07252953, 8.67890827,
             9.2993879, 9.96817595, 10.72086969, 11.51208664, 12.12722017,
             18.70628001, 18.70647883])
        assert np.allclose(HG12.evaluate(pha_test, 7.121, 0.68),
                          phi_test)

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
        assert np.allclose(np.array(HG12.fit_deriv(
            pha_test, 7.121, 0.68)), phi_test)
        assert np.allclose(HG12.fit_deriv(0.2, 7.121, 0.68),
                [1.0, -0.07693564214597949])
        assert np.allclose(HG12.fit_deriv(0.2, 7.121, 0.1),
                [1.0, 0.6739785181393765])

    def test_props(self):
        themis = HG12(7.121 * u.mag, 0.68, radius=100 * u.km, wfb='V')
        assert np.isclose(themis.geomalb, 0.06215139)
        assert np.isclose(themis.bondalb, 0.02454096)
        assert np.isclose(themis.phaseint, 0.3948577512)
        assert np.isclose(themis.phasecoeff, -1.6777182566684201)
        assert np.isclose(themis.oe_amp, 0.23412300750840437)

    def test_from_obs(self):
        pha = [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789,
             31.57894737, 37.89473684, 44.21052632, 50.52631579, 56.84210526,
             63.15789474, 69.47368421, 75.78947368, 82.10526316, 88.42105263,
             94.73684211, 101.05263158, 107.36842105, 113.68421053,
             120.] * u.deg
        data = [6.95036472, 7.71609702, 8.04175457, 7.88226545, 8.28192813,
             8.50954834, 8.36880691, 8.73216696, 8.90742914, 9.05696656,
             9.20869753, 9.52578025, 9.8427691, 9.91588852, 10.3636637,
             10.26459992, 10.79316978, 10.79202241, 11.36950747,
             11.61018708] * u.mag
        from astropy.modeling.fitting import LevMarLSQFitter
        fitter = LevMarLSQFitter()
        m = HG12.from_obs({'alpha': pha, 'mag': data}, fitter)
        assert isinstance(m, HG12)
        assert isinstance(m.H, Parameter) & np.isclose(
            m.H.value, 7.13939) & (m.H.unit == u.mag)
        assert isinstance(m.G12, Parameter) & np.isclose(
            m.G12.value, 0.44872) & (m.G12.unit == u.dimensionless_unscaled)

    def test_to_mag(self):
        pha_test = np.linspace(0, np.pi, 10) * u.rad
        mag_test = [7.121, 8.07252953, 8.67890827, 9.2993879, 9.96817595,
                10.72086969, 11.51208664, 12.12722017, 18.70628001,
                18.70647883] * u.mag
        themis = HG12(7.121 * u.mag, 0.68)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(themis.to_mag(pha_test), mag_test)

    def test_to_ref(self):
        pha_test = np.linspace(0, np.pi, 10) * u.rad
        ref_test = [1.97834009e-02, 8.23548424e-03, 4.71126618e-03,
            2.66039298e-03, 1.43691333e-03, 7.18378086e-04, 3.46630119e-04,
            1.96703860e-04, 4.59397839e-07, 4.59313722e-07] / u.sr
        themis = HG12(7.121 * u.mag, 0.68, radius=100 * u.km, wfb='V')
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(themis.to_ref(pha_test), ref_test)

    def test_g_validator(self):
        with pytest.warns(NonmonotonicPhaseFunctionWarning):
            m = HG12(0, -0.71)
            m = HG12(0, 1.31)

    def test_G1_G2(self):
        themis = HG12(7.121 * u.mag, 0.68, radius=100 * u.km, wfb='V')
        assert np.isclose(themis._G1, 0.669592)
        assert np.isclose(themis._G2, 0.1407)


class TestHG12_Pen16:
    def test_init(self):
        themis = HG12_Pen16(7.121 * u.mag, 0.68, radius=100 * u.km, wfb='V')
        assert themis._unit == 'mag'
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.isclose(themis.radius, 100*u.km)
        assert isinstance(themis.H, Parameter)
        assert np.isclose(themis.H.value, 7.121)
        assert isinstance(themis.G12, Parameter)
        assert np.isclose(themis.G12.value, 0.68)
        assert themis.wfb == 'V'

    def test__G1_G2(self):
        themis = HG12_Pen16(7.121 * u.mag, 0.68)
        assert np.isclose(themis._G1, 0.5731968132)
        assert np.isclose(themis._G2, 0.17124272)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array(
            [7.121, 8.07252953, 8.67890827, 9.2993879, 9.96817595,
             10.72086969, 11.51208664, 12.12722017, 18.70628001, 18.70647883])
        assert np.allclose(HG12_Pen16.evaluate(pha_test, 7.121, 0.68),
                          phi_test)

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
        assert np.allclose(np.array(HG12_Pen16.fit_deriv(
            pha_test, 7.121, 0.68)), phi_test)
        assert np.allclose(HG12_Pen16.fit_deriv(0.2, 7.121, 0.68),
                [1., -0.08302351])

    def test_props(self):
        themis = HG12_Pen16(7.121 * u.mag, 0.68, radius=100 * u.km, wfb='V')
        assert np.isclose(themis.geomalb, 0.06215139)
        assert np.isclose(themis.bondalb, 0.02364406)
        assert np.isclose(themis.phaseint, 0.38042683486452)

    def test_from_obs(self):
        pha = [0., 6.31578947, 12.63157895, 18.94736842, 25.26315789,
             31.57894737, 37.89473684, 44.21052632, 50.52631579, 56.84210526,
             63.15789474, 69.47368421, 75.78947368, 82.10526316, 88.42105263,
             94.73684211, 101.05263158, 107.36842105, 113.68421053,
             120.] * u.deg
        data = [7.15663893, 7.4389134, 8.00006177, 7.9044872, 8.16865497,
             8.51010016, 8.63386712, 8.65893367, 8.84895152, 9.24495642,
             9.16195702, 9.54770054, 9.60599559, 10.06129054, 10.22544773,
             10.49122575, 10.78544483, 11.12145723, 11.18055954,
             11.40468613] * u.mag
        from astropy.modeling.fitting import LevMarLSQFitter
        fitter = LevMarLSQFitter()
        m = HG12_Pen16.from_obs({'alpha': pha, 'mag': data}, fitter)
        assert isinstance(m, HG12_Pen16)
        assert isinstance(m.H, Parameter) & np.isclose(
            m.H.value, 7.091456) & (m.H.unit == u.mag)
        assert isinstance(m.G12, Parameter) & np.isclose(
            m.G12.value, 0.631243) & (m.G12.unit == u.dimensionless_unscaled)

    def test_to_mag(self):
        pha_test = np.linspace(0, np.pi, 10) * u.rad
        mag_test = [7.121, 8.07252953, 8.67890827, 9.2993879, 9.96817595,
            10.72086969, 11.51208664, 12.12722017, 18.70628001,
            18.70647883] * u.mag
        themis = HG12_Pen16(7.121 * u.mag, 0.68)
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(themis.to_mag(pha_test), mag_test)

    def test_to_ref(self):
        pha_test = np.linspace(0, np.pi, 10) * u.rad
        ref_test = [1.97834009e-02, 8.23548424e-03, 4.71126618e-03,
            2.66039298e-03, 1.43691333e-03, 7.18378086e-04, 3.46630119e-04,
            1.96703860e-04, 4.59397839e-07, 4.59313722e-07] / u.sr
        themis = HG12_Pen16(7.121 * u.mag, 0.68, radius=100 * u.km, wfb='V')
        if LooseVersion(astropy.__version__) >= req_ver:
            assert u.allclose(themis.to_ref(pha_test), ref_test)

    def test_g_validator(self):
        with pytest.warns(NonmonotonicPhaseFunctionWarning):
            m = HG12(0, -0.71)
            m = HG12(0, 1.31)
