# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from astropy import units as u
from astropy.modeling import Parameter
from ...calib import solar_fluxd
from ...data import Ephem
from ..core import *


def setup_module(module):
    module.solar_fluxd_default = solar_fluxd.get()
    solar_fluxd.set({'V': -26.77 * u.mag})


def teardown_module(module):
    solar_fluxd.set(module.solar_fluxd_default)


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
                    1.70003727e-03, 8.13688994e-04, 3.89456039e-04,
                    1.86405380e-04, 8.92192241e-05, 4.27030055e-05,
                    2.04389434e-05] / u.sr
        ref_norm_test = np.array(
            [1., 0.47863009, 0.22908677, 0.10964782, 0.05248075,
             0.02511886, 0.01202264, 0.0057544, 0.00275423,
             0.00131826]) * u.dimensionless_unscaled
        assert u.allclose(eph['ref'], ref_test)
        assert u.allclose(eph['alpha'], pha_test)
        assert set(eph.field_names) == {'alpha', 'ref'}
        eph_norm = linphase.to_ref(pha_test, normalized=0 * u.deg,
                                   append_results=True)
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
        pytest.importorskip("scipy")
        linphase = LinearPhaseFunc(5 * u.mag, 2.29 * u.mag/u.rad,
                                   radius=300 * u.km, wfb='V')
        assert np.isclose(linphase.geomalb, 0.0487089)
        assert np.isclose(linphase.bondalb, 0.01790315)
        assert np.isclose(linphase.phaseint, 0.36755394203990327)

    def test__distance_module(self):
        r = [0.5, 1, 1.2, 2] * u.au
        delta = [0.3, 1, 1, 2] * u.au
        m = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag / u.deg)
        module_test = [0.0225, 1., 1.44, 16.]
        module = m._distance_module(Ephem.from_dict({'r': r, 'delta': delta}))
        assert np.allclose(module, module_test)

    def test_fit(self):
        pytest.importorskip("scipy")
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
