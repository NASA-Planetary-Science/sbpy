# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy import units as u
from astropy.modeling import Parameter
from ..core import *

def test_ref2mag():
    assert np.isclose(ref2mag(0.1/np.pi, 460), 3.3208379018205285)
    assert u.isclose(ref2mag(0.1/np.pi, 460*u.km), 3.3208379018205285*u.mag)


def test_mag2ref():
    assert np.isclose(mag2ref(3.32, 460), 0.03185556322261901)
    assert u.isclose(mag2ref(3.32, 460*u.km), 0.03185556322261901/u.sr)


def test_spline():
    nodes = np.array([10, 20, 30, 40, 50])
    values = np.array([1.2, 3.5, 8.9, 15.3, 30])
    dy = np.array([0.5, 2])
    y = spline(nodes, values, dy)
    x_test = np.linspace(5,55,10)
    y_test = np.array([-1.3       ,  1.45926293,  2.73994342,  4.12536376,
        7.223823  , 10.4137113 , 13.8713955 , 19.75606139, 28.90522609,
        40.        ])
    assert np.isclose(y(x_test), y_test).all()


class TestHG:
    def test_init(self):
        ceres = HG(3.4, 0.12, radius=480*u.km, M_sun=-26.74)
        assert u.isclose(ceres.radius, 480*u.km)
        assert np.isclose(ceres.M_sun, -26.74)
        assert isinstance(ceres.H, Parameter)
        assert np.isclose(ceres.H.value, 3.4)
        assert isinstance(ceres.G, Parameter)
        assert np.isclose(ceres.G.value, 0.12)

    def test_evaluate(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array([ 3.4       ,  4.43957695,  5.03935011,
                5.61797517,  6.26811269, 7.09560828,  8.2895693 , 10.25445176,
                14.35255427,         np.inf])
        assert np.isclose(HG.evaluate(pha_test, 3.4, 0.12), phi_test).all()

    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi*0.999, 10)
        deriv_test = np.array([[ 1.        ,  1.        ,  1.        ,
            1.        ,  1.        ,  1.        ,  1.        ,  1.        ,
            1.        ,  1.        ],
            [ 0.        , -1.33897222, -2.00745805, -2.43171634, -2.6155218 ,
            -2.46808816, -1.71295215,  0.06121683,  1.22716693,  1.23379114]])
        assert np.isclose(np.array(HG.fit_deriv(pha_test, 3.4, 0.12)), deriv_test).all()

    def test__check_unit(self):
        ceres = HG(3.4, 0.12, radius=480*u.km, M_sun=-26.74)
        assert ceres._unit == 'mag'

    def test_props(self):
        ceres = HG(3.4, 0.12, radius=480)
        assert np.isclose(ceres.geoalb, 0.08538239826749969)
        assert np.isclose(ceres.bondalb, 0.031109125948834842)
        assert np.isclose(ceres.phaseint, 0.36435057552929323)
        ceres.radius = 480*u.km
        assert u.isclose(ceres.geoalb, 0.0853824*u.dimensionless_unscaled)
        assert u.isclose(ceres.bondalb, 0.03110913*u.dimensionless_unscaled)
        assert np.isclose(ceres.phaseint, 0.36435057552929323)

    def test_fit(self):
        pass

    def test_distance_module(self):
        pass

    def test_mag(self):
        ceres = HG(3.4, 0.12, radius=480)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array([ 3.4       ,  4.43957695,  5.03935011,
            5.61797517,  6.26811269, 7.09560828,  8.2895693 , 10.25445176,
            14.35255427,         np.inf])
        assert np.isclose(ceres.mag(pha_test), mag_test).all()

    def test_ref(self):
        ceres = HG(3.4, 0.12, radius=480)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array([2.71780615e-02, 1.04324833e-02, 6.00452188e-03,
            3.52393922e-03, 1.93630337e-03, 9.03597797e-04, 3.00878224e-04,
            4.92535767e-05, 1.13030792e-06, 0.00000000e+00])
        ref_norm_test = np.array([1.00000000e+00, 3.83856785e-01,
            2.20932677e-01, 1.29661169e-01, 7.12450876e-02, 3.32473233e-02,
            1.10706286e-02, 1.81225496e-03, 4.15889823e-05, 0.00000000e+00])
        assert np.isclose(ceres.ref(pha_test), ref_test).all()
        assert np.isclose(ceres.ref(pha_test, normalized=0), ref_norm_test).all()


class TestHG1G2:
    def test_init(self):
        themis = HG1G2(7.063,0.62,0.14,radius=100.*u.km, M_sun=-26.74)
        assert themis._unit == 'mag'
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
        phi_test = np.array([ 7.063     ,  8.07436233,  8.68048572,
            9.29834638,  9.96574599,  10.72080704, 11.52317465, 12.15094612,
            18.65369516, 18.65389398])
        assert np.isclose(HG1G2.evaluate(pha_test, 7.063, 0.62, 0.14), phi_test).all()

    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi*0.999, 10)
        deriv_test = np.array([[ 1.00000000e+00,  1.00000000e+00,
            1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,
            1.00000000e+00,  1.00000000e+00,  1.00000000e+00, 1.00000000e+00],
            [ 2.43068538e-09, -1.21097348e+00, -1.20062175e+00,
            -1.14122248e+00, -1.10961252e+00, -1.16263148e+00, -1.41529341e+00,
            -1.70796793e+00,  0.00000000e+00, 0.00000000e+00],
            [-6.92019713e-10, -2.07410290e+00, -2.43821941e+00,
            -2.70127334e+00, -2.84126028e+00, -2.60646207e+00, -1.48753067e+00,
            -1.91400644e-01, -7.75525861e+00, -7.75525861e+00]])
        assert np.isclose(np.array(HG1G2.fit_deriv(pha_test, 7.063, 0.62, 0.14)), deriv_test).all()

    def test_props(self):
        themis = HG1G2(7.063,0.62,0.14,radius=100)
        assert np.isclose(themis.geoalb, 0.06739859193616704)
        assert np.isclose(themis.bondalb, 0.02521731797010077)
        assert np.isclose(themis.phaseint, 0.374152)
        themis.radius = 100*u.km
        assert u.isclose(themis.geoalb, 0.06739859*u.dimensionless_unscaled)
        assert u.isclose(themis.bondalb, 0.02521732*u.dimensionless_unscaled)
        assert np.isclose(themis.phaseint, 0.374152)

    def test_fit(self):
        pass

    def test_distance_module(self):
        pass

    def test_mag(self):
        themis = HG1G2(7.063,0.62,0.14,radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array([ 7.063     ,  8.07436233,  8.68048572,
            9.29834638,  9.96574599,  10.72080704, 11.52317465, 12.15094612,
            18.65369516, 18.65389398])
        assert np.isclose(themis.mag(pha_test), mag_test).all()

    def test_ref(self):
        themis = HG1G2(7.063,0.62,0.14,radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array([2.14536381e-02, 8.45193252e-03, 4.83622683e-03,
            2.73755213e-03,  1.48048003e-03, 7.38546998e-04, 3.52720817e-04,
            1.97843827e-04,  4.95704528e-07, 4.95613763e-07])
        assert np.isclose(themis.ref(pha_test), ref_test).all()


class TestHG12:
    def test_init(self):
        themis = HG12(7.121, 0.68, radius=100*u.km, M_sun=-26.74)
        assert themis._unit == 'mag'
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
        phi_test = np.array([ 7.121     ,  8.07252953,  8.67890827,
            9.2993879 ,  9.96817595,  10.72086969, 11.51208664, 12.12722017,
            18.70628001, 18.70647883])
        assert np.isclose(HG12.evaluate(pha_test, 7.121, 0.68), phi_test).all()


    def test_fit_deriv(self):
        pha_test = np.linspace(0, np.pi, 10)
        phi_test = np.array([[ 1.00000000e+00,  1.00000000e+00,
            1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,
            1.00000000e+00,  1.00000000e+00,  1.00000000e+00, 1.00000000e+00],
            [ 8.27405304e-10, -3.90147030e-01, -3.86830716e-01,
            -3.68566572e-01, -3.58865952e-01, -3.75398096e-01, -4.52786343e-01,
            -5.39541077e-01,  0.00000000e+00,  0.00000000e+00]])
        assert np.isclose(np.array(HG12.fit_deriv(pha_test, 7.121, 0.68)), phi_test).all()

    def test_props(self):
        themis = HG12(7.121,0.68,radius=100)
        assert np.isclose(themis.geoalb, 0.06389263856216909)
        assert np.isclose(themis.bondalb, 0.02522850358089249)
        assert np.isclose(themis.phaseint, 0.3948577512)
        themis.radius = 100*u.km
        assert u.isclose(themis.geoalb, 0.06389264*u.dimensionless_unscaled)
        assert u.isclose(themis.bondalb, 0.0252285*u.dimensionless_unscaled)
        assert np.isclose(themis.phaseint, 0.3948577512)

    def test_fit(self):
        pass

    def test_distance_module(self):
        pass

    def test_mag(self):
        themis = HG12(7.121,0.68,radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        mag_test = np.array([ 7.121     ,  8.07252953,  8.67890827,
            9.2993879 ,  9.96817595,  10.72086969, 11.51208664, 12.12722017,
            18.70628001, 18.70647883])
        assert np.isclose(themis.mag(pha_test), mag_test).all()

    def test_ref(self):
        themis = HG12(7.121,0.68,radius=100)
        pha_test = np.linspace(0, np.pi, 10)
        ref_test = np.array([2.03376585e-02, 8.46621202e-03, 4.84325842e-03,
            2.73492734e-03,  1.47717032e-03, 7.38504380e-04, 3.56341412e-04,
            2.02214774e-04,  4.72268466e-07, 4.72181992e-07])
        assert np.isclose(themis.ref(pha_test), ref_test).all()

