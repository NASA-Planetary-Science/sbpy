# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

import astropy.units as u
from astropy.time import Time

from ... import exceptions as sbe
from .. import orbit as sbo
from ..orbit import Orbit

# CERES and CERES2 retrieved from Horizons on 24 Sep 2021
CERES = {
    'targetname': '1 Ceres (A801 AA)',
    'H': u.Quantity(3.53, 'mag'),
    'G': 0.12,
    'e': 0.07842518340409492,
    'q': u.Quantity(2.548847914875325, 'au'),
    'incl': u.Quantity(10.58812123343471, 'deg'),
    'Omega': u.Quantity(80.26788752671405, 'deg'),
    'w': u.Quantity(73.70547052298863, 'deg'),
    'n': u.Quantity(0.21428120879034, 'deg / d'),
    'M': u.Quantity(266.0066794881785, 'deg'),
    'nu': u.Quantity(257.137950098325, 'deg'),
    'a': u.Quantity(2.765752567209041, 'AU'),
    'Q': u.Quantity(2.982657219542757, 'AU'),
    'P': u.Quantity(1680.035323826441, 'd'),
    'epoch': Time(2459482.454237461, scale='tdb', format='jd'),
    'Tp': Time(2459921.098955971, scale='tdb', format='jd'),
    'orbtype': 'KEP',
}

CERES2 = {  # CERES epoch + 100 days
    'targetname': "1 Ceres (A801 AA)",
    'H': u.Quantity('3.53 mag'),
    'G': 0.12,
    'e': 0.0784882058716986,
    'q': u.Quantity('2.54890665 AU'),
    'incl': u.Quantity('10.58775349 deg'),
    'Omega': u.Quantity('80.26856263 deg'),
    'w': u.Quantity('73.64403121 deg'),
    'n': u.Quantity('0.21425182 deg / d'),
    'M': u.Quantity('287.5012003 deg'),
    'nu': u.Quantity('278.6978862 deg'),
    'a': u.Quantity('2.76600546 AU'),
    'Q': u.Quantity('2.98310427 AU'),
    'P': u.Quantity('1680.26575599 d'),
    'epoch': Time(2459582.45425436, scale='tdb', format='jd'),
    'Tp': Time(2459920.8355057, scale='tdb', format='jd')
}


class TestOOTransform:
    def test_oo_transform_kep(self):
        """ test oo_transform method"""

        pytest.importorskip("pyoorb")

        orbit = Orbit.from_dict(CERES)

        # transform to cart and back
        cart_orbit = orbit.oo_transform('CART')
        kep_orbit = cart_orbit.oo_transform('KEP')
        elements = ['a', 'e', 'i', 'Omega', 'w', 'M']
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd,
                         kep_orbit['epoch'][0].utc.jd)

        # transform to com
        com_orbit = orbit.oo_transform('COM')
        assert u.isclose(CERES['q'], com_orbit['q'])
        assert u.isclose(CERES['Tp'].utc.jd, com_orbit['Tp'].utc.jd)

        # and back
        kep_orbit = com_orbit.oo_transform('KEP')
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd,
                         kep_orbit['epoch'][0].utc.jd)
        assert kep_orbit['epoch'].scale == 'tdb'

    def test_oo_transform_com(self):
        """ test oo_transform method"""

        pytest.importorskip("pyoorb")

        orbit = Orbit.from_dict(CERES)

        # transform to cart and back
        cart_orbit = orbit.oo_transform('CART')
        com_orbit = cart_orbit.oo_transform('COM')
        elements = ['q', 'e', 'i', 'Omega', 'w']
        assert all([u.isclose(orbit[k][0], com_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd,
                         com_orbit['epoch'][0].utc.jd)
        assert u.isclose(orbit['Tp'][0].utc.jd,
                         com_orbit['Tp'][0].utc.jd)

        # transform to kep
        kep_orbit = orbit.oo_transform('KEP')
        assert u.isclose(CERES['a'], kep_orbit['a'])
        assert u.isclose(CERES['M'], kep_orbit['M'])

        # and back
        com_orbit = kep_orbit.oo_transform('COM')
        assert all([u.isclose(orbit[k][0], com_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd,
                         com_orbit['epoch'][0].utc.jd)
        assert u.isclose(orbit['Tp'][0].utc.jd,
                         com_orbit['Tp'][0].utc.jd)

        assert com_orbit['epoch'].scale == 'tdb'
        assert com_orbit['Tp'].scale == 'tdb'

    def test_timescales(self):
        """ test with input in UTC scale """
        pytest.importorskip("pyoorb")

        orbit = Orbit.from_dict(CERES)
        orbit['epoch'] = orbit['epoch'].utc

        # transform to cart and back
        cart_orbit = orbit.oo_transform('CART')
        kep_orbit = cart_orbit.oo_transform('KEP')
        elements = ['a', 'e', 'i', 'Omega', 'w', 'M']
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd,
                         kep_orbit['epoch'][0].utc.jd)

        # transform to com and back
        com_orbit = orbit.oo_transform('COM')
        kep_orbit = com_orbit.oo_transform('KEP')
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd,
                         kep_orbit['epoch'][0].utc.jd)
        assert kep_orbit['epoch'].scale == 'utc'


class TestOOPropagate:

    def test_oo_propagate(self):
        """ test oo_propagate method"""

        pytest.importorskip("pyoorb")

        orbit = Orbit.from_dict(CERES)
        future_orbit = Orbit.from_dict(CERES2)

        oo_orbit = orbit.oo_propagate(CERES2['epoch'])

        elements = ['a', 'e', 'i', 'Omega', 'w', 'M']
        assert all([u.isclose(oo_orbit[k][0], future_orbit[k][0])
                   for k in elements])
        assert u.isclose(oo_orbit['epoch'][0].utc.jd,
                         future_orbit['epoch'][0].utc.jd)
        assert oo_orbit['epoch'].scale == 'tdb'


# Halley retrieved from Horizons on 26 Apr 2022
HALLEY = {
    'targetname': '1P/Halley',
    'M1': u.Quantity(5.5, 'mag'),
    'M2': u.Quantity(13.6, 'mag'),
    'k1': 8.0,
    'k2': 5.0,
    'phasecoeff': u.Quantity(0.03, 'mag / deg'),
    'e': 0.9663712597171566,
    'q': u.Quantity(0.60100362, 'au'),
    'incl': u.Quantity(162.21355489, 'deg'),
    'Omega': u.Quantity(59.06301115, 'deg'),
    'w': u.Quantity(111.98929651, 'deg'),
    'n': u.Quantity(0.01304531, 'deg / d'),
    'M': u.Quantity(172.4302269, 'deg'),
    'nu': u.Quantity(179.49619633, 'deg'),
    'a': u.Quantity(17.87172549, 'au'),
    'Q': u.Quantity(35.14244737, 'au'),
    'P': u.Quantity(27596.12824787, 'd'),
    'epoch': Time(2459696.54292208, scale='tdb', format='jd'),
    'Tp': Time(2446478.74665706, scale='tdb', format='jd'),
    'orbtype': 'KEP'
}

# 252P and 2016 BA14 retrieved from Horizons on 26 Apr 2022
COMETS = {
    'targetname': ['252P/LINEAR', 'PANSTARRS (P/2016 BA14)'],
    'M1': u.Quantity([16.4, 21.3], 'mag'),
    'e': [0.671813410097529, 0.6650394891797625],
    'k1': [16.5, 4.5],
    'q': u.Quantity([1.00097003, 1.01345464], 'au'),
    'incl': u.Quantity([10.41045416, 18.89502003], 'deg'),
    'Omega': u.Quantity([190.93492189, 180.51874455], 'deg'),
    'w': u.Quantity([343.32963355, 351.92333244], 'deg'),
    'n': u.Quantity([0.18503493, 0.18727854], 'deg / d'),
    'M': u.Quantity([53.61784042, 58.76185268], 'deg'),
    'nu': u.Quantity([133.73260496, 136.4433599], 'deg'),
    'a': u.Quantity([3.05000283, 3.02559439], 'au'),
    'Q': u.Quantity([5.09903563, 5.03773413], 'au'),
    'P': u.Quantity([1945.5786428, 1922.27043206], 'd'),
    'M2': u.Quantity([np.nan, 23.1], 'mag'),
    'k2': [np.nan, 5.0],
    'phasecoeff': u.Quantity([np.nan, 0.03], 'mag / deg'),
    'epoch': Time([2459696.55177132, 2459696.55177132], scale='tdb',
                  format='jd'),
    'Tp': Time([2459406.78031243, 2459382.78462704], scale='tdb', format='jd'),
    'orbtype': ['KEP', 'KEP']
 }

# Jupiter retrieved from Horizons on 26 Apr 2022
JUPITER = {
    'targetname': 'Jupiter (599)',
    'e': 0.04833017971265399,
    'q': u.Quantity(4.95081216, 'au'),
    'incl': u.Quantity(1.3036833, 'deg'),
    'Omega': u.Quantity(100.52981374, 'deg'),
    'w': u.Quantity(273.5645208, 'deg'),
    'n': u.Quantity(0.08310479, 'deg / d'),
    'M': u.Quantity(337.58379998, 'deg'),
    'nu': u.Quantity(335.34786485, 'deg'),
    'a': u.Quantity(5.20223722, 'au'),
    'Q': u.Quantity(5.45366228, 'au'),
    'P': u.Quantity(4331.8804878, 'd'),
    'epoch': Time(2459696.54405449, scale='tdb', format='jd'),
    'Tp': Time(2459966.27821974, scale='tdb', format='jd'),
    'orbtype': 'KEP'
}


class TestOrbitTisserandDCriterion:
    def testTisserand(self):
        """ test tisserand method"""

        halley = Orbit.from_dict(HALLEY)
        jupiter = Orbit.from_dict(JUPITER)
        assert u.isclose(halley.tisserand(jupiter), u.Quantity(-0.61659744))
        comets = Orbit.from_dict(COMETS)
        assert u.allclose(comets.tisserand(jupiter),
                          u.Quantity([2.82130706, 2.79709686]))

    @pytest.mark.parametrize('e, q, w, Omega, i, D_SH, D_D', (
        (0.7457, 0.95 * u.au, 152.000 * u.deg, 45.0 * u.deg,
            30.0 * u.deg, 0.0, 0.0),
        (0.9335, 0.95 * u.au, 153.694 * u.deg, 42.0 * u.deg,
            32.0 * u.deg, 0.193, 0.113),
        (0.5954, 0.95 * u.au, 150.000 * u.deg, 46.0 * u.deg,
            32.0 * u.deg, 0.155, 0.113),
        (0.5568, 0.95 * u.au, 149.339 * u.deg, 46.0 * u.deg,
            32.0 * u.deg, 0.193, 0.146)
    ))
    def testDCriterion(self, e, q, w, Omega, i, D_SH_test, D_D_test):
        """test D_criterion method against Table 1 in Jopek 1993, Icarus
        106, 603
        """

        orb0 = Orbit.from_dict({
            'e': 0.7457,
            'q': 0.95 * u.au,
            'w': 152 * u.deg,
            'Omega': 45.0 * u.deg,
            'incl': 30.0 * u.deg
            })
        orb = Orbit.from_dict({
            'e': e,
            'q': q,
            'w': w,
            'Omega': Omega,
            'incl': i
            })
        D_SH = orb0.D_criterion(orb, version='sh')
        D_D = orb0.D_criterion(orb, version='d')
        assert u.isclose(D_SH, D_SH_test, atol=1e-3)
        assert u.isclose(D_D, D_D_test, atol=1e-3)

    def testDCriterionException(self):
        orb1 = Orbit.from_dict(HALLEY)
        orb2 = Orbit.from_dict(JUPITER)
        with pytest.raises(ValueError):
            orb1.D_criterion(orb2, version='test')

    def testDCriterion(self):
        """test D_criterion method with higher precision"""
        comets = Orbit.from_dict(COMETS)
        assert u.isclose(comets[0].D_criterion(comets[1]),
                         u.Quantity(0.15598291))
        assert u.isclose(comets[0].D_criterion(comets[1], version='d'),
                         u.Quantity(0.0502142))
        assert u.isclose(comets[0].D_criterion(comets[1], version='h'),
                         u.Quantity(0.15560596))
        halley = Orbit.from_dict(HALLEY)
        assert u.allclose(halley.D_criterion(comets),
                          u.Quantity([2.24077662, 2.50705989]))
        assert u.allclose(halley.D_criterion(comets, version='d'),
                          u.Quantity([1.14268444, 1.1268579]))
        assert u.allclose(halley.D_criterion(comets, version='h'),
                          u.Quantity([2.21888313, 2.48606116]))
