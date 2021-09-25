# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import astropy.units as u
from astropy.time import Time

from ... import exceptions as sbe
from .. import orbit as sbo
from ..orbit import Orbit

try:
    import pyoorb  # noqa
except ImportError:
    pyoorb = None

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


@pytest.mark.skipif('pyoorb is None')
class TestOOTransform:
    def test_missing_pyoorb(self, monkeypatch):
        monkeypatch.setattr(sbo, 'pyoorb', None)
        with pytest.raises(sbe.RequiredPackageUnavailable):
            Orbit.from_dict(CERES).oo_transform('CART')

    def test_oo_transform_kep(self):
        """ test oo_transform method"""

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


@pytest.mark.skipif('pyoorb is None')
class TestOOPropagate:
    def test_missing_pyoorb(self, monkeypatch):
        monkeypatch.setattr(sbo, 'pyoorb', None)
        with pytest.raises(sbe.RequiredPackageUnavailable):
            Orbit.from_dict(CERES).oo_transform('CART')

    def test_oo_propagate(self):
        """ test oo_propagate method"""

        orbit = Orbit.from_dict(CERES)
        future_orbit = Orbit.from_dict(CERES2)

        oo_orbit = orbit.oo_propagate(CERES2['epoch'])

        elements = ['a', 'e', 'i', 'Omega', 'w', 'M']
        assert all([u.isclose(oo_orbit[k][0], future_orbit[k][0])
                   for k in elements])
        assert u.isclose(oo_orbit['epoch'][0].utc.jd,
                         future_orbit['epoch'][0].utc.jd)
        assert oo_orbit['epoch'].scale == 'tdb'
