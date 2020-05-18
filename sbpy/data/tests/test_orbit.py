# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from numpy.testing import assert_allclose
import astropy.units as u
from astropy.time import Time
import warnings

from ... import exceptions as sbe
from .. import orbit as sbo
from ..orbit import Orbit, QueryError
from ..names import TargetNameParseError
from ... import bib

try:
    import pyoorb
    HAS_PYOORB = True
except ImportError:
    HAS_PYOORB = False

# retreived from Horizons on 23 Apr 2020
CERES = {
    'targetname': '1 Ceres',
    'H': u.Quantity(3.4, 'mag'),
    'G': 0.12,
    'e': 0.07741102040801928,
    'q': u.Quantity(2.55375156, 'au'),
    'incl': u.Quantity(10.58910839, 'deg'),
    'Omega': u.Quantity(80.29081558, 'deg'),
    'w': u.Quantity(73.7435117, 'deg'),
    'n': u.Quantity(0.21401711, 'deg / d'),
    'M': u.Quantity(154.70418799, 'deg'),
    'nu': u.Quantity(158.18663933, 'deg'),
    'a': u.Quantity(2.76802739, 'AU'),
    'Q': u.Quantity(2.98230321, 'AU'),
    'P': u.Quantity(1682.10848349, 'd'),
    'epoch': Time(2458963.26397076, scale='tdb', format='jd'),
    'Tp': Time(2458240.40500675, scale='tdb', format='jd'),
    'orbtype': 'KEP',
}


@pytest.mark.skipif('not HAS_PYOORB')
class TestOOTransform:
    def test_missing_pyoorb(self, monkeypatch):
        monkeypatch.setattr(sbo, 'pyoorb', None)
        with pytest.raises(sbe.RequiredPackageUnavailable):
            Orbit.from_dict(CERES).oo_transform('CART')

    def test_oo_transform(self):
        """ test oo_transform method"""

        orbit = Orbit.from_dict(CERES)

        # transform to cart and back
        cart_orbit = orbit.oo_transform('CART')
        kep_orbit = cart_orbit.oo_transform('KEP')
        elements = ['a', 'e', 'i', 'Omega', 'w', 'M']
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

        # transform to com and back
        com_orbit = orbit.oo_transform('COM')
        kep_orbit = com_orbit.oo_transform('KEP')
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

    def test_timescales(self):
        """ test with input in UTC scale """
        orbit = Orbit.from_dict(CERES)
        orbit['epoch'] = orbit['epoch'].utc

        # transform to cart and back
        cart_orbit = orbit.oo_transform('CART')
        kep_orbit = cart_orbit.oo_transform('KEP')
        elements = ['a', 'e', 'i', 'Omega', 'w', 'M']
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

        # transform to com and back
        com_orbit = orbit.oo_transform('COM')
        kep_orbit = com_orbit.oo_transform('KEP')
        assert all([u.isclose(orbit[k][0], kep_orbit[k][0]) for k in elements])
        assert u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

        assert kep_orbit['epoch'].scale == 'tdb'


@pytest.mark.skipif('not HAS_PYOORB')
class TestOOPropagate:
    def test_missing_pyoorb(self, monkeypatch):
        monkeypatch.setattr(sbo, 'pyoorb', None)
        with pytest.raises(sbe.RequiredPackageUnavailable):
            Orbit.from_dict(CERES).oo_transform('CART')
