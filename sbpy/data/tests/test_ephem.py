# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import assert_quantity_allclose

from ... import exceptions as sbe
from ... import bib
from .. import ephem
from .. import Ephem, Orbit

try:
    import pyoorb
except ImportError:
    pyoorb = None

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
    'Tp': Time(2458240.40500675, scale='tdb', format='jd')
}


@pytest.mark.skipif('pyoorb is None')
class TestEphemFromOorb:
    def test_missing_pyoorb(self, monkeypatch):
        monkeypatch.setattr(ephem, 'pyoorb', None)
        with pytest.raises(sbe.RequiredPackageUnavailable):
            Ephem.from_oo(CERES)

    def test_units(self):
        orbit1 = Orbit.from_dict(CERES)
        eph1 = Ephem.from_oo(orbit1)

        orbit2 = Orbit.from_dict({
            'targetname': orbit1['targetname'][0],
            'a': orbit1['a'].value[0] * u.au,
            'e': orbit1['e'][0],
            'i': orbit1['i'].value[0] * u.deg,
            'w': orbit1['w'].value[0] * u.deg,
            'Omega': orbit1['Omega'].value[0] * u.deg,
            'epoch': Time(orbit1['epoch'][0], format='jd'),
            'M': orbit1['M'].value[0] * u.deg,
            'H': orbit1['H'].value[0] * u.mag,
            'G': orbit1['G'][0]
        })
        eph2 = Ephem.from_oo(orbit2)

        for k in ['ra', 'dec', 'RA*cos(Dec)_rate', 'dec_rate', 'alpha', 'r',
                  'delta', 'V', 'hlon', 'hlat', 'EL']:
            assert u.isclose(eph1[k], eph2[k])

    def test_basic(self):
        orbit = Orbit.from_dict(CERES)
        oo_ephem = Ephem.from_oo(orbit, scope='basic')
        assert 'dec_rate' not in oo_ephem.field_names

    def test_timescale(self):
        orbit = Orbit.from_dict(CERES)
        epoch = Time.now()
        oo_ephem = Ephem.from_oo(orbit, epochs=epoch, scope='basic')
        assert oo_ephem['epoch'].scale == epoch.scale

    def test_bib(self):
        with bib.Tracking():
            orbit = Orbit.from_dict(CERES)
            oo_ephem = Ephem.from_oo(orbit, scope='basic')
            assert 'sbpy.data.ephem.Ephem.from_oo' in bib.show()
        bib.reset()


class TestFillDeltaAndPhase:
    def test_co_orbital(self):
        """Test Earth co-orbital"""
        eph = Ephem.from_dict({'rh': np.ones(91) * u.au,
                               'solar_elongation': np.arange(91) * u.deg})
        observer_rh = 1 * u.au
        eph.fill_delta_and_phase(observer_rh=observer_rh)

        # verify they are valid triangles
        phase = np.arccos((eph['delta']**2 + eph['rh']**2 - observer_rh**2)
                          / 2 / eph['delta'] / eph['rh'])
        assert_quantity_allclose(phase, eph['phase'])

    def test_interior(self):
        """Test interior object"""
        observer_rh = 1 * u.au

        # given solar elongation, what rh are possible?
        for selong in [1, 10, 30, 45, 60, 80, 89] * u.deg:
            min_rh = observer_rh * np.sin(selong)
            rh = np.linspace(min_rh, observer_rh)[
                :-1]  # only use interior values

            eph = Ephem.from_dict({'rh': rh})
            eph['solar_elongation'] = selong
            eph.fill_delta_and_phase(closest=True)

            # verify they are valid triangles
            phase = np.arccos((eph['delta']**2 + eph['rh']**2 - observer_rh**2)
                              / 2 / eph['delta'] / eph['rh'])
            assert_quantity_allclose(phase, eph['phase'])

            # repeat for farthest solutions
            closest = eph['delta'].copy()
            eph.fill_delta_and_phase(closest=False, overwrite=True)
            phase = np.arccos((eph['delta']**2 + eph['rh']**2 - observer_rh**2)
                              / 2 / eph['delta'] / eph['rh'])
            assert_quantity_allclose(phase, eph['phase'])

            # first should be the same distance
            assert_quantity_allclose(closest[0], eph['delta'][0])
            # the rest should be different
            assert (closest[1:] < eph['delta'][1:]).all()

    def test_exterior(self):
        """Text exterior object"""
        observer_rh = 1 * u.au
        rh = (1 + np.logspace(-5, 2)) * u.au
        for selong in np.linspace(0.01, np.pi) * u.rad:
            eph = Ephem.from_dict({'rh': rh})
            eph['solar_elongation'] = selong
            eph.fill_delta_and_phase()

            # verify they are valid triangles
            if selong == np.pi * u.rad:
                phase = 0 * u.deg
            else:
                phase = np.arccos((eph['delta']**2 + eph['rh']**2 - observer_rh**2)
                                  / 2 / eph['delta'] / eph['rh'])
            assert_quantity_allclose(eph['phase'], phase)

    def test_exceptions(self):
        eph = Ephem.from_dict({'ra': 100 * u.deg, 'dec': -50 * u.deg})
        with pytest.raises(ValueError):
            eph.fill_delta_and_phase()

        eph['rh'] = 1 * u.au
        with pytest.raises(ValueError):
            eph.fill_delta_and_phase()

        eph['delta'] = 1 * u.au
        with pytest.raises(ValueError):
            eph.fill_delta_and_phase()

        del eph.table['delta']
        assert 'delta' not in eph
        eph['phase'] = 1 * u.deg
        with pytest.raises(ValueError):
            eph.fill_delta_and_phase()
