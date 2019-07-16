# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from copy import deepcopy

from numpy.testing import assert_allclose
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation

from sbpy.data import Ephem, Orbit
from sbpy import bib
from ..core import conf


@pytest.mark.remote_data
class TestEphemFromHorizons:

    def test_current_epoch(self):
        now = Time.now()
        data = Ephem.from_horizons('Ceres')
        assert_allclose(data['datetime_jd'], now.jd*u.d)

    def test_daterange_Time(self):
        epochs = {'start': Time('2018-01-02', format='iso'),
                  'stop': Time('2018-01-05', format='iso'),
                  'step': '6h'}
        data = Ephem.from_horizons('Ceres', epochs=epochs)
        assert len(data.table) == 13

    def test_daterange_strings(self):
        # (this is not really supported by sbpy)
        epochs = {'start': '2018-01-02',
                  'stop': '2018-01-05',
                  'step': '6h'}
        data = Ephem.from_horizons('Ceres', epochs=epochs)
        assert len(data.table) == 13

    def test_daterange_number(self):
        epochs = {'start': Time('2018-01-02', format='iso'),
                  'stop': Time('2018-01-05', format='iso'),
                  'number': 10}
        data = Ephem.from_horizons('Ceres', epochs=epochs)
        assert len(data.table) == 10

    def test_Time_list(self):
        epochs = Time(['2018-01-01', '2018-01-02'])
        data = Ephem.from_horizons('Ceres', epochs=epochs)
        assert len(data.table) == 2

    def test_multiple_objects(self):
        data = Ephem.from_horizons(['Ceres', 'Pallas'])
        assert len(data.table) == 2

    def test_Earthlocation(self):
        lowell = EarthLocation.of_site('Lowell Observatory')
        data = Ephem.from_horizons(1, epochs=Time('2018-01-01',
                                                  format='iso'),
                                   location=lowell)
        assert len(data.table) == 1

    def test_bib(self):
        bib.track()
        data = Ephem.from_horizons(['Ceres', 'Pallas'])
        assert 'sbpy.data.Ephem' in bib.to_text()


@pytest.mark.remote_data
class TestEphemFromMPC:
    def test_single_epoch_now(self):
        eph = Ephem.from_mpc('Ceres')
        assert len(eph.table) == 1

    def test_single_epoch(self):
        eph = Ephem.from_mpc('Ceres', epochs='2018-10-01')
        assert len(eph.table) == 1

    def test_multiple_epochs(self):
        eph = Ephem.from_mpc('Ceres', epochs=['2018-10-01', '2019-10-01'])
        assert len(eph.table) == 2

    def test_start_stop_step(self):
        epochs = dict(start='2018-10-01', stop='2018-10-31', step='1d')
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert len(eph.table) == 31

    def test_minute_steps_pr88(self):
        """https://github.com/NASA-Planetary-Science/sbpy/pull/88"""
        epochs = dict(start='2018-10-01', step='1min', number=10)
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert len(eph.table) == 10

    def test_start_stop_no_step(self):
        with pytest.raises(ValueError):
            eph = Ephem.from_mpc('Ceres', epochs={'start': '2018-10-01',
                                                  'stop': '2018-10-31'})

    def test_start_step_number(self):
        epochs = dict(start='2018-10-01', step='1d', number=31)
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert len(eph.table) == 31
        assert eph['Date'][-1] == '2018-10-31 00:00:00.000'

    def test_start_stop_jd(self):
        epochs = {'start': 2458396.5, 'stop': 2458397.5, 'step': '1d'}
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert eph['Date'][0] == '2018-10-05 00:00:00.000'
        assert eph['Date'][1] == '2018-10-06 00:00:00.000'

    def test_epochs_jd(self):
        epochs = ['2018-10-05', 2458397.5]
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert eph['Date'][1] == '2018-10-06 00:00:00.000'

    def test_step_unit(self):
        with pytest.raises(ValueError):
            eph = Ephem.from_mpc('Ceres', epochs={'start': '2018-10-01',
                                                  'step': '1yr'})

    def test_ra_dec_format(self):
        epochs = dict(start='2018-10-01', step='1d', number=31)
        ra_format = {'sep': ':', 'unit': 'hourangle', 'precision': 1}
        dec_format = {'sep': ':', 'precision': 1}
        eph = Ephem.from_mpc('Ceres', epochs=epochs, ra_format=ra_format,
                             dec_format=dec_format)
        assert isinstance(eph['RA'][0], str)
        assert isinstance(eph['Dec'][0], str)

    def test_multiple_targets(self):
        eph = Ephem.from_mpc(['Ceres', 'Pallas', 'Vesta'])
        assert all(eph['Targetname'].data == ['Ceres', 'Pallas', 'Vesta'])


@pytest.mark.remote_data
def test_from_oo():
    """test from_oo method"""

    try:
        import pyoorb
    except ImportError:
        return None

    orbit = Orbit.from_horizons('Ceres')
    horizons_ephem = Ephem.from_horizons('Ceres', location='500')
    oo_ephem = Ephem.from_oo(orbit)

    u.isclose(horizons_ephem['ra'][0], oo_ephem['ra'][0])
    u.isclose(horizons_ephem['dec'][0], oo_ephem['dec'][0])
    u.isclose(horizons_ephem['RA*cos(Dec)_rate'][0],
              oo_ephem['RA*cos(Dec)_rate'][0])
    u.isclose(horizons_ephem['dec_rate'][0], oo_ephem['dec_rate'][0])
    u.isclose(horizons_ephem['alpha'][0], oo_ephem['alpha'][0])
    u.isclose(horizons_ephem['r'][0], oo_ephem['r'][0])
    u.isclose(horizons_ephem['delta'][0], oo_ephem['delta'][0])
    u.isclose(horizons_ephem['V'][0], oo_ephem['V'][0])
    u.isclose(horizons_ephem['hlon'][0], oo_ephem['hlon'][0])
    u.isclose(horizons_ephem['hlat'][0], oo_ephem['hlat'][0])
    u.isclose(horizons_ephem['EL'][0], oo_ephem['EL'][0])

    # test manual orbit definition lacking units
    manorbit = Orbit.from_dict({
        'targetname': orbit['targetname'][0],
        'a': orbit['a'].value[0],
        'e': orbit['e'][0],
        'i': orbit['i'].value[0],
        'w': orbit['w'].value[0],
        'Omega': orbit['Omega'].value[0],
        'datetime_jd': orbit['datetime_jd'][0].value,
        'M': orbit['M'].value[0],
        'H': orbit['H'].value[0],
        'G': orbit['G'][0],
        'timescale': orbit['timescale'][0]})

    oo_ephem = Ephem.from_oo(manorbit)

    u.isclose(horizons_ephem['ra'][0], oo_ephem['ra'][0])
    u.isclose(horizons_ephem['dec'][0], oo_ephem['dec'][0])
    u.isclose(horizons_ephem['RA*cos(Dec)_rate'][0],
              oo_ephem['RA*cos(Dec)_rate'][0])
    u.isclose(horizons_ephem['dec_rate'][0], oo_ephem['dec_rate'][0])
    u.isclose(horizons_ephem['alpha'][0], oo_ephem['alpha'][0])
    u.isclose(horizons_ephem['r'][0], oo_ephem['r'][0])
    u.isclose(horizons_ephem['delta'][0], oo_ephem['delta'][0])
    u.isclose(horizons_ephem['V'][0], oo_ephem['V'][0])
    u.isclose(horizons_ephem['hlon'][0], oo_ephem['hlon'][0])
    u.isclose(horizons_ephem['hlat'][0], oo_ephem['hlat'][0])
    u.isclose(horizons_ephem['EL'][0], oo_ephem['EL'][0])
