# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from numpy.testing import assert_allclose
import astropy.units as u
from astropy.time import Time

from sbpy.data import Ephem, Orbit
from sbpy import bib
from ..core import conf


@pytest.mark.remote_data
def test_from_horizons():
    """ test from_horizons method"""

    # current epoch
    now = Time.now()
    data = Ephem.from_horizons('Ceres')
    assert_allclose(data.datetime_jd, now.jd*u.d)

    # date range - astropy.time.Time objects
    epochs = {'start': Time('2018-01-02', format='iso'),
              'stop': Time('2018-01-05', format='iso'),
              'step': '6h'}
    data = Ephem.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 13

    # date range - strings
    epochs = {'start': '2018-01-02',
              'stop': '2018-01-05',
              'step': '6h'}
    data = Ephem.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 13

    # discrete epochs - astropy.time.Time objects
    epochs = [Time('2018-01-02', format='iso'),
              Time('2018-01-05', format='iso')]
    data = Ephem.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 2

    # discrete epochs - Julian Dates
    epochs = [Time('2018-01-02', format='iso').jd,
              Time('2018-01-05', format='iso').jd]
    data = Ephem.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 2

    # query two objects
    data = Ephem.from_horizons(['Ceres', 'Pallas'])
    assert len(data.table) == 2

    # test bib service
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
    u.isclose(horizons_ephem['ra_rate'][0], oo_ephem['ra_rate'][0])
    u.isclose(horizons_ephem['dec_rate'][0], oo_ephem['dec_rate'][0])
    u.isclose(horizons_ephem['alpha'][0], oo_ephem['alpha'][0])
    u.isclose(horizons_ephem['r'][0], oo_ephem['r'][0])
    u.isclose(horizons_ephem['delta'][0], oo_ephem['delta'][0])
    u.isclose(horizons_ephem['V'][0], oo_ephem['V'][0])
    u.isclose(horizons_ephem['hlon'][0], oo_ephem['hlon'][0])
    u.isclose(horizons_ephem['hlat'][0], oo_ephem['hlat'][0])
    u.isclose(horizons_ephem['EL'][0], oo_ephem['EL'][0])
