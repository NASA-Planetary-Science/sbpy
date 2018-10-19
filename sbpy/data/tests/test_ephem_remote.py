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
