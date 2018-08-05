# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from numpy.testing import assert_allclose
import astropy.units as u
from astropy.time import Time

from sbpy.data import Orbit
from sbpy import bib


@pytest.mark.remote_data
def test_from_horizons():
    """ test from_horizons method"""

    # current epoch
    now = Time.now()
    data = Orbit.from_horizons('Ceres')
    assert_allclose(data.datetime_jd, now.jd*u.d)

    # # date range - astropy.time.Time objects
    # epochs = {'start': Time('2018-01-02', format='iso'),
    #           'stop': Time('2018-01-05', format='iso'),
    #           'step': '6h'}
    # data = Orbit.from_horizons('Ceres', epochs=epochs)
    # assert len(data.table) == 13

    # # date range - strings
    # epochs = {'start': '2018-01-02',
    #           'stop': '2018-01-05',
    #           'step': '6h'}
    # data = Orbit.from_horizons('Ceres', epochs=epochs)
    # assert len(data.table) == 13

    # discrete epochs - astropy.time.Time objects
    epochs = [Time('2018-01-02', format='iso'),
              Time('2018-01-05', format='iso')]
    data = Orbit.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 2

    # discrete epochs - Julian Dates
    epochs = [Time('2018-01-02', format='iso').jd,
              Time('2018-01-05', format='iso').jd]
    data = Orbit.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 2

    # query two objects
    data = Orbit.from_horizons(['Ceres', 'Pallas'])
    assert len(data.table) == 2

    # test bib service
    bib.track()
    data = Orbit.from_horizons(['Ceres', 'Pallas'])
    assert bib.to_text() == ('sbpy.data.Orbit:\n  '
                             'data service: 1996DPS....28.2504G\n')
