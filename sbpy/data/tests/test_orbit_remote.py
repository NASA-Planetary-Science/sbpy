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
    assert_allclose(data['datetime_jd'], now.jd*u.d)

    # date range - astropy.time.Time objects
    epochs = {'start': Time('2018-01-02', format='iso'),
              'stop': Time('2018-01-05', format='iso'),
              'step': '6h'}
    data = Orbit.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 13

    # date range - strings
    epochs = {'start': '2018-01-02',
              'stop': '2018-01-05',
              'step': '6h'}
    data = Orbit.from_horizons('Ceres', epochs=epochs)
    assert len(data.table) == 13

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
    assert 'sbpy.data.Orbit' in bib.to_text()


@pytest.mark.remote_data
def test_oo_transform():
    """ test oo_transform method"""

    try:
        import pyoorb
    except ImportError:
        return None

    orbit = Orbit.from_horizons('Ceres')

    cart_orbit = orbit.oo_transform('CART')

    kep_orbit = cart_orbit.oo_transform('KEP')
    u.isclose(orbit['a'][0], kep_orbit['a'][0])
    u.isclose(orbit['e'][0], kep_orbit['e'][0])
    u.isclose(orbit['i'][0], kep_orbit['i'][0])
    u.isclose(orbit['Omega'][0], kep_orbit['Omega'][0])
    u.isclose(orbit['w'][0], kep_orbit['w'][0])
    u.isclose(orbit['M'][0], kep_orbit['M'][0])
    u.isclose(orbit['epoch'][0], kep_orbit['epoch'][0])

    com_orbit = orbit.oo_transform('COM')
    kep_orbit = com_orbit.oo_transform('KEP')
    u.isclose(orbit['a'][0], kep_orbit['a'][0])
    u.isclose(orbit['e'][0], kep_orbit['e'][0])
    u.isclose(orbit['i'][0], kep_orbit['i'][0])
    u.isclose(orbit['Omega'][0], kep_orbit['Omega'][0])
    u.isclose(orbit['w'][0], kep_orbit['w'][0])
    u.isclose(orbit['M'][0], kep_orbit['M'][0])
    u.isclose(orbit['epoch'][0], kep_orbit['epoch'][0])


@pytest.mark.remote_data
def test_oo_propagate():
    """ test oo_propagate method"""

    try:
        import pyoorb
    except ImportError:
        return None

    orbit = Orbit.from_horizons('Ceres')
    epoch = Time.now().jd + 100

    future_orbit = Orbit.from_horizons('Ceres', epochs=epoch)

    oo_orbit = orbit.oo_propagate(epoch)

    u.isclose(oo_orbit['a'][0], future_orbit['a'][0])
    u.isclose(oo_orbit['e'][0], future_orbit['e'][0])
    u.isclose(oo_orbit['i'][0], future_orbit['i'][0])
    u.isclose(oo_orbit['Omega'][0], future_orbit['Omega'][0])
    u.isclose(oo_orbit['w'][0], future_orbit['w'][0])
    u.isclose(oo_orbit['M'][0], future_orbit['M'][0])
    u.isclose(oo_orbit['epoch'][0], future_orbit['epoch'][0])
