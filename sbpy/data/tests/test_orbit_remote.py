# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from numpy.testing import assert_allclose
import astropy.units as u
from astropy.time import Time
import warnings

from ..orbit import Orbit, QueryError
from ..names import TargetNameParseError
from ... import bib


@pytest.mark.remote_data
class TestOrbitFromHorizons:

    def test_now(self):
        # current epoch
        now = Time.now()
        data = Orbit.from_horizons('Ceres')
        assert_allclose(data['epoch'].utc.jd, now.utc.jd)

    def test_range_step(self):
        # date range - astropy.time.Time objects
        epochs = {'start': Time('2018-01-02', format='iso'),
                  'stop': Time('2018-01-05', format='iso'),
                  'step': 6*u.h}
        data = Orbit.from_horizons('Ceres', epochs=epochs)
        assert len(data.table) == 13

    def test_range_step(self):
        # date range - astropy.time.Time objects
        epochs = {'start': Time('2018-01-02', format='iso'),
                  'stop': Time('2018-01-05', format='iso'),
                  'number': 10}
        data = Orbit.from_horizons('Ceres', epochs=epochs)
        assert len(data.table) == 10

    def test_mult_epochs(self):
        # discrete epochs - astropy.time.Time objects
        epochs = Time(['2018-01-02', '2018-01-05'], format='iso')
        data = Orbit.from_horizons('Ceres', epochs=epochs)
        assert len(data.table) == 2

    def test_mult_objects(self):
        # query two objects
        data = Orbit.from_horizons(['Ceres', 'Pallas'])
        assert len(data.table) == 2

    def test_bib(self):
        # test bib service
        bib.track()
        Orbit.from_horizons(['Ceres', 'Pallas'])
        assert 'sbpy.data.orbit.Orbit.from_horizons' in bib.to_text()

    def test_queryfail(self):
        with pytest.raises(QueryError):
            Orbit.from_horizons('target does not exist')

    def test_timescale(self):
        # test same timescale
        time = Time('2020-01-01 00:00', format='iso', scale='tdb')
        a = Orbit.from_horizons(1, epochs=time)
        assert a['epoch'].scale == 'tdb'
        assert a['epoch'].tdb.jd == 2458849.5

        # test different timescale (and check for warning)
        time = Time('2020-01-01 00:00', format='iso', scale='utc')
        with warnings.catch_warnings(record=True) as w:
            a = Orbit.from_horizons(1, epochs=time)
            for i in range(len(w)):
                print(w[i])
            assert any(["astroquery.jplhorizons" in str(w[i].message)
                        for i in range(len(w))])
        assert a['epoch'].scale == 'tdb'
        assert_allclose(a['epoch'][0].tdb.jd, 2458849.49919926)

    def test_output_epochs(self):
        # all tdb
        epochs = Time(2437675.5, format='jd', scale='tdb')
        a = Orbit.from_horizons(1, epochs=epochs)
        assert a['epoch'][0].scale == 'tdb'

        # request tai, expect tdb and warning
        epochs = Time([2437655.5, 2437675.5], format='jd', scale='tai')
        with warnings.catch_warnings(record=True) as w:
            a = Orbit.from_horizons(1, epochs=epochs)
            assert any(["astroquery.jplhorizons" in str(w[i].message)
                        for i in range(len(w))])
        assert a['epoch'].scale == 'tdb'


@pytest.mark.remote_data
class TestOrbitFromMPC:

    def test_single(self):
        a = Orbit.from_mpc('Ceres')
        assert len(a) == 1

    def test_multiple(self):
        a = Orbit.from_mpc(['1P', '2P', '3P'])
        assert len(a) == 3

    def test_break(self):
        with pytest.raises(TargetNameParseError):
            a = Orbit.from_mpc('does not exist')


@pytest.mark.remote_data
class TestOOTransform:
    def test_oo_transform(self):
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
        u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

        com_orbit = orbit.oo_transform('COM')
        kep_orbit = com_orbit.oo_transform('KEP')
        u.isclose(orbit['a'][0], kep_orbit['a'][0])
        u.isclose(orbit['e'][0], kep_orbit['e'][0])
        u.isclose(orbit['i'][0], kep_orbit['i'][0])
        u.isclose(orbit['Omega'][0], kep_orbit['Omega'][0])
        u.isclose(orbit['w'][0], kep_orbit['w'][0])
        u.isclose(orbit['M'][0], kep_orbit['M'][0])
        u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

    def test_timescales(self):

        orbit = Orbit.from_horizons('Ceres')
        orbit['epoch'] = orbit['epoch'].tdb

        cart_orbit = orbit.oo_transform('CART')

        kep_orbit = cart_orbit.oo_transform('KEP')
        u.isclose(orbit['a'][0], kep_orbit['a'][0])
        u.isclose(orbit['e'][0], kep_orbit['e'][0])
        u.isclose(orbit['i'][0], kep_orbit['i'][0])
        u.isclose(orbit['Omega'][0], kep_orbit['Omega'][0])
        u.isclose(orbit['w'][0], kep_orbit['w'][0])
        u.isclose(orbit['M'][0], kep_orbit['M'][0])
        u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

        com_orbit = orbit.oo_transform('COM')
        kep_orbit = com_orbit.oo_transform('KEP')
        u.isclose(orbit['a'][0], kep_orbit['a'][0])
        u.isclose(orbit['e'][0], kep_orbit['e'][0])
        u.isclose(orbit['i'][0], kep_orbit['i'][0])
        u.isclose(orbit['Omega'][0], kep_orbit['Omega'][0])
        u.isclose(orbit['w'][0], kep_orbit['w'][0])
        u.isclose(orbit['M'][0], kep_orbit['M'][0])
        u.isclose(orbit['epoch'][0].utc.jd, kep_orbit['epoch'][0].utc.jd)

        assert kep_orbit['epoch'].scale == 'tdb'


@pytest.mark.remote_data
class TestOOPropagate:

    def test_oo_propagate(self):
        """ test oo_propagate method"""

        try:
            import pyoorb
        except ImportError:
            return None

        orbit = Orbit.from_horizons('Ceres')
        epoch = Time(Time.now().jd + 100, format='jd', scale='utc')

        future_orbit = Orbit.from_horizons('Ceres',
                                           epochs=Time(epoch, format='jd',
                                                       scale='utc'))

        oo_orbit = orbit.oo_propagate(epoch)

        u.isclose(oo_orbit['a'][0], future_orbit['a'][0])
        u.isclose(oo_orbit['e'][0], future_orbit['e'][0])
        u.isclose(oo_orbit['i'][0], future_orbit['i'][0])
        u.isclose(oo_orbit['Omega'][0], future_orbit['Omega'][0])
        u.isclose(oo_orbit['w'][0], future_orbit['w'][0])
        u.isclose(oo_orbit['M'][0], future_orbit['M'][0])
        u.isclose(oo_orbit['epoch'][0].utc.jd,
                  future_orbit['epoch'][0].utc.jd)

        assert oo_orbit['epoch'].scale == 'utc'
