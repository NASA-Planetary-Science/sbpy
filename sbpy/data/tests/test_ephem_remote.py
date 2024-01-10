# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from numpy import abs
import warnings

from numpy.testing import assert_allclose
import astropy.units as u
from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.tests.helper import assert_quantity_allclose

from ... import bib
from .. import Ephem, Orbit, QueryError

pytest.importorskip("astroquery")

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


@pytest.mark.remote_data
class TestEphemFromHorizons:

    def test_current_epoch(self):
        now = Time.now()
        data = Ephem.from_horizons('Ceres')
        assert_allclose(data['epoch'].jd, now.jd)

    def test_daterange_step(self):
        """Test start/stop/step.

        Also, verify that the epochs parameter was not modified in-place,
        sbpy PR#247.

        """
        epochs = {'start': Time('2018-01-02', format='iso'),
                  'stop': Time('2018-01-05', format='iso'),
                  'step': 6*u.hour}
        epochs_saved = epochs.copy()

        data = Ephem.from_horizons('Ceres', epochs=epochs)

        assert len(data.table) == 13
        assert (epochs == epochs_saved) and isinstance(epochs['start'], Time)

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

    def test_epochs_order(self):
        """sbpy GitHub issue #317: Horizons returns sorted dates."""
        epochs = Time(['2010-01-02', '2010-01-01'])
        assert epochs[0] > epochs[1]
        data = Ephem.from_horizons('Ceres', epochs=epochs)
        # verify that from_horizons returns the dates in the order they were
        # requested
        assert_allclose(epochs.mjd, data['date'].mjd)

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
        with bib.Tracking():
            Ephem.from_horizons(['Ceres', 'Pallas'])
            assert 'sbpy.data.ephem.Ephem.from_horizons' in bib.show()
        bib.reset()

    def test_timescale(self):
        # test same timescale
        time = Time('2020-01-01 00:00', format='iso', scale='utc')
        a = Ephem.from_horizons(1, epochs=time)
        assert a['epoch'].scale == 'utc'
        assert a['epoch'].utc.jd == 2458849.5

        # test different timescale (and check for warning)
        time = Time('2020-01-01 00:00', format='iso', scale='tdb')
        with warnings.catch_warnings(record=True) as w:
            a = Ephem.from_horizons(1, epochs=time)
            assert any(["astroquery.jplhorizons" in str(w[i].message)
                        for i in range(len(w))])
        assert a['epoch'].scale == 'utc'
        assert_allclose(a['epoch'][0].utc.jd, 2458849.49919926)

    def test_queryfail(self):
        with pytest.raises(QueryError):
            Ephem.from_horizons('target does not exist')

    def test_output_epochs(self):
        from astropy.utils.iers import IERSRangeError

        # transformed from ut1 to utc
        # this test sometimes causes trouble in remote tests, potentially
        # because it has to download some data, which might fail;
        # if that happens, skip the test
        try:
            epochs = Time(2437655.5, format='jd', scale='utc')
            a = Ephem.from_horizons(1, epochs=epochs)
            assert a['epoch'][0].scale == 'utc'
        except IERSRangeError:
            pass

        # all utc
        epochs = Time(2437675.5, format='jd', scale='utc')
        a = Ephem.from_horizons(1, epochs=epochs)
        assert a['epoch'][0].scale == 'utc'

        # mixed ut1 and utc
        # same as above
        try:
            epochs = Time([2437655.5, 2437675.5], format='jd', scale='utc')
            a = Ephem.from_horizons(1, epochs=epochs)
            assert a['epoch'].scale == 'utc'
        except IERSRangeError:
            pass

    def test_minimum_working_columns(self):
        """Avoid assuming any particular columns are returned.

        Regression test for #242.

        """
        e = Ephem.from_horizons(1, quantities='1')  # just RA, Dec
        assert len(e) == 1
        e = Ephem.from_horizons(1, quantities='19')  # just rh and delta-rh
        assert len(e) == 1

    def test_invalid_epochs(self):
        with pytest.raises(ValueError):
            Ephem.from_horizons('Ceres', epochs='2020-04-01')


@pytest.mark.remote_data
class TestEphemFromMPC:
    def test_single_epoch_now(self):
        eph = Ephem.from_mpc('Ceres')
        assert len(eph.table) == 1

    def test_single_epoch(self):
        # utc
        eph1 = Ephem.from_mpc('Ceres',
                              epochs=Time('2018-10-01', scale='utc'))
        assert len(eph1.table) == 1

        # tdb
        with warnings.catch_warnings(record=True) as w:
            eph2 = Ephem.from_mpc('Ceres',
                                  epochs=Time('2018-10-01', scale='tdb'))
            assert any(["astroquery.mpc" in str(w[i].message)
                        for i in range(len(w))])
        assert all(eph2['epoch'] != eph1['epoch'])

    def test_multiple_epochs(self):
        # utc
        eph1 = Ephem.from_mpc('Ceres', epochs=Time(['2018-10-01',
                                                    '2019-10-01'],
                                                   scale='utc'))
        assert len(eph1.table) == 2

        # tai
        with warnings.catch_warnings(record=True) as w:
            eph2 = Ephem.from_mpc('Ceres', epochs=Time(['2018-10-01',
                                                        '2019-10-01'],
                                                       scale='tai'))
            assert any(["astroquery.mpc" in str(w[i].message)
                        for i in range(len(w))])
        assert all(eph2['epoch'] != eph1['epoch'])

    def test_invalid_epochs(self):
        with pytest.raises(ValueError):
            Ephem.from_mpc('Ceres', epochs='2020-04-01')

    def test_start_stop_step(self):
        """Test start/stop/step.

        Also, verify that the epochs parameter was not modified in-place,
        sbpy PR#247.

        """

        # utc
        epochs = dict(start=Time('2018-10-01'), stop=Time('2018-10-31'),
                      step=1*u.d)
        epochs_saved = epochs.copy()
        eph1 = Ephem.from_mpc('Ceres', epochs=epochs)
        assert len(eph1.table) == 31
        assert (epochs == epochs_saved) and isinstance(epochs['start'], Time)

        # tt
        epochs = dict(start=Time('2018-10-01', scale='tt'),
                      stop=Time('2018-10-31', scale='tt'),
                      step=1*u.d)
        with warnings.catch_warnings(record=True) as w:
            eph2 = Ephem.from_mpc('Ceres', epochs=epochs)
            assert any(["astroquery.mpc" in str(w[i].message)
                        for i in range(len(w))])
        assert all(eph2['epoch'] != eph1['epoch'])

    def test_minute_steps_pr88(self):
        """https://github.com/NASA-Planetary-Science/sbpy/pull/88"""
        epochs = dict(start=Time('2018-10-01'), step=1*u.min, number=10)
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert len(eph.table) == 10

    def test_start_stop_no_step(self):
        with pytest.raises(QueryError):
            eph = Ephem.from_mpc('Ceres', epochs={
                'start': Time('2018-10-01'),
                'stop': Time('2018-10-31')})

    def test_start_step_number(self):
        epochs = dict(start=Time('2018-10-01'), step=1*u.d, number=31)
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert len(eph.table) == 31
        assert eph['Date'][-1] == '2018-10-31 00:00:00.000'

    def test_start_stop_number(self):
        epochs = dict(start=Time('2018-10-01'), stop=Time('2018-10-02'),
                      number=31)
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert len(eph.table) == 31

    def test_start_stop_jd(self):
        epochs = {'start': Time(2458396.5, format='jd'),
                  'stop': Time(2458397.5, format='jd'),
                  'step': 1*u.d}
        eph = Ephem.from_mpc('Ceres', epochs=epochs)
        assert eph['Date'][0] == '2018-10-05 00:00:00.000'
        assert eph['Date'][1] == '2018-10-06 00:00:00.000'

    def test_step_unit(self):
        with pytest.raises(QueryError):
            Ephem.from_mpc('Ceres', epochs={
                'start': Time('2018-10-01'),
                'step': 1*u.year})

    def test_ra_dec_format(self):
        epochs = dict(start=Time('2018-10-01'),
                      step=1*u.d, number=31)
        ra_format = {'sep': ':', 'unit': 'hourangle', 'precision': 1}
        dec_format = {'sep': ':', 'precision': 1}
        eph = Ephem.from_mpc('Ceres', epochs=epochs, ra_format=ra_format,
                             dec_format=dec_format)

        # round-trip through a string, then inspect RA, Dec
        tab = ascii.read(
            '\n'.join(eph.table.pformat(
                show_unit=False,
                max_lines=-1,
                max_width=-1
            )),
            format='fixed_width_two_line'
        )
        assert isinstance(tab['RA'][0], str)
        assert isinstance(tab['Dec'][0], str)

    def test_multiple_targets(self):
        eph = Ephem.from_mpc(['Ceres', 'Pallas', 'Vesta'])
        assert all(eph['Targetname'].data == ['Ceres', 'Pallas', 'Vesta'])

    def test_queryfail(self):
        with pytest.raises(QueryError):
            Ephem.from_mpc('target does not exist')

    def test_bib(self):
        with bib.Tracking():
            Ephem.from_mpc(['Ceres', 'Pallas'])
            assert 'sbpy.data.ephem.Ephem.from_mpc' in bib.show()
        bib.reset()


@pytest.mark.remote_data
class TestEphemFromMiriade:
    def test_singleobj_now(self):
        eph = Ephem.from_miriade("Ceres")
        assert_quantity_allclose(eph['epoch'][0].jd, Time.now().jd)

    def test_multiobj_now(self):
        eph = Ephem.from_miriade(["Ceres", 'Pallas'],
                                 epochs=Time.now())
        assert len(eph.table) == 2

    def test_epochTime(self):
        # utc
        eph1 = Ephem.from_miriade(["Ceres", 'Pallas'],
                                  epochs=Time('2018-10-01', scale='utc'))
        assert_quantity_allclose(eph1['epoch'][0].jd, Time('2018-10-01').jd)

        # tt
        with warnings.catch_warnings(record=True) as w:
            eph2 = Ephem.from_miriade(["Ceres", 'Pallas'],
                                      epochs=Time('2018-10-01', scale='tt'))
            assert any(["astroquery.imcce" in str(w[i].message)
                        for i in range(len(w))])
        assert all(eph1['epoch'] != eph2['epoch'])

    def test_epochTimemult(self):
        eph = Ephem.from_miriade(["Ceres", 'Pallas'],
                                 epochs=Time(['2018-10-01', '2018-10-02']))
        assert_quantity_allclose(eph['epoch'][0].jd, Time('2018-10-01').jd)

    def test_epochstart(self):
        eph = Ephem.from_miriade(["Ceres", 'Pallas'],
                                 epochs={'start': Time('2018-10-01')})
        assert_quantity_allclose(eph['epoch'][0].jd, Time('2018-10-01').jd)

    def test_epochrange_number(self):
        # utc
        eph1 = Ephem.from_miriade(["Ceres", 'Pallas'],
                                  epochs={'start': Time('2018-10-01'),
                                          'stop': Time('2018-11-01'),
                                          'number': 10})
        assert len(eph1.table) == 20

        # tdb
        with warnings.catch_warnings(record=True) as w:
            eph2 = Ephem.from_miriade(["Ceres", 'Pallas'],
                                      epochs={'start': Time('2018-10-01',
                                                            scale='tdb'),
                                              'stop': Time('2018-11-01',
                                                           scale='tdb'),
                                              'number': 10})
            assert any(["astroquery.imcce" in str(w[i].message)
                        for i in range(len(w))])
        assert all(eph1['epoch'] != eph2['epoch'])

    def test_epochrange_step(self):
        """Test start/stop/step.

        Also, verify that the epochs parameter was not modified in-place,
        sbpy PR#247.

        """
        epochs = {
            'start': Time('2018-10-01'),
            'stop': Time('2018-10-10'),
            'step': 0.5*u.d
        }
        epochs_saved = epochs.copy()
        eph = Ephem.from_miriade(["Ceres", 'Pallas'], epochs=epochs)
        assert len(eph.table) == 38
        assert (epochs == epochs_saved) and isinstance(epochs['start'], Time)

    def test_invalid_epochs(self):
        with pytest.raises(ValueError):
            Ephem.from_miriade('Ceres', epochs='2020-04-01')

    def test_iaulocation(self):
        eph1 = Ephem.from_miriade(
            "Ceres", location='G37', epochs=Time('2019-01-01'))
        eph2 = Ephem.from_miriade("Ceres", epochs=Time('2019-01-01'))
        assert abs(eph1['RA'][0]-eph2['RA'][0]) > 0.00001*u.deg

    def test_EarthLocation(self):
        lowell = EarthLocation.of_site('Lowell Observatory')
        eph1 = Ephem.from_miriade(
            "Ceres", location=lowell, epochs=Time('2019-01-01'))
        eph2 = Ephem.from_miriade("Ceres", epochs=Time('2019-01-01'))
        assert abs(eph1['RA'][0]-eph2['RA'][0]) > 0.0001*u.deg

    def test_queryfail(self):
        with pytest.raises(QueryError):
            Ephem.from_miriade('target does not exist')

    def test_bib(self):
        with bib.Tracking():
            data = Ephem.from_miriade(['Ceres', 'Pallas'])
            assert 'sbpy.data.ephem.Ephem.from_miriade' in bib.show()
        bib.reset()


@pytest.mark.remote_data
class TestEphemFromOorb:
    def test_by_comparison(self):
        """test from_oo method"""

        pytest.importorskip("pyoorb")

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
