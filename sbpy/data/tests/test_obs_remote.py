# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.time import Time
import astropy.units as u

from .. import Obs
from ... import bib
from ..core import QueryError

pytest.importorskip("astroquery")


@pytest.mark.remote_data
class TestObsfromMPC:

    def test_simple(self):
        # asteroid
        data = Obs.from_mpc('12893')
        assert len(data) >= 1605

        # comet
        data = Obs.from_mpc('235P')
        assert len(data) >= 274

    def test_manual_idtypes(self):
        data = Obs.from_mpc('12893', id_type='asteroid number')
        assert len(data) >= 1605

        data = Obs.from_mpc('1998 QS55', id_type='asteroid designation')
        assert len(data) >= 46

        data = Obs.from_mpc('2019 AA', id_type='asteroid designation')
        assert len(data) >= 33

        # comet
        data = Obs.from_mpc('235P', id_type='comet number')
        assert len(data) >= 274

        data = Obs.from_mpc('P/2010 F2', id_type='comet designation')
        assert len(data) >= 35

    def test_break(self):
        with pytest.raises(QueryError):
            Obs.from_mpc('2019 AA345', id_type='asteroid designation')

    def test_bib(self):
        bib.track()
        Obs.from_mpc('235P')
        assert 'sbpy.data.obs.Obs.from_mpc' in bib.to_text()


@pytest.mark.remote_data
class TestSupplement:

    def test_jplhorizons(self):
        bib.track()
        obs = Obs.from_dict({'epoch': Time([2451200, 2451201], format='jd'),
                             'mag': [12, 13]*u.mag,
                             'targetname': ['3552', '3552']})
        data = obs.supplement(service='jplhorizons', modify_fieldnames='obs')
        assert len(data.field_names) > len(obs.field_names)
        assert 'targetname_obs' in data.field_names

        data = obs.supplement(service='jplhorizons', modify_fieldnames='eph')
        assert len(data.field_names) > len(obs.field_names)
        assert 'targetname_eph' in data.field_names

        assert 'sbpy.data.ephem.Ephem.from_horizons' in bib.to_text()

    def test_mpc(self):
        bib.track()
        obs = Obs.from_dict({'epoch': Time([2451200, 2451201], format='jd'),
                             'mag': [12, 13]*u.mag,
                             'targetname': ['3552', '3552']})
        data = obs.supplement(service='mpc')
        assert len(data.field_names) > len(obs.field_names)

        assert 'sbpy.data.ephem.Ephem.from_mpc' in bib.to_text()

    def test_miriade(self):
        bib.track()
        obs = Obs.from_dict({'epoch': Time([2451200, 2451201], format='jd'),
                             'mag': [12, 13]*u.mag,
                             'targetname': ['3552', '3552']})
        data = obs.supplement(service='miriade')
        assert len(data.field_names) > len(obs.field_names)

        assert 'sbpy.data.ephem.Ephem.from_miriade' in bib.to_text()

    def test_breaks(self):
        obs = Obs.from_dict({'epoch': Time([2451200, 2451201], format='jd'),
                             'mag': [12, 13]*u.mag,
                             'targetname': ['3552', '3552']})

        with pytest.raises(QueryError):
            obs.supplement(service='this will not work')

    def test_multiple_jplhorizons(self):
        obs = Obs.from_dict({'epoch': Time([2451200, 2451201,
                                            2451200, 2451201], format='jd'),
                             'mag': [12, 13, 16, 17]*u.mag,
                             'targetname': ['3552', '3552',
                                            '12893', '12893']})
        data = obs.supplement(service='jplhorizons', modify_fieldnames='obs')
        assert len(set(data['targetname'])) == 2
