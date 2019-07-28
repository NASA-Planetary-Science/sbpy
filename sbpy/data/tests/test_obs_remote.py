# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from numpy.testing import assert_allclose

from sbpy.data import Obs
from sbpy import bib


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

        # comet
        data = Obs.from_mpc('235P', id_type='comet number')
        assert len(data) >= 274

        data = Obs.from_mpc('P/2010 F2', id_type='comet designation')
        assert len(data) >= 35

    def test_bib(self):
        bib.track()
        Obs.from_mpc('235P')
        assert 'sbpy.data.obs.from_mpc' in bib.to_text()
