# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from sbpy.data import Phys
from sbpy import bib


@pytest.mark.remote_data
def test_from_sbdb():
    """ test from_horizons method"""

    # query one object
    data = Phys.from_sbdb('Ceres')
    assert len(data.table) == 1

    # query several objects
    with bib.Tracking():
        data = Phys.from_sbdb([n+1 for n in range(5)])
    assert len(data.table) == 5

    assert data['H'].unit == 'mag'
    assert data['G'].unit is None

    # assert 'sbpy.data.Phys.from_sbdb' in bib.to_text()
