# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import astropy.units as u
from sbpy.data import Phys

pytest.importorskip("astroquery")


@pytest.mark.remote_data
def test_from_sbdb():
    """ test from_sbdb method"""

    # query one object
    data = Phys.from_sbdb('Ceres')
    assert len(data.table) == 1

    # query several objects
    data = Phys.from_sbdb([n+1 for n in range(5)])
    assert len(data.table) == 5


@pytest.mark.remote_data
def test_from_sbdb_comet():
    """Regression test for issues #349 and #358.

    As of June 2022, astroquery does not assign units to M1, M2 and their
    uncertainties, nor H.

    """
    # need a comet with both M1 and M2, and their uncertainties:
    data = Phys.from_sbdb('147P')
    for k in ('M1', 'M2', 'M1_sig', 'M2_sig'):
        assert isinstance(data[k], u.Quantity) and data[k].unit == u.mag

    data = Phys.from_sbdb('1')
    assert isinstance(data['H'], u.Quantity) and data['H'].unit == u.mag
