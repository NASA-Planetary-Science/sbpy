# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from astropy import units as u

from ..transformations import twovec, xyz2sph, sph2xyz


def test_twovec():
    x = [1, 0, 0] * u.au
    y = [0, 1, 0] * u.au
    z = [0, 0, 1] * u.au

    # x -> y'
    # y -> z'
    # z -> x'
    M = twovec(z, 0, y, 2)
    assert u.allclose(np.dot(M, x), y)
    assert u.allclose(np.dot(M, y), z)
    assert u.allclose(np.dot(M, z), x)

    # x -> y'
    # y -> -x'
    # z -> z'
    M = twovec(x, 1, z, 2)
    assert u.allclose(np.dot(M, x), y)
    assert u.allclose(np.dot(M, y), -x)
    assert u.allclose(np.dot(M, z), z)

    with pytest.raises(ValueError):
        M = twovec(z, 0, 2 * z, 1)


def test_xyz2sph():
    x = [1, 0, 0] * u.au
    y = [0, 1, 0] * u.au
    z = [0, 0, 1] * u.au

    assert u.allclose(xyz2sph(*x), [0, 0] * u.deg)
    assert u.allclose(xyz2sph(*y), [90, 0] * u.deg)
    assert u.allclose(xyz2sph(*z), [0, 90] * u.deg)


def test_sph2xyz():
    x = [1, 0, 0] * u.au
    y = [0, 1, 0] * u.au
    z = [0, 0, 1] * u.au

    assert u.allclose(sph2xyz(0 * u.deg, 0 * u.deg), x.value)
    assert u.allclose(sph2xyz(90 * u.deg, 0 * u.deg), y.value, atol=1e-16)
    assert u.allclose(sph2xyz(0 * u.deg, 90 * u.deg), z.value, atol=1e-16)

    r = 1 * u.au
    assert u.allclose(sph2xyz(0 * u.deg, 0 * u.deg, r), x)
    assert u.allclose(sph2xyz(90 * u.deg, 0 * u.deg, r), y, atol=1e-16 * u.au)
    assert u.allclose(sph2xyz(0 * u.deg, 90 * u.deg, r), z, atol=1e-16 * u.au)
