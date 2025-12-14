# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy.coordinates import Angle
from astropy import units as u

from ..surface import Surface, min_zero_cos


def test_min_zero_cos():
    a = Angle([-91, -90, 0, 30, 45, 60, 90, 91], "deg")
    result = min_zero_cos(a)
    expected = [0, 0, 1, np.sqrt(3) / 2, 1 / np.sqrt(2), 0.5, 0, 0]
    assert np.allclose(result, expected)

    # test scalars
    for i in range(len(a)):
        assert np.isclose(min_zero_cos(a[i]), expected[i])


def test_vectors_to_angles():
    n = [1, 0, 0]
    rs = [1, 1, 0] * u.au
    ro = [1, -1, 0] * u.au
    i, e, phi = Surface._vectors_to_angles(n, rs, ro)
    assert u.isclose(i, 45 * u.deg)
    assert u.isclose(e, 45 * u.deg)
    assert u.isclose(phi, 90 * u.deg)

    n = [1, 0, 0]
    rs = [1 / np.sqrt(3), 1, 0] * u.au
    ro = [np.sqrt(3), -1, 0] * u.au
    i, e, phi = Surface._vectors_to_angles(n, rs, ro)
    assert u.isclose(i, 60 * u.deg)
    assert u.isclose(e, 30 * u.deg)
    assert u.isclose(phi, 90 * u.deg)

    ro = [1 / np.sqrt(3), -1, 0] * u.au
    i, e, phi = Surface._vectors_to_angles(n, rs, ro)
    assert u.isclose(e, 60 * u.deg)
    assert u.isclose(phi, 120 * u.deg)
