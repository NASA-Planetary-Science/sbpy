# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

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


@pytest.mark.parametrize(
    "a,b,theta",
    [
        ([1, 0, 0], [1, 1, 0], 45),
        ([1, 0, 0], [1, -1, 0], 45),
        ([1, 1, 0], [1, -1, 0], 90),
        ([1, 0, 0], [1 / np.sqrt(3), 1, 0], 60),
        ([1, 0, 0], [np.sqrt(3), -1, 0], 30),
        ([1 / np.sqrt(3), 1, 0], [np.sqrt(3), -1, 0], 90),
        ([1, 0, 0], [1 / np.sqrt(3), -1, 0], 60),
        ([1 / np.sqrt(3), 1, 0], [1 / np.sqrt(3), -1, 0], 120),
    ],
)
def test_angle(a, b, theta):
    assert u.isclose(Surface._angle(a, b * u.au), theta * u.deg)
