# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import mock
import pytest
import numpy as np
from ..bandpass import *


@pytest.mark.parametrize('name, avgwave', (
    ('2mass j', 12410.52630476),
    ('2mass h', 16513.66475736),
    ('2mass ks', 21656.32078498),
    ('cousins r', 6499.91478190),
    ('cousins i', 7884.10581303),
    ('johnson u', 3598.54452094),
    ('johnson b', 4385.9244053),
    ('johnson v', 5490.55520036),
    ('ps1 g', 486.64578708),
    ('ps1 r', 621.4623038),
    ('ps1 i', 754.4570357),
    ('ps1 w', 638.93518241),
    ('ps1 y', 963.32481028),
    ('ps1 z', 867.946803),
    ('sdss u', 3561.78873418),
    ('sdss g', 4718.87224631),
    ('sdss r', 6185.19447698),
    ('sdss i', 7499.70417489),
    ('sdss z', 8961.48833667),
    ('wfc3 f438w', 4324.44438089),
    ('wfc3 f606w', 5946.7429129),
    ('wise w1', 34002.59750555),
    ('wise w2', 46520.16384937),
    ('wise w3', 128108.72187073),
    ('wise w4', 223752.74423983),
))
def test_bandpass(name, avgwave):
    bp = bandpass(name)
    assert np.isclose(bp.avgwave().value, avgwave)


def test_bandpass_synphot():
    with mock.patch.dict(sys.modules, {'synphot': None}):
        with pytest.raises(ImportError):
            bandpass('sdss u')
