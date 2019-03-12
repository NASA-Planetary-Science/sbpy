# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from ..core import *


@pytest.mark.parametrize('name, avgwave', (
    ('f438w', 4324.44438089),
    ('f606w', 5946.7429129),
    ('cousins i', 7884.10581303),
    ('sdss u', 3561.78873418),
    ('sdss g', 4718.87224631),
    ('sdss r', 6185.19447698),
    ('sdss i', 7499.70417489),
    ('sdss z', 8961.48833667)
))
def test_get_bandpass(name, avgwave):
    bp = get_bandpass(name)
    assert np.isclose(bp.avgwave().value, avgwave)
