# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..core import *


def test_predefined_reinitialisation():
    assert u.mag('VegaFluxDensity') == VegaMag


def test_predefined_string_roundtrip():
    assert u.Unit(VegaFluxDensity.to_string()) == VegaFluxDensity
    assert u.Unit(VegaMag.to_string()) == VegaMag
    assert u.Unit(hundred_nm.to_string()) == hundred_nm
