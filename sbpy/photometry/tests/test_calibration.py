# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import pytest
import astropy.units as u
import synphot
from ..calibration import *
from ..calibration import validate_mag
from ...units import VegaMag
from ...spectroscopy.vega import Vega
from ...spectroscopy.sun import Sun


@pytest.mark.parametrize('standard, unit, m0', (
    (Vega.from_default(), 'W/(m2 Hz)', 0.119 * u.ABmag),
    (Vega.from_default(), 'W/(m2 um)', 0.380 * u.STmag),
    (Sun.from_default(), 'W/(m2 um)', -27.04 * VegaMag),
    (Sun.from_default(), 'W/(m2 Hz)', -26.93 * u.ABmag),
    (Sun.from_default(), 'W/(m2 um)', -26.66 * u.STmag)
))
def test_fluxd_to_mag(standard, unit, m0):
    """Compare sbpy's Sun and Vega apparent magnitudes to Willmer 2018."""
    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    wave, fluxd = standard.filt(bp, unit=unit)

    m = fluxd_to_mag(fluxd, m0.unit, bandpass=bp)
    assert u.isclose(m, m0)

    m = fluxd_to_mag(fluxd, m0.unit, wave=wave)
    assert u.isclose(m, m0)
