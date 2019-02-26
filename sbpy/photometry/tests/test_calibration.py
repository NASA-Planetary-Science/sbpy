# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import pytest
import astropy.units as u
import numpy as np
import synphot
from ..calibration import *
from ..calibration import validate_mag
from ...units import VegaMag
from ...spectroscopy.vega import Vega
from ...spectroscopy.sun import Sun


@pytest.mark.parametrize('standard, unit, m0', (
    (Sun.from_default(), 'W/(m2 um)', -27.04 * VegaMag),
    (Sun.from_default(), 'W/(m2 Hz)', -26.93 * u.ABmag),
    (Sun.from_default(), 'W/(m2 um)', -26.66 * u.STmag)
))
def test_fluxd_to_mag_sun(standard, unit, m0):
    """Compare sbpy's Sun apparent magnitudes to Willmer 2018.

    Test r bandpass and r effective wavelength, since the two are
    fairly close.

    """
    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    wave, fluxd = standard.filt(bp, unit=unit)

    m = fluxd_to_mag(fluxd, m0.unit, bandpass=bp)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit

    fluxd = standard(wave, unit=unit)
    m = fluxd_to_mag(fluxd, m0.unit, wave=wave)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit


@pytest.mark.parametrize('standard, unit, m0', (
    (Vega.from_default(), 'W/(m2 Hz)', 0.03 * VegaMag),
    (Vega.from_default(), 'W/(m2 Hz)', 0.149 * u.ABmag),
    (Vega.from_default(), 'W/(m2 um)', 0.410 * u.STmag)
))
def test_fluxd_to_mag_vega(standard, unit, m0):
    """Compare sbpy's Vega apparent magnitudes to Willmer 2018.

    Note: AB and ST mags in Table 3 of Willmer 2018 are for conversion
    to/from Vegamag, i.e., they include a 0.03 mag offset from the
    spectrum of Vega: AB(Vega) = 0.119, ST(Vega) = 0.380.

    Unlike the Sun, fluxd(r bandpass) and fluxd(r effective
    wavelength) are not similar.

    """

    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    wave, fluxd = standard.filt(bp, unit=unit)
    print(wave, fluxd)

    m = fluxd_to_mag(fluxd, m0.unit, bandpass=bp)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit
