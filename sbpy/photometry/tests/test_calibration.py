# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import pytest
import mock
import importlib
import astropy.units as u
import numpy as np
import synphot
from ...units import VegaMag
from ...spectroscopy.vega import Vega
from ...spectroscopy.sun import Sun


def test_synphot_import_fail():
    # importlib rigamarole is supposedly not needed, but this didn't work without it.
    with mock.patch.dict(sys.modules, {'synphot': None}):
        from .. import calibration
        importlib.reload(calibration)
        assert calibration.synphot is None
    importlib.reload(calibration)


@pytest.mark.parametrize('unit, m0', (
    ('W/(m2 um)', -27.04 * VegaMag),
    ('W/(m2 Hz)', -26.93 * u.ABmag),
    ('W/(m2 um)', -26.66 * u.STmag)
))
def test_fluxd_to_mag_sun(unit, m0):
    """Compare sbpy's Sun apparent magnitudes to Willmer 2018.

    Tests SDSS r bandpass.

    """
    from ..calibration import fluxd_to_mag

    sun_fn = os.path.join(os.path.dirname(__file__),
                          'data', 'sun-r-haberreiter17.txt')
    sun = Sun.from_file(sun_fn, wave_unit=u.AA, flux_unit='W/(m2 um)')

    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    wave, fluxd = sun.filt(bp, unit=unit)

    m = fluxd_to_mag(fluxd, m0.unit, bandpass=bp)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit

    # also test wave
    m = fluxd_to_mag(fluxd, m0.unit, wave=wave)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit


@pytest.mark.parametrize('unit, m0', (
    ('W/(m2 Hz)', 0.03 * VegaMag),
    ('W/(m2 Hz)', 0.149 * u.ABmag),
    ('W/(m2 um)', 0.410 * u.STmag)
))
def test_fluxd_to_mag_vega(unit, m0):
    """Compare sbpy's Vega apparent magnitudes to Willmer 2018.

    Note: AB and ST mags in Table 3 of Willmer 2018 are for conversion
    to/from Vegamag, i.e., they include a 0.03 mag offset from the
    spectrum of Vega: AB(Vega) = 0.119, ST(Vega) = 0.380.

    Unlike the Sun, fluxd(r bandpass) and fluxd(r effective
    wavelength) are not similar.

    """
    from ..calibration import fluxd_to_mag

    vega = Vega.from_default()
    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    wave, fluxd = vega.filt(bp, unit=unit)
    print(wave, fluxd)

    m = fluxd_to_mag(fluxd, m0.unit, bandpass=bp)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit


def test_fluxd_to_mag_bp_wave_error():
    from ..calibration import fluxd_to_mag
    with pytest.raises(ValueError):
        fluxd_to_mag(1 * u.Jy, VegaMag)


def test_fluxd_to_mag_mag_unit_error():
    from ..calibration import fluxd_to_mag
    with pytest.raises(ValueError):
        fluxd_to_mag(1 * u.Jy, u.mag)


def test_fluxd_to_mag_fluxd_unit_error():
    from ..calibration import fluxd_to_mag
    with pytest.raises(ValueError):
        fluxd_to_mag(1 * u.m, VegaMag, wave=1 * u.um)


@pytest.mark.parametrize('m, unit', (
    (VegaMag, VegaMag),
    (synphot.units.VEGAMAG, VegaMag),
    (u.ABmag, u.ABmag),
    (u.STmag, u.STmag),
))
def test_validate_mag(m, unit):
    from ..calibration import validate_mag
    assert validate_mag(m) == unit
    if not isinstance(m, str):
        assert validate_mag(0 * m).unit == unit


def test_validate_mag_not_magnitude():
    from ..calibration import validate_mag
    with pytest.raises(ValueError):
        validate_mag(0 * u.m)
