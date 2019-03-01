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
from .. import calibration as cal


def test_synphot_import_fail():
    # importlib rigamarole is supposedly not needed, but this
    # didn't work without it. MSK / 2019 Mar 01 / Python 3.6.7
    with mock.patch.dict(sys.modules, {'synphot': None}):
        from .. import calibration
        importlib.reload(calibration)
        assert calibration.synphot is None
    importlib.reload(calibration)


@pytest.mark.parametrize('m, unit', (
    (VegaMag, VegaMag),
    (synphot.units.VEGAMAG, VegaMag),
    (u.ABmag, u.ABmag),
    (u.STmag, u.STmag),
))
def test_validate_mag(m, unit):
    assert cal.validate_mag(m) == unit
    if not isinstance(m, str):
        assert cal.validate_mag(0 * m).unit == unit


def test_validate_mag_not_magnitude():
    with pytest.raises(ValueError):
        cal.validate_mag(0 * u.Jy)


@pytest.mark.parametrize('unit, m0', (
    ('W/(m2 um)', -27.04 * VegaMag),
    ('W/(m2 Hz)', -26.93 * u.ABmag),
    ('W/(m2 um)', -26.66 * u.STmag)
))
def test_fluxd_to_mag_sun(unit, m0):
    """Compare sbpy's Sun apparent magnitudes to Willmer 2018.

    Tests SDSS r bandpass and solar spectrum of Haberreiter et
    al. 2017, JGR Space Physics, 122, 5910-5930, 10.1002/2016JA023492.

    """

    sun_fn = os.path.join(os.path.dirname(__file__),
                          'data', 'sun-r-haberreiter17.txt')
    sun = Sun.from_file(sun_fn, wave_unit=u.AA, flux_unit='W/(m2 um)')

    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    wave, fluxd = sun.filt(bp, unit=unit)

    m = cal.fluxd_to_mag(fluxd, m0.unit, bandpass=bp)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit

    # also test wavelength
    m = cal.fluxd_to_mag(fluxd, m0.unit, wave=wave)
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
    wavelength) are not similar, so this does not test wave keyword.

    """
    vega = Vega.from_default()
    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    wave, fluxd = vega.filt(bp, unit=unit)
    print(wave, fluxd)

    m = cal.fluxd_to_mag(fluxd, m0.unit, bandpass=bp)
    assert np.isclose(m.value, m0.value, atol=0.01)
    assert m.unit == m0.unit


def test_fluxd_to_mag_bp_wave_error():
    with pytest.raises(ValueError):
        cal.fluxd_to_mag(1 * u.Jy, VegaMag)


def test_fluxd_to_mag_fluxd_unit_error():
    with pytest.raises(ValueError):
        cal.fluxd_to_mag(1 * u.m, VegaMag, wave=1 * u.um)


@pytest.mark.parametrize('m', (
    (-27.04 * VegaMag),
    (-26.93 * u.ABmag),
    (-26.66 * u.STmag)
))
def test_mag_to_fluxd_sun(m):
    """Inverse of test_fluxd_to_mag_sun."""

    sun_fn = os.path.join(os.path.dirname(__file__),
                          'data', 'sun-r-haberreiter17.txt')
    sun = Sun.from_file(sun_fn, wave_unit=u.AA, flux_unit='W/(m2 um)')

    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    fluxd0 = 1674.943 * u.Unit('W/(m2 um)')

    fluxd = cal.mag_to_fluxd(m, fluxd0.unit, bandpass=bp)
    assert np.isclose(fluxd.value, fluxd0.value, rtol=0.01)
    assert fluxd.unit == fluxd0.unit

    # also test wavelength
    fluxd = cal.mag_to_fluxd(m, fluxd0.unit, wave=bp.pivot())
    assert np.isclose(fluxd.value, fluxd0.value, rtol=0.01)
    assert fluxd.unit == fluxd0.unit


def test_mag_to_fluxd_bp_wave_error():
    with pytest.raises(ValueError):
        cal.mag_to_fluxd(0 * VegaMag, u.Jy)


def test_mag_to_fluxd_fluxd_unit_error():
    with pytest.raises(ValueError):
        cal.mag_to_fluxd(1 * VegaMag, u.m, wave=1 * u.um)


@pytest.mark.parametrize('from_value, to_unit, to_value', (
    (0 * VegaMag, u.Jy, 3229.07),
    (3630.78 * u.Jy, u.ABmag, 1.6e-7),
    (0 * u.ABmag, u.STmag, 0.26133661),
    (0 * u.STmag, u.ABmag, -0.26133661),
))
def test_convert_mag(from_value, to_unit, to_value):
    """These tests are not necessarily for accuracy.  Accuracy tests are
    reserved for the fluxd_to_mag and mag_to_fluxd tests.  Except, AB
    to ST conversion:

        bp = synphot.SpectralElement.from_file(
            'sbpy/photometry/tests/data/sdss-r.fits')
        m = -2.5 * np.log10((3630.78 * u.Jy).to(
            'erg/(s cm2 AA)', u.spectral_density(bp.pivot())).value
            / 3.63078e-9)

    Out[8]: 0.2613366094776436

    """

    fn = os.path.join(os.path.dirname(__file__), 'data', 'sdss-r.fits')
    bp = synphot.SpectralElement.from_file(fn)
    to = cal.convert_mag(from_value, to_unit, bandpass=bp)
    assert np.isclose(to.value, to_value)
