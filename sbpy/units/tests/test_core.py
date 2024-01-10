# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import pytest
import numpy as np
import astropy.units as u
from astropy.io import ascii
from astropy.utils.data import get_pkg_data_path

from .. import core
from ..core import *
from ...photometry import bandpass
from ...calib import (vega_spectrum, vega_fluxd, solar_fluxd,
                      solar_spectrum, Sun, Vega)
from ...exceptions import RequiredPackageUnavailable


@pytest.mark.parametrize('unit,test', (
    ('VEGA', 'VEGA'),
    ('VEGAflux', 'VEGA'),
    ('mag(VEGA)', 'mag(VEGA)'),
    ('mag(VEGAflux)', 'mag(VEGA)'),
    ('JM', 'JM'),
    ('JMflux', 'JM'),
    ('mag(JM)', 'mag(JM)'),
    ('mag(JMflux)', 'mag(JM)')))
def test_enable(unit, test):
    with core.enable():
        assert str(u.Unit(unit)) == test


def test_hundred_nm():
    assert (1 * hundred_nm).to(u.nm).value == 100


@pytest.mark.parametrize('wf, fluxd, to', (
    (5557.5 * u.AA, 3.44e-9 * u.Unit('erg/(cm2 s AA)'), 0 * VEGAmag),
    (5557.5 * u.AA, 3.44e-9 * u.Unit('erg/(cm2 s AA)'), 0.03 * JMmag),
    (5557.5 * u.AA, 0 * VEGAmag, 3.44e-9 * u.Unit('erg/(cm2 s AA)')),
    (5557.5 * u.AA, 0.03 * JMmag, 3.44e-9 * u.Unit('erg/(cm2 s AA)')),
    (5557.5 * u.AA, 3.544e-23 * u.Unit('W/(m2 Hz)'), 0 * VEGAmag),
    (5557.5 * u.AA, 3.544e-23 * u.Unit('W/(m2 Hz)'), 0.03 * JMmag),
    (5557.5 * u.AA, 0 * VEGAmag, 3.544e-23 * u.Unit('W/(m2 Hz)')),
    (5557.5 * u.AA, 0.03 * JMmag, 3.544e-23 * u.Unit('W/(m2 Hz)')),
    (539.44 * u.THz, 3.544e-23 * u.Unit('W/(m2 Hz)'), 0 * VEGAmag),
    (539.44 * u.THz, 3.544e-23 * u.Unit('W/(m2 Hz)'), 0.03 * JMmag),
    (539.44 * u.THz, 0 * VEGAmag, 3.544e-23 * u.Unit('W/(m2 Hz)')),
    (539.44 * u.THz, 0.03 * JMmag, 3.544e-23 * u.Unit('W/(m2 Hz)')),
))
def test_spectral_density_vega_wf(wf, fluxd, to):
    """Test vega magnitude system conversions for wavelength / frequency.

    Flux density at 5557.5 AA is from Bohlin 2014 (0.5% uncertainty).

    """

    pytest.importorskip("synphot")

    v = fluxd.to(to.unit, spectral_density_vega(wf))
    assert v.unit == to.unit
    if to.unit in (VEGAmag, JMmag):
        assert np.isclose(v.value, to.value, atol=0.001)
    else:
        assert np.isclose(v.value, to.value, rtol=0.001)


@pytest.mark.parametrize('filename, fluxd, to, tol', (
    ('sdss-r.fits', 2.55856e-9 * u.Unit('erg/(s cm2 AA)'), 0 * JMmag,
     0.005),
    ('sdss-r.fits', 0 * JMmag, 2.55856e-9 * u.Unit('erg/(s cm2 AA)'),
     0.005),
    ('wfc3_uvis_f438w_004_syn.fits',
     4278.69 * u.Unit('Jy'), 0 * JMmag, 0.005),
    ('wfc3_uvis_f438w_004_syn.fits',
     0 * JMmag, 4278.69 * u.Unit('Jy'), 0.005),
    ('wfc3_uvis_f606w_004_syn.fits',
     2.97294e-9 * u.Unit('erg/(s cm2 AA)'), 0 * JMmag, 0.012),
    ('wfc3_uvis_f606w_004_syn.fits',
     0 * JMmag, 2.97294e-9 * u.Unit('erg/(s cm2 AA)'), 0.005),
))
def test_spectral_density_vega_bp(filename, fluxd, to, tol):
    """Test VEGAmag conversions for bandpasses.

    Compare to Willmer 2018 Vega-mag zeropoints.  According to Eq. 13,
    Table 3 assumes Vega is 0 mag, but only the Cousins I filter
    tested here agrees with that definition.  The rest have better
    agreement with 0.03 mag.

    """

    synphot = pytest.importorskip("synphot")

    fn = get_pkg_data_path(os.path.join(
        '..', '..', 'photometry', 'data', filename))
    bp = synphot.SpectralElement.from_file(fn)

    v = fluxd.to(to.unit, spectral_density_vega(bp))
    assert v.unit == to.unit
    if to.unit in (VEGAmag, JMmag):
        assert np.isclose(v.value, to.value, atol=tol)
    else:
        assert np.isclose(v.value, to.value, rtol=tol)


def test_spectral_density_vega_undefinedsourceerror():
    pytest.importorskip("synphot")
    with vega_spectrum.set(Vega(None)):
        assert spectral_density_vega([1, 2, 3] * u.um) == []


@pytest.mark.parametrize('fluxd, wfb, f_sun, ref', (
    (3.4 * VEGAmag, "Johnson V", None, 0.02865984),
    (3.4 * VEGAmag, 5500 * u.AA, None, 0.02774623),
    (1.56644783e-09 * u.Unit('W/(m2 um)'), 'V',
        1839.93273227 * u.Unit('W/(m2 um)'), 0.02865984),
    (1.55728147e-24 * u.Unit('W/(m2 Hz)'), 'V',
        1.86599755e-12 * u.Unit('W/(m2 Hz)'), 0.02809415),
    (3.4 * VEGAmag, 'V', -26.77471503 * VEGAmag, 0.02865984),
    (3.4 * u.ABmag, 'V', -26.77471503 * u.ABmag, 0.02865984),
    (3.4 * u.mag, 'V', -26.77471503 * u.mag, 0.02865984)
))
def test_reflectance_ref(fluxd, wfb, f_sun, ref):
    """Test conversion from flux to reflectance

    Use Ceres as the reference: H = 3.4 mag, radius = 460 km, average
    bidirectional reflectance at zero phase angle 0.029 (~1/pi of geometric
    albedo).

    """

    pytest.importorskip("synphot")

    try:
        # passes for "Johnson V", but fails for "V" band
        # we will set "V" band to a specific value below
        wfb = bandpass(wfb) if type(wfb) is str else wfb
    except KeyError:
        pass

    xsec = 6.648e5 * u.km**2

    with vega_fluxd.set({'V': u.Quantity(3.589e-9, 'erg/(s cm2 AA)')}):
        with solar_fluxd.set({wfb: f_sun}):
            r = fluxd.to('1/sr', reflectance(wfb, cross_section=xsec))
    assert r.unit == u.sr**-1
    assert np.isclose(r.value, ref)


@pytest.mark.parametrize('fluxd, wfb, f_sun, radius', (
    (3.4 * VEGAmag, "Johnson V", None, 460.01351274),
    (3.4 * VEGAmag, 5500 * u.AA, None, 452.62198065),
    (1.56644783e-09 * u.Unit('W/(m2 um)'), 'V',
        1839.93273227 * u.Unit('W/(m2 um)'), 460.01351274),
    (1.55728147e-24 * u.Unit('W/(m2 Hz)'), 'V',
        1.86599755e-12 * u.Unit('W/(m2 Hz)'), 455.45095634),
    (3.4 * VEGAmag, 'V', -26.77471503 * VEGAmag, 460.01351274),
    (3.4 * u.ABmag, 'V', -26.77471503 * u.ABmag, 460.01351274),
    (3.4 * u.mag, 'V', -26.77471503 * u.mag, 460.01351274)
))
def test_reflectance_xsec(fluxd, wfb, f_sun, radius):
    """Test conversion from flux to reflectance

    Use Ceres as the reference: H = 3.4 mag, radius = 460 km, average
    bidirectional reflectance at zero phase angle 0.029 (~1/pi of geometric
    albedo 0.09).

    """

    pytest.importorskip("synphot")

    try:
        # passes for "Johnson V", but fails for "V" band
        # we will set "V" band to a specific value below
        wfb = bandpass(wfb) if type(wfb) is str else wfb
    except KeyError:
        pass

    ref = 0.02865984 / u.sr
    with vega_fluxd.set({'V': u.Quantity(3.589e-9, 'erg/(s cm2 AA)')}):
        with solar_fluxd.set({wfb: f_sun}):
            xs = fluxd.to('km2', reflectance(wfb, reflectance=ref))
            ra = np.sqrt(xs / np.pi)
    assert ra.unit == u.km
    assert np.isclose(ra.value, radius)


def test_reflectance_exception():
    assert reflectance('B', reflectance=1/u.sr) == []


def test_reflectance_spec():
    """Test conversion from flux spectrum to reflectance spectrum

    Use Tempel 1 spectrum collected by Deep Impact HRI-IR at
    2005-07-04 05:32 UT as the standard.  Data file
    'data/hi05070405_9000036-avg-spec.txt', data from McLaughlin et al (2014).
    """

    pytest.importorskip("synphot")

    ifov = 1e-5 * u.rad
    delta = 15828 * u.km
    rh = 1.5 * u.au

    # Tempel 1 spectrum, includes reference solar spectrum
    fn = get_pkg_data_path(
        os.path.join('data', 'hi05070405_9000036-avg-spec.txt'))
    t1 = ascii.read(fn)
    sun = Sun.from_array(t1['wave'] * u.um,
                         t1['flux_sun_nu'] * u.Unit('W/(m2 Hz)'))
    with solar_spectrum.set(sun):
        wave = t1['wave'] * u.um
        spec = ((t1['spec'] * u.Unit('W/(m2 um sr)') * ifov**2)
                * (delta * rh / u.au**2)**2).to('W/(m2 um)')
        spec_nu = spec.to('W/(m2 Hz)', u.spectral_density(wave))

        xsec = (ifov * delta)**2 / u.sr
        ref1 = spec.to('1/sr', reflectance(
            wave, cross_section=xsec, interpolate=True))
        ref2 = spec_nu.to('1/sr', reflectance(
            wave, cross_section=xsec, interpolate=True))

    # use built-in solar spectrum
    ref3 = spec.to('1/sr', reflectance(wave, cross_section=xsec))

    for ref in (ref1, ref2, ref3):
        assert ref.unit == '1/sr'
        assert np.allclose(ref.value, t1['ref'])


@pytest.mark.parametrize('value, delta, test', (
    (1 * u.arcsec, 1 * u.au, np.tan(1 * u.arcsec) * u.au),
    (725.27 * u.km, 1 * u.au, np.tan(1 * u.arcsec) * u.rad),
))
def test_projected_size(value, delta, test):
    test = test.decompose()
    result = value.to(test.unit, projected_size(delta))
    assert np.isclose(result.value, test.value)
