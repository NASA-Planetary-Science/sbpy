# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import mock
import importlib
import pytest
import numpy as np
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
import synphot
from ..core import *


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
    with enable():
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

    Compare to Willmer 2018 Vega-mag zerpoints.  According to Eq. 13,
    Table 3 assumes Vega is 0 mag, but only the Cousins I filter
    tested here agrees with that definition.  The rest have better
    agreement with 0.03 mag.

    """
    fn = get_pkg_data_filename(os.path.join(
        '..', '..', 'photometry', 'data', filename))
    bp = synphot.SpectralElement.from_file(fn)

    v = fluxd.to(to.unit, spectral_density_vega(bp))
    assert v.unit == to.unit
    if to.unit in (VEGAmag, JMmag):
        assert np.isclose(v.value, to.value, atol=tol)
    else:
        assert np.isclose(v.value, to.value, rtol=tol)


def test_spectral_density_vega_synphot_import_fail():
    with mock.patch.dict(sys.modules, {'synphot': None}):
        assert spectral_density_vega(1 * u.um) == []


@pytest.mark.parametrize('mag, wfb, M_sun, ref', (
    (3.4 * VEGAmag, None, None, 0.02865984),
    (3.4 * VEGAmag, 'johnson_v', None, 0.02865984),
    (3.4 * VEGAmag, 'johnson_i', None, 0.01387486),
    (3.4 * VEGAmag, 5500 * u.AA, None, 0.02774623),
    (3.4 * VEGAmag, None, -25.0 * VEGAmag, 0.14694669),
    (3.4 * u.ABmag, None, -25.0 * u.ABmag, 0.14694669),
))
def test_magnitude_reflectance(mag, wfb, M_sun, ref):
    xsec = 6.648e5 * u.km**2
    r = mag.to('1/sr', reflectance(cross_section=xsec, wfb=wfb, M_sun=M_sun))
    assert r.unit == u.sr**-1
    assert np.isclose(r.value, ref)


@pytest.mark.parametrize('mag, wfb, M_sun, xsec', (
    (3.4 * VEGAmag, None, None, 664800.062),
    (3.4 * VEGAmag, 'johnson_v', None, 664800.062),
    (3.4 * VEGAmag, 'johnson_i', None, 321844.372),
    (3.4 * VEGAmag, 5500 * u.AA, None, 643607.586),
    (3.4 * VEGAmag, None, -25.0 * VEGAmag, 3408608.03),
    (3.4 * u.ABmag, None, -25.0 * u.ABmag, 3408608.03),
))
def test_magnitude_xsection(mag, wfb, M_sun, xsec):
    ref = 0.02865984 / u.sr
    xs = mag.to('km2', reflectance(reflectance=ref, wfb=wfb, M_sun=M_sun))
    assert xs.unit == u.km**2
    assert np.isclose(xs.value, xsec)
