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


def test_synphot_import_fail():
    with mock.patch.dict(sys.modules, {'synphot': None}):
        assert spectral_density_vega(1 * u.um) == []


def test_enable():
    with enable():
        assert str(u.Unit('VEGA')) == 'VEGA'
        assert str(u.Unit('VEGAflux')) == 'VEGA'
        assert str(u.Unit('mag(VEGA)')) == 'mag(VEGA)'
        assert str(u.Unit('mag(VEGAflux)')) == 'mag(VEGA)'
        assert str(u.Unit('JM')) == 'JM'
        assert str(u.Unit('JMflux')) == 'JM'
        assert str(u.Unit('mag(JM)')) == 'mag(JM)'
        assert str(u.Unit('mag(JMflux)')) == 'mag(JM)'


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

    Compare to Willmer 2018 Vega-mag zerpoints (Vega assumed to be
    0.03 mag).

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
