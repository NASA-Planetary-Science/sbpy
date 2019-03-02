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
from .. import core


def test_synphot_import_fail():
    # importlib rigamarole is supposedly not needed, but this
    # didn't work without it. MSK / 2019 Mar 02 / Python 3.6.7
    with mock.patch.dict(sys.modules, {'synphot': None}):
        importlib.reload(core)
        assert core.synphot is None

        with pytest.raises(ImportError):
            core.spectral_density_vega(1 * u.um)

    importlib.reload(core)


@pytest.mark.parametrize('wf, fluxd, to', (
    (5557.5 * u.AA, 3.44e-9 * u.Unit('erg/(cm2 s AA)'), 0.03 * core.VEGAMAG),
    (5557.5 * u.AA, 0.03 * core.VEGAMAG, 3.44e-9 * u.Unit('erg/(cm2 s AA)')),
    (5557.5 * u.AA, 3.544e-23 * u.Unit('W/(m2 Hz)'), 0.03 * core.VEGAMAG),
    (5557.5 * u.AA, 0.03 * core.VEGAMAG, 3.544e-23 * u.Unit('W/(m2 Hz)')),
    (539.44 * u.THz, 3.544e-23 * u.Unit('W/(m2 Hz)'), 0.03 * core.VEGAMAG),
    (539.44 * u.THz, 0.03 * core.VEGAMAG, 3.544e-23 * u.Unit('W/(m2 Hz)')),
))
def test_spectral_density_vega_wf(wf, fluxd, to):
    """Test VEGAMAG conversions at wavelength / frequency.

    Flux density at 5557.5 AA is from Bohlin 2014 (0.5% uncertainty).

    """
    v = fluxd.to(to.unit, core.spectral_density_vega(wf))
    assert v.unit == to.unit
    if to.unit == core.VEGAMAG:
        assert np.isclose(v.value, to.value, atol=0.005)
    else:
        assert np.isclose(v.value, to.value, rtol=0.005)


@pytest.mark.parametrize('filename, fluxd, to, tol', (
    ('cousins_i_004_syn.fits', 2478.76 * u.Unit('Jy'), 0 * core.VEGAMAG,
     0.014),
    ('cousins_i_004_syn.fits', 0 * core.VEGAMAG, 2478.76 * u.Unit('Jy'),
     0.013),
    ('sdss-r.fits', 2.55856e-9 * u.Unit('erg/(s cm2 AA)'), 0 * core.VEGAMAG,
     0.005),
    ('sdss-r.fits', 0 * core.VEGAMAG, 2.55856e-9 * u.Unit('erg/(s cm2 AA)'),
     0.005),
    ('wfc3_uvis_f438w_004_syn.fits',
     4278.69 * u.Unit('Jy'), 0 * core.VEGAMAG, 0.005),
    ('wfc3_uvis_f438w_004_syn.fits',
     0 * core.VEGAMAG, 4278.69 * u.Unit('Jy'), 0.005),
    ('wfc3_uvis_f606w_004_syn.fits',
     2.97294e-9 * u.Unit('erg/(s cm2 AA)'), 0 * core.VEGAMAG, 0.012),
    ('wfc3_uvis_f606w_004_syn.fits',
     0 * core.VEGAMAG, 2.97294e-9 * u.Unit('erg/(s cm2 AA)'), 0.005),
))
def test_spectral_density_vega_bp(filename, fluxd, to, tol):
    """Test VEGAMAG conversions for bandpasses.

    Compare to Willmer 2018 Vega-mag zerpoints.

    """
    fn = get_pkg_data_filename(os.path.join(
        '..', '..', 'photometry', 'data', filename))
    bp = synphot.SpectralElement.from_file(fn)

    v = fluxd.to(to.unit, core.spectral_density_vega(bp))
    assert v.unit == to.unit
    if to.unit == core.VEGAMAG:
        assert np.isclose(v.value, to.value, atol=tol)
    else:
        assert np.isclose(v.value, to.value, rtol=tol)
