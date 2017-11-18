# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..sun import solar_fluxd

def test_solar_fluxd_singleton_K93():
    fluxd = solar_fluxd(1 * u.um, 'W/(m2 um)', source='K93')
    assert np.isclose(fluxd.value, 730.629964229039)

def test_solar_fluxd_singleton_C96():
    fluxd = solar_fluxd(1 * u.um, 'W/(m2 um)', source='C96')
    assert np.isclose(fluxd.value, 724.7101409359462)

def test_solar_fluxd_binning_K93():
    wave = np.logspace(-0.5, -0.7, 10) * u.um
    fluxd = solar_fluxd(wave, 'W/(m2 um)', source='K93')
    assert np.allclose(fluxd.value, [ 733.33690989,  633.60332033,  347.81298649,  273.89028631,        113.28964751,   67.61440112,   52.98876192,   54.78504952,         32.16428268,    5.05709748])

def test_solar_fluxd_binning_C96():
    wave = np.logspace(-0.5, -0.7, 10) * u.um
    fluxd = solar_fluxd(wave, 'W/(m2 um)', source='C96')
    assert np.allclose(fluxd.value, [ 738.32457762,  645.39256746,  355.15314878,  280.57101803,        109.40310067,   49.89236658,   47.65261149,   51.73931787,         25.5313468 ,    6.17967145])

