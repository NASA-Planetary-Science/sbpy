# Licensed under a 3-clause BSD style license - see LICENSE.rst

from numpy import pi
from astropy import units as u
from astropy.modeling.models import BlackBody

from ..neatm import neatm


def test_neatm():
    wave = 10 * u.um
    D = 1 * u.km
    A = 0.1
    epsilon = 0.95
    eta = 1.0
    eph = {"rh": 1 * u.au, "delta": 1 * u.au, "phase": 45 * u.deg}

    fnu, err = neatm(wave, D, A, epsilon, eta, eph, unit="Jy", epsrel=1e-5)
    assert err / fnu < 1e-5

    # assert u.isclose(fnu, 0.01776793 * u.Jy)
