# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..dust import *

def test_phase_HalleyMarcus():
    assert np.isclose(phase_HalleyMarcus(0 * u.deg), 1.0)
    assert np.isclose(phase_HalleyMarcus(15 * u.deg), 5.8720e-01)
    assert np.isclose(phase_HalleyMarcus(14.5 * u.deg), 0.5959274462322928)

class TestAfrho:
    def test_init(self):
        afrho = Afrho(1000 * u.cm)
        assert afrho.value == 1000
        assert afrho.unit == u.cm

    def test_scaler_ops(self):
        afrho = Afrho(1000 * u.cm)
        afrho = afrho / 2
        assert afrho == 500 * u.cm

    def test_quantity_ops(self):
        afrho = Afrho(1000 * u.cm)
        v = afrho * 2 * u.cm
        assert v == 2000 * u.cm**2
        assert not isinstance(v, Afrho)

    def test_from_fluxd(self):
        fluxd = 6.730018324465526e-14 * u.W / u.m**2 / u.um
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        S = 1869 * u.W / u.m**2 / u.um
        afrho = Afrho.from_fluxd(None, fluxd, aper, eph=eph, S=S)
        assert np.isclose(afrho.cm, 1000)

    def test_flux(self):
        afrho = Afrho(1000, 'cm')
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        S = 1869 * u.W / u.m**2 / u.um
        fluxd = afrho.fluxd(None, aper, eph=eph, S=S)
        assert fluxd.unit.is_equivalent(S.unit)
        assert np.isclose(fluxd.value, 6.730018324465526e-14)

    def test_to_phase(self):
        afrho = Afrho(10 * u.cm).to_phase(15 * u.deg, 0 * u.deg)
        assert np.isclose(afrho.cm, 5.8720)
    
