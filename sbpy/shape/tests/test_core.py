# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy import units as u

from ...calib import Sun, solar_spectrum
from ...data.ephem import Ephem
from ..core import Shape


def test_Shape_incident_sunlight():
    # Test branching on wfb input, scaling by rh
    wave = np.linspace(0.4, 0.9) * u.um
    unit = u.Jy
    eph = Ephem.from_dict({"rh": 2 * u.au})
    with solar_spectrum.set(Sun.from_array(wave, 1000 * unit)):
        S = Shape._incident_sunlight(wave, eph, unit, True)
        assert u.allclose(S, 250 * u.Jy)

        S = Shape._incident_sunlight(550 * u.nm, eph, unit, False)
        assert u.isclose(S, 250 * u.Jy, atol=0.02 * unit)
