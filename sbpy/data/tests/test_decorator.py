# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from .. import quantity_to_dataclass, Ephem, Orbit, Phys


def test_quantity_to_dataclass_single():
    @quantity_to_dataclass('eph', 'rh', Ephem)
    def temperature(eph):
        return 278 * u.K / np.sqrt(eph['rh'] / u.au)

    rh = 1 * u.au
    eph = Ephem.from_dict({'rh': 1 * u.au})
    assert temperature(rh) == temperature(eph)
    assert np.isclose(temperature(rh).value, 278.0)


def test_quantity_to_dataclass_double():
    @quantity_to_dataclass('eph', 'rh', Ephem)
    @quantity_to_dataclass('orbit', 'a', Orbit)
    def contrived(eph, orbit):
        return (eph['rh'] / orbit['a']).decompose()

    rh = 1 * u.au
    a = 2 * u.au
    assert np.isclose(contrived(rh, a).value, 0.5)


def test_quantity_to_dataclass_triple():
    @quantity_to_dataclass('eph', 'rh', Ephem)
    @quantity_to_dataclass('orbit', 'a', Orbit)
    @quantity_to_dataclass('phys', 'R', Phys)
    def contrived(eph, orbit, phys):
        return (eph['rh'] / orbit['a']).decompose() * phys['R']

    rh = 1 * u.au
    a = 2 * u.au
    R = 100 * u.km
    assert np.isclose(contrived(rh, a, R).value, 50)
