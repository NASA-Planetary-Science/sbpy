# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from astropy.table import Table, QTable
import pytest
from ... import data as sbd
from ..decorators import *


def test_quantity_to_dataclass_single():
    @quantity_to_dataclass(eph=('rh', sbd.Ephem))
    def temperature(eph, a, b, c, d=5):
        return 278 * u.K / np.sqrt(eph['rh'] / u.au)

    rh = 1 * u.au
    eph = sbd.Ephem.from_dict({'rh': 1 * u.au})
    assert temperature(rh, 1, 2, 3) == temperature(eph, 7, 8, 9, d=10)
    assert np.isclose(temperature(rh, 4, 5, 6).value, 278.0)
    with pytest.raises(TypeError):
        temperature(1, 1, 2, 3)
    with pytest.raises(u.UnitsError):
        temperature(1 * u.s, 1, 2, 3)


def test_quantity_to_dataclass_multiple():
    @quantity_to_dataclass(eph=('rh', sbd.Ephem), orbit=('a', sbd.Orbit))
    def contrived(eph, orbit):
        return (eph['rh'] / orbit['a']).decompose()

    rh = 1 * u.au
    a = 2 * u.au
    assert np.isclose(contrived(rh, a).value, 0.5)
    with pytest.raises(TypeError):
        contrived(1, a)
    with pytest.raises(u.UnitsError):
        contrived(rh, 2 * u.s)


def test_quantity_to_dataclass_stacked():
    """Note, this is not a preferred use case."""
    @quantity_to_dataclass(eph=('rh', sbd.Ephem))
    @quantity_to_dataclass(orbit=('a', sbd.Orbit))
    @quantity_to_dataclass(phys=('R', sbd.Phys))
    def contrived(eph, orbit, phys):
        return (eph['rh'] / orbit['a']).decompose() * phys['R']

    rh = 1 * u.au
    a = 2 * u.au
    R = 100 * u.km
    assert np.isclose(contrived(rh, a, R).value, 50)
    with pytest.raises(TypeError):
        contrived(1, a, R)
    with pytest.raises(u.UnitsError):
        contrived(rh, 2 * u.s, R)


def test_quantity_to_dataclass_optional():
    @quantity_to_dataclass(eph=('rh', sbd.Ephem))
    def temperature(eph=None):
        if eph is None:
            rh = 1 * u.au
        else:
            rh = eph['rh']
        return 278 / np.sqrt(rh / u.au)

    assert temperature(1 * u.au) == temperature()
    assert np.isclose(temperature(2 * u.au).value, 278 / np.sqrt(2))


@pytest.mark.parametrize('eph, orb', (
    ({'rh': 1 * u.au},
     Table([[2] * u.au], names=['a'])),

    (sbd.Obs.from_dict({'rh': 1 * u.au, 'fluxd': 1 * u.Jy}),
     QTable([[2] * u.au], names=['a']))
))
def test_dataclass_input_kwargs(eph, orb):
    @dataclass_input(eph=sbd.Ephem, orb=sbd.Orbit)
    def func(eph, orb):
        return eph['rh'] / orb['a']

    assert np.isclose(func(eph, orb).value, 0.5)


def test_dataclass_input_annotation():
    @dataclass_input
    def func(eph: sbd.Ephem, orb: sbd.Orbit):
        return eph['rh'] / orb['a']

    eph = {'rh': 1 * u.au}
    orb = Table([[2] * u.au], names=['a'])
    assert np.isclose(func(eph, orb).value, 0.5)


def test_quantity_input_optional():
    @dataclass_input
    def temperature(eph: sbd.Ephem=None):
        if eph is None:
            rh = 1 * u.au
        else:
            rh = eph['rh']
        return 278 * u.K / np.sqrt(rh / u.au)

    eph = sbd.Ephem.from_dict({'rh': 1 * u.au})
    assert temperature(eph) == temperature()

    eph = sbd.Ephem.from_dict({'rh': 2 * u.au})
    assert np.isclose(temperature(eph).value, 278 / np.sqrt(2))


def test_dataclass_input_after_quantity_to_dataclass():
    @dataclass_input
    @quantity_to_dataclass(eph=('rh', sbd.Ephem))
    def func(eph: sbd.Ephem):
        return eph['rh']

    assert func(1 * u.au) == (1 * u.au)
