# Licensed under a 3-clause BSD style license - see LICENSE.rst

from sbpy.data.core import FieldError
import pytest
import numpy as np
import astropy.units as u
from astropy.table import Table, QTable
from astropy.time import Time
from ... import data as sbd
from ..decorators import quantity_to_dataclass, dataclass_input


def test_quantity_to_dataclass_single():
    @quantity_to_dataclass(eph=('rh', sbd.Ephem))
    def temperature(eph, a, b, c, d=5):
        return 278 * u.K / np.sqrt(eph['rh'] / u.au)

    rh = 1 * u.au
    eph = sbd.Ephem.from_dict({'rh': 1 * u.au})
    assert temperature(rh, 1, 2, 3) == temperature(eph, 7, 8, 9, d=10)
    assert np.isclose(temperature(rh, 4, 5, 6).value, 278.0)
    with pytest.raises(FieldError):
        temperature(1, 1, 2, 3)
    with pytest.raises(FieldError):
        temperature(1 * u.s, 1, 2, 3)


def test_quantity_to_dataclass_multiple():
    @quantity_to_dataclass(eph=('rh', sbd.Ephem), orbit=('a', sbd.Orbit))
    def contrived(eph, orbit):
        return (eph['rh'] / orbit['a']).decompose()

    rh = 1 * u.au
    a = 2 * u.au
    assert np.isclose(contrived(rh, a).value, 0.5)
    with pytest.raises(FieldError):
        contrived(1, a)
    with pytest.raises(FieldError):
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
    with pytest.raises(FieldError):
        contrived(1, a, R)
    with pytest.raises(FieldError):
        contrived(rh, 2 * u.s, R)


u.imperial.enable()


@pytest.mark.parametrize('field, arg, test_invalid_unit', (
    ('a', 1 * u.au, 's'),
    ('e', 0.5, 's'),
    ('inc', 10 * u.deg, 's'),
    ('delta_v', 0.5 * u.km/u.s, 'm'),
    ('P', 9. * u.hour, 'm'),
    ('ra_rate', 10 * u.rad / u.day, 's'),
    ('V', 3 * u.mag, 'm'),
    ('lun_illum', 10 * u.percent, 'm'),
    ('area_3sigma', 10 * u.deg**2, 'm'),
    ('area_3sigma', 10 * u.sr, 'm'),
    ('sband_3sigma', 5 * u.GHz, 's'),
    ('temp', 300 * u.K, 'm'),  # looks like astropy has problem handling
    # temperature units
    # ('lgint300', 300 * u.Unit('W/(m**2 sr)'), 'm/s'),
    ('eup_j', 1 * u.J, 'm'),
    ('eup_j', 1 * u.eV, 'm')
))
def test_quantity_to_dataclass_dimensions(field, arg, test_invalid_unit):
    @quantity_to_dataclass(x=(field, sbd.DataClass))
    def test(x):
        return x[field]

    assert all(test(arg) == arg)
    with pytest.raises(FieldError):
        test(1 * u.Unit(test_invalid_unit))


@pytest.mark.parametrize('arg, test_arg', (
    (Time('2019-01-01'), Time('2019-01-01')),
    (Time(['2019-01-01', '2019-01-02']), Time(['2019-01-01', '2019-01-02'])),
    ([Time('2019-01-01'), Time('2019-01-02')],
     Time(['2019-01-01', '2019-01-02']))
))
def test_quantity_to_dataclass_timetype(arg, test_arg):
    @quantity_to_dataclass(x=('epoch', sbd.Obs))
    def test(x):
        return x['epoch']

    assert np.all(test(arg) == test_arg)


@pytest.mark.parametrize('arg', (
    [1],
    ['2015-01-01'],
    [1 * u.s],
    [1, 2],
    [1, 2] * u.s
))
def test_quantity_to_dataclass_timetype_error(arg):
    @quantity_to_dataclass(x=('epoch', sbd.Obs))
    def test(x):
        return x['epoch']

    with pytest.raises(FieldError):
        test(arg)


def test_quantity_to_dataclass_equivalencies():
    @quantity_to_dataclass(
        x=('sband_3sigma', sbd.Obs), equivalencies=u.spectral()
    )
    def test(x):
        return x['sband_3sigma']

    assert test(1 * u.mm) == (1 * u.mm)
    assert test(300 * u.GHz) == (300 * u.GHz)


def test_quantity_to_dataclass_dataclasserror():
    @quantity_to_dataclass(x=('test', sbd.Obs))
    def test(x):
        return x['test']

    with pytest.raises(sbd.DataClassError):
        test(1)


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
    def temperature(eph: sbd.Ephem = None):
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
