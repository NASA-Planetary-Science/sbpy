# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from collections import OrderedDict
from sbpy.data import dimensions
import pytest
from copy import deepcopy
import astropy.units as u
from astropy.coordinates import Angle
from astropy.table import QTable, Column, Row
from astropy.time import Time
from ..core import DataClass, Conf, DataClassError, FieldError


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


def test_creation_single_row():
    """ test the creation of DataClass objects from dicts or arrays;
    single row only"""

    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1], [2], [3]], names=('aa', 'bb', 'cc'))
    ground_truth_2 = QTable([[1]*u.m, [2]*u.kg, [3]*u.cm/u.s],
                            names=('aa', 'bb', 'cc'))
    ground_truth_3 = QTable([[1]*u.m, [2], [3]*u.kg], names=('aa', 'bb', 'cc'))
    ground_truth_4 = QTable([[1], ['stuff'], [3]], names=('aa', 'bb', 'cc'))
    ground_truth_5 = QTable([[1]*u.km, [2], ['test']],
                            names=('aa', 'bb', 'cc'))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict(OrderedDict([('aa', 1), ('bb', 2),
                                                   ('cc', 3)]))
    test_dict_2 = DataClass.from_dict(OrderedDict([('aa', 1*u.m),
                                                   ('bb', 2*u.kg),
                                                   ('cc', 3*u.cm/u.s)]))
    test_dict_3 = DataClass.from_dict(OrderedDict([('aa', 1*u.m), ('bb', 2),
                                                   ('cc', 3*u.kg)]))
    test_dict_4 = DataClass.from_dict(OrderedDict([('aa', 1), ('bb', 'stuff'),
                                                   ('cc', 3)]))
    test_dict_5 = DataClass.from_dict(OrderedDict([('aa', 1*u.km), ('bb', 2),
                                                   ('cc', 'test')]))

    assert test_dict_1.table == ground_truth_1
    assert test_dict_2.table == ground_truth_2
    assert test_dict_3.table == ground_truth_3
    assert test_dict_4.table == ground_truth_4
    assert test_dict_5.table == ground_truth_5

    # test DataClass.from_rows for different cases
    test_array_1 = DataClass.from_rows([1, 2, 3],
                                       names=('aa', 'bb', 'cc'))
    test_array_2 = DataClass.from_rows([1*u.m, 2*u.kg, 3*u.cm/u.s],
                                       names=('aa', 'bb', 'cc'))
    test_array_3 = DataClass.from_rows([1*u.m, 2, 3*u.kg],
                                       names=('aa', 'bb', 'cc'))
    test_array_4 = DataClass.from_rows([1, 'stuff', 3],
                                       names=('aa', 'bb', 'cc'))
    test_array_5 = DataClass.from_rows([1*u.km, 2, 'test'],
                                       names=('aa', 'bb', 'cc'))

    assert test_array_1.table == ground_truth_1
    assert test_array_2.table == ground_truth_2
    assert test_array_3.table == ground_truth_3
    assert test_array_4.table == ground_truth_4
    assert test_array_5.table == ground_truth_5

    # test single row, single column
    ground_truth_6 = QTable([[1]], names=('aa',))

    test_dict_6 = DataClass.from_dict({'aa': 1})
    assert test_dict_6.table == ground_truth_6

    test_array_6 = DataClass.from_rows([1], names='aa')
    assert test_array_6.table == ground_truth_6

    # test units parameter
    test_array_2b = DataClass.from_rows([1, 2, 3],
                                        names=('aa', 'bb', 'cc'),
                                        units=(u.m, u.kg, u.cm/u.s))
    test_array_3b = DataClass.from_rows([1, 2, 3],
                                        names=('aa', 'bb', 'cc'),
                                        units=(u.m, None, u.kg))

    assert test_array_2b.table == ground_truth_2
    assert test_array_3b.table == ground_truth_3

    ground_truth_7 = QTable([[1]*u.kg], names=('aa',))
    test_array_7 = DataClass.from_rows([1], names='aa', units='kg')
    assert test_array_7.table == ground_truth_7

    with pytest.raises(DataClassError):
        DataClass.from_rows([1], names='aa', units=('m', 'kg'))

    # test single row starting with string
    ground_truth_8 = QTable(rows=[['a', 1]], names=('aa', 'bb'))
    test_array_8 = DataClass.from_rows(['a', 1], names=('aa', 'bb'))

    assert ground_truth_8 == test_array_8.table


def test_creation_multi_rows():
    """test the creation of DataClass objects from dicts or arrays;
    multiple rows"""
    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            names=('aa', 'bb', 'cc'))
    ground_truth_2 = QTable([[1, 2, 3]*u.m, [4, 5, 6]*u.kg,
                             [7, 8, 9]*u.cm/u.s], names=('aa', 'bb', 'cc'))
    ground_truth_3 = QTable([[1, 2, 3]*u.m, [4, 5, 6],
                             [7, 8, 9]*u.kg], names=('aa', 'bb', 'cc'))
    ground_truth_4 = QTable([[1, 2, 3], ['a', 'b', 'c'], [7, 8, 9]],
                            names=('aa', 'bb', 'cc'))
    ground_truth_5 = QTable([[1, 2, 3]*u.km, [4, 5, 6], ['a', 'b', 'c']],
                            names=('aa', 'bb', 'cc'))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict(OrderedDict([('aa', [1, 2, 3]),
                                                   ('bb', [4, 5, 6]),
                                                   ('cc', [7, 8, 9])]))
    test_dict_2 = DataClass.from_dict(OrderedDict(
        [('aa', [1, 2, 3]*u.m), ('bb', [4, 5, 6]*u.kg),
         ('cc', [7, 8, 9]*u.cm/u.s)]))
    test_dict_3 = DataClass.from_dict(OrderedDict([('aa', [1, 2, 3]*u.m),
                                                   ('bb', [4, 5, 6]),
                                                   ('cc', [7, 8, 9]*u.kg)]))
    test_dict_4 = DataClass.from_dict(OrderedDict([('aa', [1, 2, 3]),
                                                   ('bb', ['a', 'b', 'c']),
                                                   ('cc', [7, 8, 9])]))
    test_dict_5 = DataClass.from_dict(OrderedDict([('aa', [1, 2, 3]*u.km),
                                                   ('bb', [4, 5, 6]),
                                                   ('cc', ['a', 'b', 'c'])]))

    assert all(test_dict_1.table == ground_truth_1)
    assert all(test_dict_2.table == ground_truth_2)
    assert all(test_dict_3.table == ground_truth_3)
    assert all(test_dict_4.table == ground_truth_4)
    assert all(test_dict_5.table == ground_truth_5)

    # test DataClass.from_rows for different cases
    test_array_1 = DataClass.from_rows([[1, 4, 7], [2, 5, 8], [3, 6, 9]],
                                       names=('aa', 'bb', 'cc'))
    test_array_2 = DataClass.from_rows([[1*u.m, 4*u.kg, 7*u.cm/u.s],
                                        [2*u.m, 5*u.kg, 8*u.cm/u.s],
                                        [3*u.m, 6*u.kg, 9*u.cm/u.s]],
                                       names=('aa', 'bb', 'cc'))
    test_array_3 = DataClass.from_rows([[1*u.m, 4, 7*u.kg],
                                        [2*u.m, 5, 8*u.kg],
                                        [3*u.m, 6, 9*u.kg]],
                                       names=('aa', 'bb', 'cc'))
    test_array_4 = DataClass.from_rows([[1, 'a', 7],
                                        [2, 'b', 8],
                                        [3, 'c', 9]],
                                       names=('aa', 'bb', 'cc'))
    test_array_5 = DataClass.from_rows([[1*u.km, 4, 'a'],
                                        [2*u.km, 5, 'b'],
                                        [3*u.km, 6, 'c']],
                                       names=('aa', 'bb', 'cc'))

    assert all(test_array_1.table == ground_truth_1)
    assert all(test_array_2.table == ground_truth_2)
    assert all(test_array_3.table == ground_truth_3)
    assert all(test_array_4.table == ground_truth_4)
    assert all(test_array_5.table == ground_truth_5)


def test_creation_single_column():
    """test the creation of DataClass objects from dicts or arrays;
    single column only"""

    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1, 2, 3]], names=('aa',))
    ground_truth_2 = QTable([[1, 2, 3]*u.kg], names=('aa',))
    ground_truth_3 = QTable([['a', 'b', 'c']], names=('aa',))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict({'aa': [1, 2, 3]})
    test_dict_2 = DataClass.from_dict({'aa': [1, 2, 3]*u.kg})
    test_dict_3 = DataClass.from_dict({'aa': ['a', 'b', 'c']})

    assert all(test_dict_1.table == ground_truth_1)
    assert all(test_dict_2.table == ground_truth_2)
    assert all(test_dict_3.table == ground_truth_3)

    # test DataClass.from_columns for different cases
    test_array_1 = DataClass.from_columns([1, 2, 3], names='aa')
    test_array_2 = DataClass.from_columns([1, 2, 3]*u.kg, names='aa')
    test_array_3 = DataClass.from_columns(['a', 'b', 'c'], names='aa')

    assert all(test_array_1.table == ground_truth_1)
    assert all(test_array_2.table == ground_truth_2)
    assert all(test_array_3.table == ground_truth_3)

    # test single row, single column
    ground_truth_4 = QTable([[1]], names=('aa',))

    test_dict_4 = DataClass.from_dict({'aa': 1})
    assert test_dict_4.table == ground_truth_4

    test_array_4 = DataClass.from_columns([1], names='aa')
    assert test_array_4.table == ground_truth_4


def test_creation_multi_column():
    """test the creation of DataClass objects from dicts or arrays;
    multiple columns"""

    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            names=('aa', 'bb', 'cc'))
    ground_truth_2 = QTable([[1, 2, 3]*u.kg, [4, 5, 6]*u.m/u.s],
                            names=('aa', 'bb'))
    ground_truth_3 = QTable([[1, 2, 3], [4, 5, 6]*u.m/u.s],
                            names=('aa', 'bb'))
    ground_truth_4 = QTable([[1, 2, 3], ['a', 'b', 'c']],
                            names=('aa', 'bb'))
    ground_truth_5 = QTable([[1, 2, 3], [4, 5, 6]*u.s/u.kg, ['a', 'b', 'c']],
                            names=('aa', 'bb', 'cc'))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict({'aa': [1, 2, 3], 'bb': [4, 5, 6],
                                       'cc': [7, 8, 9]})
    test_dict_2 = DataClass.from_dict({'aa': [1, 2, 3]*u.kg,
                                       'bb': [4, 5, 6]*u.m/u.s})
    test_dict_3 = DataClass.from_dict({'aa': [1, 2, 3],
                                       'bb': [4, 5, 6]*u.m/u.s})
    test_dict_4 = DataClass.from_dict({'aa': [1, 2, 3],
                                       'bb': ['a', 'b', 'c']})
    test_dict_5 = DataClass.from_dict({'aa': [1, 2, 3],
                                       'bb': [4, 5, 6]*u.s/u.kg,
                                       'cc': ['a', 'b', 'c']})
    assert all(test_dict_1.table == ground_truth_1)
    assert all(test_dict_2.table == ground_truth_2)
    assert all(test_dict_3.table == ground_truth_3)
    assert all(test_dict_4.table == ground_truth_4)
    assert all(test_dict_5.table == ground_truth_5)

    # test DataClass.from_columns for different cases
    test_array_1 = DataClass.from_columns(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]], names=('aa', 'bb', 'cc'))
    test_array_2 = DataClass.from_columns([[1, 2, 3]*u.kg, [4, 5, 6]*u.m/u.s],
                                          names=('aa', 'bb'))
    test_array_3 = DataClass.from_columns([[1, 2, 3], [4, 5, 6]*u.m/u.s],
                                          names=('aa', 'bb'))
    test_array_4 = DataClass.from_columns([[1, 2, 3], ['a', 'b', 'c']],
                                          names=('aa', 'bb'))
    test_array_5 = DataClass.from_columns([[1, 2, 3], [4, 5, 6]*u.s/u.kg,
                                           ['a', 'b', 'c']],
                                          names=('aa', 'bb', 'cc'))
    assert all(test_array_1.table == ground_truth_1)
    assert all(test_array_2.table == ground_truth_2)
    assert all(test_array_3.table == ground_truth_3)
    assert all(test_array_4.table == ground_truth_4)
    assert all(test_array_5.table == ground_truth_5)


def test_get_set():
    """ test the get and set methods"""

    data = DataClass.from_dict(
        OrderedDict((('aa', [1, 2, 3]),
                     ('bb', [4, 5, 6]*u.m),
                     ('cc', [7, 8, 8]))))

    # get a single column
    x = data['aa']
    assert len(x) == 3
    assert isinstance(x, Column)
    x = data['bb']
    assert len(x) == 3
    assert isinstance(x, u.Quantity)

    # get a list of columns
    x = data[['aa', 'cc']]
    assert len(x.field_names) == 2
    assert isinstance(x, DataClass)

    # mask rows
    masked = data[[True, False, False]]
    assert len(masked) == 1
    assert masked['bb'][0] == 4*u.m
    assert isinstance(masked, DataClass)

    # get single row
    shortened = data[1]
    assert shortened['aa'] == 2
    assert isinstance(shortened, DataClass)

    # get list of rows
    shortened = data[[0, 1]]
    assert len(shortened) == 2
    assert isinstance(shortened, DataClass)

    # get slice
    assert len(shortened) == 2
    assert isinstance(shortened, DataClass)

    # modify an existing column
    data['aa'][:] = [0, 0, 0]
    assert data['aa'][0] == 0

    with pytest.raises(KeyError):
        data['dd']

    # add non-existing column using set
    data['z'] = 3
    assert len(data['z'] == 3)
    assert isinstance(data, DataClass)

    # modify existing column using set
    data['z'] = 2
    assert data['z'][1] == 2
    assert isinstance(data, DataClass)

    # modify existing column using alternative name
    data = DataClass.from_dict(
        {'rh': [1, 2, 3] * u.au,
         'delta': [4, 5, 6] * u.au})
    data['r'] = [7, 8, 9] * u.au
    assert u.allclose(data['rh'], [7, 8, 9] * u.au)
    assert set(data.field_names) == {'rh', 'delta'}


def test_units():
    """ test units on multi-row tables """

    ground_truth = QTable([[1, 2, 3]*u.Unit('m'),
                           [4, 5, 6]*u.m/u.s,
                           ['a', 'b', 'c']],
                          names=('aa', 'bb', 'cc'))

    assert ((ground_truth['aa']**2).unit == 'm2')

    test_dict = DataClass.from_dict(
        OrderedDict((('aa', [1, 2, 3]*u.m),
                     ('bb', [4, 5, 6]*u.m/u.s),
                     ('cc', ['a', 'b', 'c']))))
    assert all(test_dict.table == ground_truth)

    test_array = DataClass.from_columns([[1, 2, 3]*u.m,
                                         [4, 5, 6]*u.m/u.s,
                                         ['a', 'b', 'c']],
                                        names=('aa', 'bb', 'cc'))
    assert all(test_array.table == ground_truth)


def test_verify_fields():
    data = DataClass.from_dict({
        'name': 'asdf',
        'RA': Angle(1 * u.deg),
        'dRA': 1 * u.deg / u.day,
        'eup_J': 1 * u.J,
        'sband_3sigma': 1 * u.Hz,
        # 'lgint': 1 * u.Hz / u.m**2,
        'col_density': 1 * u.m**-2,
        'au': 1 * u.s**-1,
        'rh': 1 * u.m,
        'V': 1 * u.mag,
        'surfbright': 1 * u.mag / u.sr,
        'frac_illum': 1 * u.percent,
        'area_3sigma': 1 * u.steradian,
        'temperature': 273 * u.K,
        'period': 1 * u.s,
        'beta': 1 * u.s * u.m**2,
        'delta-v': 1 * u.m / u.s,
        'epoch': Time('2021-09-25')
    })
    # explicitly call for verification
    data.verify_fields()


@pytest.mark.parametrize(
    'field,quantity',
    (
        ['RA', 1 * u.m],
        ['dRA', 1 * u.m / u.day],
        ['eup_J', 1 * u.kg],
        ['sband_3sigma', 1 * u.s],
        # ['lgint', 1 * u.m],
        ['col_density', 1 * u.s],
        ['au', 1 * u.radian],
        ['rh', 1 * u.radian],
        ['V', 1 * u.m],
        ['surfbright', 1 * u.s],
        ['frac_illum', 1 * u.m],
        ['area_3sigma', 1 * u.m**2],
        ['temperature', 273 * u.s],
        ['period', 1 * u.radian],
        ['beta', 1 * u.s],
        ['delta-v', 1 / u.s],
        ['epoch', Time('2021-09-25').jd]
    )
)
def test_verify_fields_error(field, quantity):
    with pytest.raises(FieldError):
        data = DataClass.from_dict({field: quantity})
        # explicitly call for verification
        data.verify_fields()


def test_alternative_name_uniqueness():
    """test the uniqueness of alternative field names"""

    assert (len(sum(Conf.fieldnames, [])) ==
            len(set(sum(Conf.fieldnames, []))))

    storage = (deepcopy(Conf.fieldnames), deepcopy(Conf.fieldname_idx))

    with pytest.raises(AssertionError):
        # repeat existing fieldname should raise Error
        Conf.fieldnames.append(['i'])
        assert (len(sum(Conf.fieldnames, [])) ==
                len(set(sum(Conf.fieldnames, []))))

    # revert changes to Conf.fieldnames
    Conf.fieldnames = storage[0]
    Conf.fieldname_idx = storage[1]


def test_translate_columns_and_contains(monkeypatch):
    """test function that translates column names and the `in` operator"""

    new_fieldnames_info = [
        {
            'fieldnames': ['zz', 'aa'],
            'dimension': dimensions.length
        },
        {
            'fieldnames': ['yy', 'dd'],
            'dimension': dimensions.time
        }
    ]
    new_fieldnames = [['zz', 'aa'], ['yy', 'dd']]
    new_fieldname_idx = {}
    for idx, field in enumerate(new_fieldnames):
        for alt in field:
            new_fieldname_idx[alt] = idx
    monkeypatch.setattr(Conf, "fieldnames_info", new_fieldnames_info)
    monkeypatch.setattr(Conf, "fieldnames", new_fieldnames)
    monkeypatch.setattr(Conf, "fieldname_idx", new_fieldname_idx)

    tab = DataClass.from_dict(
        OrderedDict((('aa', [1, 2, 3]*u.m),
                     ('bb', [4, 5, 6]*u.m/u.s),
                     ('cc', ['a', 'b', 'c']))))

    assert tab._translate_columns(['aa', 'bb', 'cc']) == ['aa', 'bb', 'cc']
    assert tab._translate_columns(['zz', 'bb', 'cc']) == ['aa', 'bb', 'cc']

    with pytest.raises(KeyError):
        tab._translate_columns(['x'])  # undefined column name
        tab._translate_columns(['dd'])  # defined column name but not in table

    trans = tab._translate_columns(['x', 'dd'], ignore_missing=True)
    assert trans == ['x', 'dd']

    assert 'aa' in tab
    assert 'bb' in tab
    assert 'zz' in tab
    assert 'x' not in tab  # undefined column name
    assert 'dd' not in tab  # defined column name but no in table


def test_indexing():
    """make sure that indexing functionality is not compromised through
    column name translation"""

    tab = DataClass.from_dict(
        OrderedDict((('aa', [1, 2, 3]*u.m),
                     ('bb', [4, 5, 6]*u.m/u.s),
                     ('cc', ['a', 'b', 'c']))))

    assert list(tab['aa'].data) == [1, 2, 3]
    assert list(tab['aa', 'bb']['aa'].data) == [1, 2, 3]
    assert len(tab[tab['aa'] < 3*u.m]) == 2


def test_field_conversion():
    """test field conversion functions"""

    tab = DataClass.from_dict(OrderedDict((('d', [1]*u.m),
                                           ('bb', [4]*u.m/u.s),
                                           ('cc', ['a']))))

    assert tab['d'] == 1*u.m
    assert tab['diameter'] == 1*u.m
    assert tab['R'] == 0.5*u.m
    assert tab['radius'] == 0.5*u.m

    tab = DataClass.from_dict(OrderedDict((('R', [1]*u.m),
                                           ('bb', [4]*u.m/u.s),
                                           ('cc', ['a']))))

    assert tab['d'] == 2*u.m
    assert tab['R'] == 1*u.m

    # test something that is not defined anywhere in data.Conf
    with pytest.raises(KeyError):
        tab = DataClass.from_dict(OrderedDict((('aa', [1]*u.m),
                                               ('bb', [4]*u.m/u.s),
                                               ('cc', ['a']))))
        tab['bullpoop']

    # test something that is defined anywhere in data.Conf
    with pytest.raises(KeyError):
        tab = DataClass.from_dict(OrderedDict((('aa', [1]*u.m),
                                               ('bb', [4]*u.m/u.s),
                                               ('cc', ['a']))))
        tab['radius']


def test_modifications():
    """test modifying tables using astropy.table methods"""

    tab = DataClass.from_dict(
        OrderedDict((('aa', [1, 2, 3]*u.m),
                     ('bb', [4, 5, 6]*u.m/u.s),
                     ('cc', ['a', 'b', 'c']))))

    # adding single rows
    tab.table.add_row([4*u.m, 7*u.m/u.s, 'd'])

    with pytest.raises(ValueError):
        tab.table.add_row([4*u.s, 7*u.m/u.s, 'd'])  # fails: wrong unit

    with pytest.raises(ValueError):
        tab.table.add_row([5*u.m, 8*u.m/u.s])  # fails: incomplete

    # adding multiple rows
    data = [[7*u.m, 10*u.m/u.s, 'g'],
            [8*u.m, 11*u.m/u.s, 'h']]
    for dat in data:
        tab.table.add_row(dat)

    assert all(tab['aa']**2 == [1, 4, 9, 16, 49, 64]*u.m*u.m)

    # adding columns
    tab['dd'] = [10, 20, 30, 40, 50, 60]*u.kg/u.um

    assert tab['dd'][5] == 60*u.kg/u.um

    assert len(tab[5]) == 1

    assert len(tab) == 6


def test_unit_apply():
    """test DataClass._unit_apply"""
    assert DataClass._unit_apply(1, None) == 1
    assert DataClass._unit_apply(1*u.km, None) == 1
    assert DataClass._unit_apply(1, u.kg) == 1*u.kg
    assert DataClass._unit_apply(1*u.kg, u.g) == 1000*u.g


def test_meta():
    """test meta data mechanisms"""
    tab = DataClass.from_dict(OrderedDict([('aa', [1, 2, 3]*u.m),
                                           ('bb', [4, 5, 6]),
                                           ('cc', [7, 8, 9]*u.kg)]))
    assert tab.meta == {}

    tab.meta['test'] = 'stuff'

    assert tab.meta == {'test': 'stuff'}


def test_io():
    """test file writing and reading capabilities"""
    tab = DataClass.from_dict(OrderedDict([('aa', [1, 2, 3]*u.m),
                                           ('bb', [4, 5, 6]),
                                           ('cc', [7, 8, 9]*u.kg)]))
    tab.meta['test'] = 'stuff'

    try:
        tab.to_file('dataclass_table.fits', format='fits', overwrite=True)

        tab2 = DataClass.from_file('dataclass_table.fits', format='fits')

        assert all(tab.table == tab2.table)
        assert tab.meta == {key.lower(): val.lower()
                            for key, val in tab2.meta.items()}
    finally:
        os.unlink("dataclass_table.fits")


def test_apply():
    """test DataClass.apply"""
    tab = DataClass.from_columns([[2451223, 2451224, 2451226]*u.d,
                                  [120.1, 121.3, 124.9]*u.deg,
                                  [12.4, 12.2, 10.8]*u.deg],
                                 names=('JD', 'RA', 'DEC'))
    tab.apply([[12.1], [12.5, 12.6], [13.5, 13.4, 13.5]],
              name='V', unit='mag')

    assert len(tab) == 6

    tab = DataClass.from_columns([[2451223, 2451224, 2451226]*u.d,
                                  [120.1, 121.3, 124.9]*u.deg,
                                  [12.4, 12.2, 10.8]*u.deg],
                                 names=('JD', 'RA', 'DEC'))
    tab.apply([12.1, 12.5, 12.6]*u.mag, name='V')

    assert len(tab) == 3

    with pytest.raises(DataClassError):
        tab.apply([12.1, 12.5, 12.6, 99]*u.mag, name='V')  # wrong size


def test_add_row():
    """test DataClass.add_row"""
    tab = DataClass.from_columns([Time([2451223, 2451224, 2451226],
                                       format='jd'),
                                  [120.1, 121.3, 124.9]*u.deg,
                                  [12.4, 12.2, 10.8]*u.deg],
                                 names=('JD', 'RA', 'DEC'))
    # add astropy Row
    r = tab.table[0]
    assert isinstance(r, Row)
    tab.add_row(r)
    assert len(tab) == 4
    assert tab.table[-1] == r

    # add a dict
    r = {'JD': 2451228 * u.d, 'RA': 130 * u.deg, 'DEC': 8 * u.deg}
    tab.add_row(r)
    assert len(tab) == 5
    assert all(tab[-1]['JD'] == Time(r['JD'], format='jd'))
    for k in ['RA', 'DEC']:
        assert u.isclose(tab[-1][k], r[k])

    # add an iterable that matches the existing columns
    r = [2451130 * u.d, 135 * u.deg, 6 * u.deg]
    tab.add_row(r)
    assert len(tab) == 6
    assert all(tab[-1]['JD'] == Time(r[0], format='jd'))
    for i, k in enumerate(tab.field_names[1:]):
        assert u.isclose(tab[-1][k], r[i+1])

    # add an iterable with specified column names
    r = [2451132 * u.d, 140 * u.deg, 3 * u.au]
    n = ['JD', 'RA', 'rh']  # with a new column and missing an existing column
    tab.add_row(r, n)
    assert len(tab) == 7
    assert set(tab.field_names) == {'JD', 'RA', 'DEC', 'rh'}
    assert all(tab[-1]['JD'] == Time(r[0], format='jd'))
    for i, k in enumerate(n[1:]):
        assert u.isclose(tab[-1][k], r[i+1])

    # time represented by a string
    r = ['1998-11-18', 120 * u.deg, 3 * u.au]
    n = ['JD', 'RA', 'rh']
    tab.add_row(r, n)
    assert len(tab) == 8
    assert all(tab[-1]['JD'] == r[0])

    # specify different names from the Mapping object
    r = {'JD': 2451228 * u.d, 'RA': 130 * u.deg, 'DEC': 8 * u.deg}
    n = ['JD', 'RA', 'phase']
    tab.add_row(r, n)
    assert len(tab) == 9
    assert set(tab.field_names) == {'JD', 'RA', 'DEC', 'rh', 'phase'}
    assert u.isclose(tab[-1]['phase'], r['DEC'])


def test_vstack():
    """test DataClass.vstack"""
    tab = DataClass.from_columns([[2451223, 2451224, 2451226]*u.d,
                                  [120.1, 121.3, 124.9]*u.deg,
                                  [12.4, 12.2, 10.8]*u.deg],
                                 names=('JD', 'RA', 'DEC'))

    # join a DataClass, same columns
    assert isinstance(tab, DataClass)
    tab.vstack(tab)
    assert len(tab) == 6
    assert set(tab.field_names) == {'JD', 'RA', 'DEC'}
    assert all(tab.table[:3] == tab.table[-3:])

    # join a Table
    delta_tab = tab.table
    assert isinstance(delta_tab, QTable)
    tab.vstack(delta_tab)
    assert len(tab) == 12
    assert set(tab.field_names) == {'JD', 'RA', 'DEC'}
    assert all(tab.table[:6] == tab.table[-6:])

    # join a dict
    delta_tab = dict(tab.table)
    assert isinstance(delta_tab, dict)
    tab.vstack(dict(delta_tab))
    assert len(tab) == 24
    assert set(tab.field_names) == {'JD', 'RA', 'DEC'}
    assert all(tab.table[:6] == tab.table[-6:])

    # join an unrecoganized object
    with pytest.raises(ValueError):
        tab.vstack([1, 2, 3])

    # join a table with different sets of columns
    tab = DataClass.from_columns([[2451223, 2451224, 2451226]*u.d,
                                  [120.1, 121.3, 124.9]*u.deg,
                                  [12.4, 12.2, 10.8]*u.deg],
                                 names=('JD', 'RA', 'DEC'))
    subtab = QTable([[1, 2, 3] * u.au,
                     [1, 2, 3] * u.au,
                     [20, 30, 40] * u.deg],
                    names=('r', 'delta', 'DEC'))
    field0 = tab.field_names
    tab.vstack(subtab)
    assert len(tab) == 6
    assert set(field0).union(set(subtab.colnames)) == set(tab.field_names)

    # join a table that has a column using alternative names
    subtab = QTable([[4, 5] * u.au,
                     [10, 20] * u.deg],
                    names=('rh', 'phase'))
    field0 = tab.field_names
    tab.vstack(subtab)
    assert len(tab) == 8
    assert 'rh' not in tab.table.colnames
    assert set(field0).union({'phase'}) == set(tab.field_names)
