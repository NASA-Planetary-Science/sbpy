# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from collections import OrderedDict
import pytest
from copy import deepcopy
from numpy import array
import astropy.units as u
from astropy.table import QTable
from ..core import DataClass, conf, DataClassError


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


def test_creation_single_row():
    """ test the creation of DataClass objects from dicts or arrays;
    single row only"""

    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1], [2], [3]], names=('a', 'b', 'c'))
    ground_truth_2 = QTable([[1]*u.m, [2]*u.kg, [3]*u.cm/u.s],
                            names=('a', 'b', 'c'))
    ground_truth_3 = QTable([[1]*u.m, [2], [3]*u.kg], names=('a', 'b', 'c'))
    ground_truth_4 = QTable([[1], ['stuff'], [3]], names=('a', 'b', 'c'))
    ground_truth_5 = QTable([[1]*u.km, [2], ['test']], names=('a', 'b', 'c'))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict(OrderedDict([('a', 1), ('b', 2),
                                                   ('c', 3)]))
    test_dict_2 = DataClass.from_dict(OrderedDict([('a', 1*u.m),
                                                   ('b', 2*u.kg),
                                                   ('c', 3*u.cm/u.s)]))
    test_dict_3 = DataClass.from_dict(OrderedDict([('a', 1*u.m), ('b', 2),
                                                   ('c', 3*u.kg)]))
    test_dict_4 = DataClass.from_dict(OrderedDict([('a', 1), ('b', 'stuff'),
                                                   ('c', 3)]))
    test_dict_5 = DataClass.from_dict(OrderedDict([('a', 1*u.km), ('b', 2),
                                                   ('c', 'test')]))

    assert test_dict_1.table == ground_truth_1
    assert test_dict_2.table == ground_truth_2
    assert test_dict_3.table == ground_truth_3
    assert test_dict_4.table == ground_truth_4
    assert test_dict_5.table == ground_truth_5

    # test DataClass.from_rows for different cases
    test_array_1 = DataClass.from_rows([1, 2, 3],
                                       names=('a', 'b', 'c'))
    test_array_2 = DataClass.from_rows([1*u.m, 2*u.kg, 3*u.cm/u.s],
                                       names=('a', 'b', 'c'))
    test_array_3 = DataClass.from_rows([1*u.m, 2, 3*u.kg],
                                       names=('a', 'b', 'c'))
    test_array_4 = DataClass.from_rows([1, 'stuff', 3],
                                       names=('a', 'b', 'c'))
    test_array_5 = DataClass.from_rows([1*u.km, 2, 'test'],
                                       names=('a', 'b', 'c'))

    assert test_array_1.table == ground_truth_1
    assert test_array_2.table == ground_truth_2
    assert test_array_3.table == ground_truth_3
    assert test_array_4.table == ground_truth_4
    assert test_array_5.table == ground_truth_5

    # test single row, single column
    ground_truth_6 = QTable([[1]], names=('a'))

    test_dict_6 = DataClass.from_dict({'a': 1})
    assert test_dict_6.table == ground_truth_6

    test_array_6 = DataClass.from_rows([1], names='a')
    assert test_array_6.table == ground_truth_6

    # test units parameter
    test_array_2b = DataClass.from_rows([1, 2, 3],
                                        names=('a', 'b', 'c'),
                                        units=(u.m, u.kg, u.cm/u.s))
    test_array_3b = DataClass.from_rows([1, 2, 3],
                                        names=('a', 'b', 'c'),
                                        units=(u.m, None, u.kg))

    assert test_array_2b.table == ground_truth_2
    assert test_array_3b.table == ground_truth_3

    ground_truth_7 = QTable([[1]*u.kg], names=('a'))
    test_array_7 = DataClass.from_rows([1], names='a', units='kg')
    assert test_array_7.table == ground_truth_7

    with pytest.raises(DataClassError):
        DataClass.from_rows([1], names='a', units=('m', 'kg'))

    # test single row starting with string
    ground_truth_8 = QTable(rows=[['a', 1]], names=('a', 'b'))
    test_array_8 = DataClass.from_rows(['a', 1], names=('a', 'b'))

    assert ground_truth_8 == test_array_8.table


def test_creation_multi_rows():
    """test the creation of DataClass objects from dicts or arrays;
    multiple rows"""
    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            names=('a', 'b', 'c'))
    ground_truth_2 = QTable([[1, 2, 3]*u.m, [4, 5, 6]*u.kg,
                             [7, 8, 9]*u.cm/u.s], names=('a', 'b', 'c'))
    ground_truth_3 = QTable([[1, 2, 3]*u.m, [4, 5, 6],
                             [7, 8, 9]*u.kg], names=('a', 'b', 'c'))
    ground_truth_4 = QTable([[1, 2, 3], ['a', 'b', 'c'], [7, 8, 9]],
                            names=('a', 'b', 'c'))
    ground_truth_5 = QTable([[1, 2, 3]*u.km, [4, 5, 6], ['a', 'b', 'c']],
                            names=('a', 'b', 'c'))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict(OrderedDict([('a', [1, 2, 3]),
                                                   ('b', [4, 5, 6]),
                                                   ('c', [7, 8, 9])]))
    test_dict_2 = DataClass.from_dict(OrderedDict(
        [('a', [1, 2, 3]*u.m), ('b', [4, 5, 6]*u.kg),
         ('c', [7, 8, 9]*u.cm/u.s)]))
    test_dict_3 = DataClass.from_dict(OrderedDict([('a', [1, 2, 3]*u.m),
                                                   ('b', [4, 5, 6]),
                                                   ('c', [7, 8, 9]*u.kg)]))
    test_dict_4 = DataClass.from_dict(OrderedDict([('a', [1, 2, 3]),
                                                   ('b', ['a', 'b', 'c']),
                                                   ('c', [7, 8, 9])]))
    test_dict_5 = DataClass.from_dict(OrderedDict([('a', [1, 2, 3]*u.km),
                                                   ('b', [4, 5, 6]),
                                                   ('c', ['a', 'b', 'c'])]))

    assert all(test_dict_1.table == ground_truth_1)
    assert all(test_dict_2.table == ground_truth_2)
    assert all(test_dict_3.table == ground_truth_3)
    assert all(test_dict_4.table == ground_truth_4)
    assert all(test_dict_5.table == ground_truth_5)

    # test DataClass.from_rows for different cases
    test_array_1 = DataClass.from_rows([[1, 4, 7], [2, 5, 8], [3, 6, 9]],
                                       names=('a', 'b', 'c'))
    test_array_2 = DataClass.from_rows([[1*u.m, 4*u.kg, 7*u.cm/u.s],
                                        [2*u.m, 5*u.kg, 8*u.cm/u.s],
                                        [3*u.m, 6*u.kg, 9*u.cm/u.s]],
                                       names=('a', 'b', 'c'))
    test_array_3 = DataClass.from_rows([[1*u.m, 4, 7*u.kg],
                                        [2*u.m, 5, 8*u.kg],
                                        [3*u.m, 6, 9*u.kg]],
                                       names=('a', 'b', 'c'))
    test_array_4 = DataClass.from_rows([[1, 'a', 7],
                                        [2, 'b', 8],
                                        [3, 'c', 9]],
                                       names=('a', 'b', 'c'))
    test_array_5 = DataClass.from_rows([[1*u.km, 4, 'a'],
                                        [2*u.km, 5, 'b'],
                                        [3*u.km, 6, 'c']],
                                       names=('a', 'b', 'c'))

    assert all(test_array_1.table == ground_truth_1)
    assert all(test_array_2.table == ground_truth_2)
    assert all(test_array_3.table == ground_truth_3)
    assert all(test_array_4.table == ground_truth_4)
    assert all(test_array_5.table == ground_truth_5)


def test_creation_single_column():
    """test the creation of DataClass objects from dicts or arrays;
    single column only"""

    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1, 2, 3]], names=('a'))
    ground_truth_2 = QTable([[1, 2, 3]*u.kg], names=('a'))
    ground_truth_3 = QTable([['a', 'b', 'c']], names=('a'))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict({'a': [1, 2, 3]})
    test_dict_2 = DataClass.from_dict({'a': [1, 2, 3]*u.kg})
    test_dict_3 = DataClass.from_dict({'a': ['a', 'b', 'c']})

    assert all(test_dict_1.table == ground_truth_1)
    assert all(test_dict_2.table == ground_truth_2)
    assert all(test_dict_3.table == ground_truth_3)

    # test DataClass.from_columns for different cases
    test_array_1 = DataClass.from_columns([1, 2, 3], names='a')
    test_array_2 = DataClass.from_columns([1, 2, 3]*u.kg, names='a')
    test_array_3 = DataClass.from_columns(['a', 'b', 'c'], names='a')

    assert all(test_array_1.table == ground_truth_1)
    assert all(test_array_2.table == ground_truth_2)
    assert all(test_array_3.table == ground_truth_3)

    # test single row, single column
    ground_truth_4 = QTable([[1]], names=('a'))

    test_dict_4 = DataClass.from_dict({'a': 1})
    assert test_dict_4.table == ground_truth_4

    test_array_4 = DataClass.from_columns([1], names='a')
    assert test_array_4.table == ground_truth_4


def test_creation_multi_column():
    """test the creation of DataClass objects from dicts or arrays;
    multiple columns"""

    # ground truth tables - compare against these tables
    ground_truth_1 = QTable([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            names=('a', 'b', 'c'))
    ground_truth_2 = QTable([[1, 2, 3]*u.kg, [4, 5, 6]*u.m/u.s],
                            names=('a', 'b'))
    ground_truth_3 = QTable([[1, 2, 3], [4, 5, 6]*u.m/u.s],
                            names=('a', 'b'))
    ground_truth_4 = QTable([[1, 2, 3], ['a', 'b', 'c']],
                            names=('a', 'b'))
    ground_truth_5 = QTable([[1, 2, 3], [4, 5, 6]*u.s/u.kg, ['a', 'b', 'c']],
                            names=('a', 'b', 'c'))

    # test DataClass.from_dict for different cases
    test_dict_1 = DataClass.from_dict({'a': [1, 2, 3], 'b': [4, 5, 6],
                                       'c': [7, 8, 9]})
    test_dict_2 = DataClass.from_dict({'a': [1, 2, 3]*u.kg,
                                       'b': [4, 5, 6]*u.m/u.s})
    test_dict_3 = DataClass.from_dict({'a': [1, 2, 3],
                                       'b': [4, 5, 6]*u.m/u.s})
    test_dict_4 = DataClass.from_dict({'a': [1, 2, 3],
                                       'b': ['a', 'b', 'c']})
    test_dict_5 = DataClass.from_dict({'a': [1, 2, 3],
                                       'b': [4, 5, 6]*u.s/u.kg,
                                       'c': ['a', 'b', 'c']})
    assert all(test_dict_1.table == ground_truth_1)
    assert all(test_dict_2.table == ground_truth_2)
    assert all(test_dict_3.table == ground_truth_3)
    assert all(test_dict_4.table == ground_truth_4)
    assert all(test_dict_5.table == ground_truth_5)

    # test DataClass.from_columns for different cases
    test_array_1 = DataClass.from_columns([[1, 2, 3], [4, 5, 6],
                                           [7, 8, 9]], names=('a', 'b', 'c'))
    test_array_2 = DataClass.from_columns([[1, 2, 3]*u.kg, [4, 5, 6]*u.m/u.s],
                                          names=('a', 'b'))
    test_array_3 = DataClass.from_columns([[1, 2, 3], [4, 5, 6]*u.m/u.s],
                                          names=('a', 'b'))
    test_array_4 = DataClass.from_columns([[1, 2, 3], ['a', 'b', 'c']],
                                          names=('a', 'b'))
    test_array_5 = DataClass.from_columns([[1, 2, 3], [4, 5, 6]*u.s/u.kg,
                                           ['a', 'b', 'c']],
                                          names=('a', 'b', 'c'))
    assert all(test_array_1.table == ground_truth_1)
    assert all(test_array_2.table == ground_truth_2)
    assert all(test_array_3.table == ground_truth_3)
    assert all(test_array_4.table == ground_truth_4)
    assert all(test_array_5.table == ground_truth_5)


def test_get_set():
    """ test the get and set methods"""

    data = DataClass.from_dict(
        OrderedDict((('a', [1, 2, 3]),
                     ('b', [4, 5, 6]),
                     ('c', [7, 8, 8]))))

    # get a single column
    assert len(data['a']) == 3

    # mask rows
    masked = data[[True, False, False]]
    assert len(masked) == 1
    assert masked['b'][0] == 4

    # get list of rows
    shortened = data[[0, 1]]
    assert len(shortened) == 2

    # modify an existing column
    data['a'][:] = [0, 0, 0]
    assert data['a'][0] == 0

    with pytest.raises(KeyError):
        data['d']

    # add non-existing column using set
    data['z'] = 3
    assert len(data['z'] == 3)

    # modify existing column using set
    data['z'] = 2
    assert data['z'][1] == 2


def test_units():
    """ test units on multi-row tables """

    ground_truth = QTable([[1, 2, 3]*u.Unit('m'),
                           [4, 5, 6]*u.m/u.s,
                           ['a', 'b', 'c']],
                          names=('a', 'b', 'c'))

    assert ((ground_truth['a']**2).unit == 'm2')

    test_dict = DataClass.from_dict(
        OrderedDict((('a', [1, 2, 3]*u.m),
                     ('b', [4, 5, 6]*u.m/u.s),
                     ('c', ['a', 'b', 'c']))))
    assert all(test_dict.table == ground_truth)

    test_array = DataClass.from_columns([[1, 2, 3]*u.m,
                                         [4, 5, 6]*u.m/u.s,
                                         ['a', 'b', 'c']],
                                        names=('a', 'b', 'c'))
    assert all(test_array.table == ground_truth)


def test_alternative_name_uniqueness():
    """test the uniqueness of alternative field names"""
    from ..core import conf

    assert (len(sum(conf.fieldnames, [])) ==
            len(set(sum(conf.fieldnames, []))))

    storage = (deepcopy(conf.fieldnames), deepcopy(conf.fieldname_idx))

    with pytest.raises(AssertionError):
        # repeat existing fieldname should raise Error
        conf.fieldnames.append(['i'])
        assert (len(sum(conf.fieldnames, [])) ==
                len(set(sum(conf.fieldnames, []))))

    # revert changes to conf.fieldnames
    conf.fieldnames = storage[0]
    conf.fieldname_idx = storage[1]


def test_translate_columns():
    """test function that translates column names"""

    storage = (deepcopy(conf.fieldnames), deepcopy(conf.fieldname_idx))
    conf.fieldnames = [['z', 'a']]
    conf.fieldname_idx = {}
    for idx, field in enumerate(conf.fieldnames):
        for alt in field:
            conf.fieldname_idx[alt] = idx

    tab = DataClass.from_dict(
        OrderedDict((('a', [1, 2, 3]*u.m),
                     ('b', [4, 5, 6]*u.m/u.s),
                     ('c', ['a', 'b', 'c']))))

    assert tab._translate_columns(['a', 'b', 'c']) == ['a', 'b', 'c']
    assert tab._translate_columns(['z', 'b', 'c']) == ['a', 'b', 'c']

    with pytest.raises(KeyError):
        tab._translate_columns(['x'])

    # revert changes to conf.fieldnames
    conf.fieldnames = storage[0]
    conf.fieldname_idx = storage[1]


def test_indexing():
    """make sure that indexing functionality is not compromised through
    column name translation"""

    tab = DataClass.from_dict(
        OrderedDict((('a', [1, 2, 3]*u.m),
                     ('b', [4, 5, 6]*u.m/u.s),
                     ('c', ['a', 'b', 'c']))))

    assert list(tab['a'].data) == [1, 2, 3]
    assert list(tab['a', 'b']['a'].data) == [1, 2, 3]
    assert len(tab[tab['a'] < 3*u.m]) == 2


def test_field_conversion():
    """test field conversion functions"""

    tab = DataClass.from_dict(OrderedDict((('d', [1]*u.m),
                                           ('b', [4]*u.m/u.s),
                                           ('c', ['a']))))

    assert tab['d'] == 1*u.m
    assert tab['diameter'] == 1*u.m
    assert tab['R'] == 0.5*u.m
    assert tab['radius'] == 0.5*u.m

    tab = DataClass.from_dict(OrderedDict((('R', [1]*u.m),
                                           ('b', [4]*u.m/u.s),
                                           ('c', ['a']))))

    assert tab['d'] == 2*u.m
    assert tab['R'] == 1*u.m

    # test something that is not defined anywhere in data.conf
    with pytest.raises(KeyError):
        tab = DataClass.from_dict(OrderedDict((('a', [1]*u.m),
                                               ('b', [4]*u.m/u.s),
                                               ('c', ['a']))))
        tab['bullpoop']

    # test something that is defined anywhere in data.conf
    with pytest.raises(KeyError):
        tab = DataClass.from_dict(OrderedDict((('a', [1]*u.m),
                                               ('b', [4]*u.m/u.s),
                                               ('c', ['a']))))
        tab['radius']


def test_modifications():
    """test modifying tables using astropy.table methods"""

    tab = DataClass.from_dict(
        OrderedDict((('a', [1, 2, 3]*u.m),
                     ('b', [4, 5, 6]*u.m/u.s),
                     ('c', ['a', 'b', 'c']))))

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

    assert all(tab['a']**2 == [1, 4, 9, 16, 49, 64]*u.m*u.m)

    # adding columns
    tab['d'] = [10, 20, 30, 40, 50, 60]*u.kg/u.um

    assert tab['d'][5] == 60*u.kg/u.um

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
    tab = DataClass.from_dict(OrderedDict([('a', [1, 2, 3]*u.m),
                                           ('b', [4, 5, 6]),
                                           ('c', [7, 8, 9]*u.kg)]))
    assert tab.meta == {}

    tab.meta['test'] = 'stuff'

    assert tab.meta == {'test': 'stuff'}


def test_io():
    """test file writing and reading capabilities"""
    tab = DataClass.from_dict(OrderedDict([('a', [1, 2, 3]*u.m),
                                           ('b', [4, 5, 6]),
                                           ('c', [7, 8, 9]*u.kg)]))
    tab.meta['test'] = 'stuff'

    tab.to_file('dataclass_table.fits', format='fits', overwrite=True)

    tab2 = DataClass.from_file('dataclass_table.fits', format='fits')

    assert all(tab.table == tab2.table)
    assert tab.meta == {key.lower(): val.lower()
                        for key, val in tab2.meta.items()}


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
