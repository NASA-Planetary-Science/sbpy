# Licensed under a 3-clause BSD style license - see LICENSE.rst


def test_creation_single():
    """ test the creation of DataClass objects from dicts or arrays;
    single row only"""

    from astropy.table import QTable
    from ..core import DataClass

    ground_truth = QTable([[1], [2], ['test']], names=('a', 'b', 'c'))

    test_init = DataClass(a=1, b=2, c='test')
    assert test_init.table == ground_truth

    test_dict = DataClass.from_dict({'a': 1, 'b': 2, 'c': 'test'})
    assert test_dict.table == ground_truth

    test_array = DataClass.from_array([1, 2, 'test'], names=('a', 'b', 'c'))
    assert test_array.table == ground_truth


def test_creation_multi():
    """ test the creation of DataClass objects from dicts or arrays;
    multiple rows"""

    import pytest
    from astropy.table import QTable
    from ..core import DataClass

    ground_truth = QTable([[1, 2, 3], [4, 5, 6], ['a', 'b', 'c']],
                          names=('a', 'b', 'c'))

    test_dict = DataClass.from_dict([{'a': 1, 'b': 4, 'c': 'a'},
                                     {'a': 2, 'b': 5, 'c': 'b'},
                                     {'a': 3, 'b': 6, 'c': 'c'}])
    assert all(test_dict.table == ground_truth)

    test_array = DataClass.from_array([[1, 2, 3], [4, 5, 6], ['a', 'b', 'c']],
                                      names=('a', 'b', 'c'))
    assert all(test_array.table == ground_truth)

    test_table = DataClass.from_table(ground_truth)
    assert all(test_table.table == ground_truth)

    # test failing if columns have different lengths
    with pytest.raises(ValueError):
        test_dict = DataClass.from_dict([{'a': 1, 'b': 4, 'c': 'a'},
                                         {'a': 2, 'b': 5, 'c': 'b'},
                                         {'a': 3, 'b': 6}])

    with pytest.raises(ValueError):
        test_array = DataClass.from_array([[1, 2, 3], [4, 5, 6], ['a', 'b']],
                                          names=('a', 'b', 'c'))


def test_units():
    """ test units on multi-row tables """

    from astropy.table import QTable
    import astropy.units as u
    from ..core import DataClass

    ground_truth = QTable([[1, 2, 3]*u.Unit('m'),
                           [4, 5, 6]*u.m/u.s,
                           ['a', 'b', 'c']],
                          names=('a', 'b', 'c'))

    assert ((ground_truth['a']**2).unit == 'm2')

    test_dict = DataClass.from_dict([{'a': 1*u.m, 'b': 4*u.m/u.s, 'c': 'a'},
                                     {'a': 2*u.m, 'b': 5*u.m/u.s, 'c': 'b'},
                                     {'a': 3*u.m, 'b': 6*u.m/u.s, 'c': 'c'}])
    assert all(test_dict.table == ground_truth)

    test_array = DataClass.from_array([[1, 2, 3]*u.m,
                                       [4, 5, 6]*u.m/u.s,
                                       ['a', 'b', 'c']],
                                      names=('a', 'b', 'c'))
    assert all(test_array.table == ground_truth)


def test_add():
    """ test adding rows and columns to an existing table """

    import pytest
    from numpy import array
    from astropy.table import QTable
    import astropy.units as u
    from ..core import DataClass

    tab = DataClass.from_dict([{'a': 1*u.m, 'b': 4*u.m/u.s, 'c': 'a'},
                               {'a': 2*u.m, 'b': 5*u.m/u.s, 'c': 'b'},
                               {'a': 3*u.m, 'b': 6*u.m/u.s, 'c': 'c'}])

    # adding single rows

    tab.add_rows([4*u.m, 7*u.m/u.s, 'd'])

    with pytest.raises(ValueError):
        tab.add_rows([4*u.s, 7*u.m/u.s, 'd'])  # fails: wrong unit

    tab.add_rows({'a': 5*u.m, 'b': 8*u.m/u.s, 'c': 'e'})

    with pytest.raises(ValueError):
        tab.add_rows({'a': 5*u.m, 'b': 8*u.m/u.s})  # fails: incomplete

    # ignore superfluent columns
    tab.add_rows({'a': 6*u.m, 'b': 9*u.m/u.s, 'c': 'f', 'd': 'not existent'})

    # adding multiple rows

    tab.add_rows(([7*u.m, 10*u.m/u.s, 'g'],
                  [8*u.m, 11*u.m/u.s, 'h']))

    tab.add_rows([{'a': 9*u.m, 'b': 12*u.m/u.s, 'c': 'i'},
                  {'a': 10*u.m, 'b': 13*u.m/u.s, 'c': 'j'}])

    assert all(tab['a']**2 == [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]*u.m*u.m)

    # adding columns

    tab.add_column(array([10, 20, 30, 40, 50,
                          60, 70, 80, 90, 100])*u.kg/u.um, name='d')

    assert tab[0]['d'] == 10*u.kg/u.um

    with pytest.raises(ValueError):
        tab.add_column(array([10, 20, 30, 40, 50,
                              60, 70, 80, 90])*u.kg/u.um, name='e')


def test_check_columns():
    """test function that checks the existing of a number of column names
    provided"""

    import pytest
    import astropy.units as u
    from ..core import DataClass

    tab = DataClass.from_dict([{'a': 1*u.m, 'b': 4*u.m/u.s, 'c': 'a'},
                               {'a': 2*u.m, 'b': 5*u.m/u.s, 'c': 'b'},
                               {'a': 3*u.m, 'b': 6*u.m/u.s, 'c': 'c'}])

    assert tab._check_columns(['a', 'b', 'c'])
    assert tab._check_columns(['a', 'b'])
    assert tab._check_columns(['a', 'b', 'f']) == False