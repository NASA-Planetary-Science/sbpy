# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
================
SBPy Data Module
================

created on June 22, 2017
"""

__all__ = ['DataClass', 'mpc_observations', 'sb_search', 'image_search',
           'pds_ferret']

from numpy import ndarray, array
from astropy.table import QTable, Column


class DataClass():
    """`~sbpy.data.DataClass` serves as the base class for all data
    container classes in ``sbpy`` in order to provide consistent
    functionality throughout all these classes.

    The core of `~sbpy.data.DataClass` is an `~astropy.table.QTable`
    object (referred to as the `data table` below) - a type of
    `~astropy.table.Table` object that supports the `~astropy.units`
    formalism on a per-column base - which already provides most of
    the required functionality. `~sbpy.data.DataClass` objects can be
    manually generated from ``dict``
    (`~sbpy.data.DataClass.from_dict`), `~numpy.array`-like
    (`~sbpy.data.DataClass.from_array`) objects, or directly from
    another `~astropy.table.QTable` object.

    A few high-level functions for table data access or modification
    are provided; other, more complex modifications can be applied to
    the table object (`~sbpy.data.DataClass.data`) directly.

    """

    def __init__(self, **kwargs):
        """``__init__``: Build data table from ``**kwargs``."""
        self.table = QTable()
        # self.altkeys = {}  # dictionary for alternative column names

        if (len(kwargs.items()) == 1 and 'table' in kwargs.keys()):
            # single item provided named 'table' -> already Table object
            self.table = kwargs['table']
        else:
            # treat kwargs as dictionary
            for key, val in kwargs.items():
                try:
                    unit = val.unit
                    val = val.value
                except AttributeError:
                    unit = None

                # check if val is already list-like
                try:
                    val[0]
                except TypeError:
                    val = [val]

                self.table[key] = Column(val, unit=unit)

    @classmethod
    def from_dict(cls, data):
        """Create `~sbpy.data.DataClass` object from dictionary or list of
        dictionaries.

        Parameters
        ----------
        data : dictionary or list (or similar) of dictionaries
             Data that will be ingested in `~sbpy.data.DataClass` object.
             Each dictionary creates a row in the data table. Dictionary
             keys are used as column names; corresponding values must be
             scalar (cannot be lists or arrays). If a list of dicitionaries
             is provided, all dictionaries have to provide the same
             set of keys (and units, if used at all).

        Returns
        -------
        `DataClass` object

        Examples
        --------
        >>> import astropy.units as u
        >>> from sbpy.data import Orbit
        >>> orb = Orbit.from_dict({'a': 2.7674*u.au,
        ...                        'e': .0756,
        ...                        'i': 10.59321*u.deg})
        >>> print(orb)
        <sbpy.data.orbit.Orbit object at 0x...>
        >>> print(orb.column_names) # doctest: +SKIP
        <TableColumns names=('a','e','i')>
        >>> print(orb.data['a', 'e', 'i'])
          a      e       i
          AU            deg
        ------ ------ --------
        2.7674 0.0756 10.59321

        """
        if isinstance(data, dict):
            return cls(**data)
        elif isinstance(data, (list, ndarray, tuple)):
            # build table from first dict and append remaining rows
            tab = cls(**data[0])
            for row in data[1:]:
                tab.add_rows(row)
            return tab
        else:
            raise TypeError('this function requires a dictionary or a '
                            'list of dictionaries')

    @classmethod
    def from_array(cls, data, names):
        """Create `~sbpy.data.DataClass` object from list or array.

        Parameters
        ----------
        data : 1d or 2d list-like
             Data that will be ingested in `DataClass` object. Each
             of the sub-lists or sub-arrays on the 0-axis creates a row
             in the data table. Dictionary
             keys are used as column names. If a list of dicitionaries
             is provided, all dictionaries have to provide the same
             set of keys (and units, if used at all).
        data : list or array, mandatory
            data that will be rearranged in astropy `Table` format, one
            array per column
        names : list, mandatory
            column names, must have n names for n `data` arrays

        Returns
        -------
        `DataClass` object

        Examples
        --------
        #>>> import astropy.units as u
        #>>> from sbpy.data import Orbit
        #>>> from numpy.random import random as r
        #>>> orb = Orbit.from_array(data=[r(100)*2*u.au,
        #>>>                              r(100),
        #>>>                              r(100)*180*u.deg],
        #>>>                        names=['a', 'e', 'i'])
        """

        return cls.from_dict(dict(zip(names, data)))

    @classmethod
    def from_table(cls, data):
        """Create `DataClass` object from astropy `Table` object.

        Parameters
        ----------
        data : astropy `Table` object, mandatory
             Data that will be ingested in `DataClass` object.

        Returns
        -------
        `DataClass` object

        """

        return cls(table=data)

    def __getattr__(self, field):
        if field in self.table.columns:
            return self.table[field]
        else:
            raise AttributeError("field '{:s}' does not exist".format(field))

    def __setattr__(self, field, value):
        """Set attribute modify attribute in `self.table`, if available, else set it for self
        """
        try:
            # if self.table exists ...
            if field in self.table.columns:
                # set value there...
                self.table[field] = value
            else:
                super().__setattr__(field, value)
        except:
            # if, not set it for self
            super().__setattr__(field, value)

    def __getitem__(self, ident):
        """Return column or row from data table"""
        return self.table[ident]

    @property
    def data(self):
        """returns the `~astropy.table.QTable` containing all data."""
        return self.table

    @property
    def column_names(self):
        """Returns a list of all column names in the data table"""
        return self.table.columns

    def add_rows(self, rows):
        """Appends additional rows to the existing data table. Must be in the
        form of a list, tuple, or `~numpy.ndarray` of rows or a single
        list, tuple, `~numpy.ndarray`, or dictionary to the current
        data table. The new data rows can each be provided in the form
        of a dictionary or a list. In case of a dictionary, all table
        column names must be provided in ``row``; additional keys that
        are not yet column names in the table will be discarded. In
        case of a list, the list elements must be in the same order as
        the table columns. In either case, `~astropy.units` must be
        provided in ``rows`` if used in the data table.

        Parameters
        ----------
        rows : list, tuple, `~numpy.ndarray`, or dict
            data to be appended to the table; required to have the same 
            length as the existing table, as well as the same units


        Returns
        -------
        n : int, the total number of rows in the data table

        """
        if isinstance(rows, dict):
            try:
                newrow = [rows[colname] for colname in self.table.columns]
            except KeyError as e:
                raise ValueError('data for column {0} missing in row {1}'.
                                 format(e, rows))
            self.add_rows(newrow)
        if isinstance(rows, (list, ndarray, tuple)):
            if len(array(rows).shape) > 1 or isinstance(rows[0], dict):
                for subrow in rows:
                    self.add_rows(subrow)
            else:
                self.table.add_row(rows)
        return len(self.table)

    def add_column(self, data, name):
        """Append a single column to the current data table. The lenght of 
        the input list, `~numpy.ndarray`, or tuple must match the current 
        number of rows in the data table.

        Parameters
        ----------
        data : list, `~numpy.ndarray`, or tuple
            data to be filled into the table; required to have the same 
            length as the existing table
        name : string, new column's name

        Returns
        -------
        n : int, the total number of columns in the data table

        """
        self.table.add_column(Column(data, name=name))
        return len(self.column_names)

    def _check_columns(self, colnames):
        """Checks whether all of the elements in colnames exist as
        column names in the data table."""

        return all([col in self.column_names for col in colnames])


def mpc_observations(targetid):
    """Obtain all available observations of a small body from the Minor
    Planet Center (http://www.minorplanetcenter.net) and provides them in
    the form of an Astropy table.

    Parameters
    ----------
    targetid : str, mandatory
        target identifier

    Returns
    -------
    `~sbpy.data.DataClass` object

    Examples
    --------
    >>> from sbpy.data import mpc_observations
    >>> obs = mpc_observations('ceres')  # doctest: +SKIP

    not yet implemented

    """


def sb_search(field):
    """Use the Skybot service (http://vo.imcce.fr/webservices/skybot/) at
    IMCCE to Identify moving objects potentially present in a registered
    FITS images.

    Parameters
    ----------
    field : string, astropy.io.fits Header object, Primary HDU, or Image HDU
      A FITS image file name, HDU data structure, or header with
      defined WCS

    Returns
    -------
    `~sbpy.data.DataClass` object

    Examples
    --------
    >>> from sbpy.data import sb_search
    >>> objects = sb_search('ceres')  # doctest: +SKIP

    not yet implemented

    """


def image_search(targetid):
    """Use the Solar System Object Image Search function of the Canadian
    Astronomy Data Centre
    (http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/ssois/) to identify
    images with a specific small body in them.

    Parameters
    ----------
    targetid : str, mandatory
        target identifier

    Returns
    -------
    `~sbpy.data.DataClass` object

    Examples
    --------
    >>> from sbpy.data import Misc  # doctest: +SKIP
    >>> images = Misc.image_search('ceres')  # doctest: +SKIP

    not yet implemented

    """


def pds_ferret(targetid):
    """Use the Small Bodies Data Ferret (http://sbntools.psi.edu/ferret/)
    at the Planetary Data System's Small Bodies Node to query for
    information on a specific small body in the PDS.

    Parameters
    ----------
    targetid : str, mandatory
        target identifier

    Returns
    -------
    data : dict
      A hierarchical data object

    Examples
    --------
    >>> from sbpy.data import pds_ferret
    >>> data = pds_ferret('ceres')  # doctest: +SKIP

    not yet implemented

    """
