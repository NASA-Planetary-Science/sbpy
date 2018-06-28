# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
================
SBPy Data Module
================

created on June 22, 2017
"""

__all__ = ['DataClass', 'mpc_observations', 'sb_search', 'image_search',
           'pds_ferret']

from astropy.table import Table, Column


class DataClass():
    """`DataClass` serves as a base class for all data container classes
    in `sbpy` in order to provide consistent functionality for all
    these classes.

    The core of `DataClass` is an `astropy.Table` object
    (`DataClass.table`), which already provides most of the required
    functionality. `DataClass` objects can be manually generated from
    `dict` (DataClass.from_dict), `array`-like
    (DataClass.from_array) objects, or astropy `Table` objects. A few 
    high-level functions for
    table modification are provided; other modifications can be
    applied to the table object (`DataClass.table`) directly.
    """

    def __init__(self, **kwargs):
        """Build data table from `**kwargs`."""
        self.table = Table()
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
        """Create `DataClass` object from dictionary or list of
        dictionaries.

        Parameters
        ----------
        data : dictionary or list of dictionaries
             Data that will be ingested in `DataClass` object. Each
             dictionary creates a row in the data table. Dictionary
             keys are used as column names. If a list of dicitionaries
             is provided, all dictionaries have to provide them same
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
        elif isinstance(data, list):
            # build table from first dict and append remaining rows
            tab = cls(**data[0])
            for row in data[1:]:
                tab.add_row(row)
            return tab
        else:
            raise TypeError('this function requires a dictionary or a '
                            'list of dictionaries')

    @classmethod
    def from_array(cls, data, names):
        """Create `DataClass` object from list or array.

        Parameters
        ----------
        data : list of lists or array
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
        """Set attribute

        modify attribute in `self.table`, if available, else set it for self
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
        """returns the Astropy Table containing all data."""
        return self.table

    @property
    def column_names(self):
        """Returns a list of column names in Table"""
        return self.table.columns

    def add_row(self, row):
        """Append a single row to the current table. The new data can be
        provided in the form of a dictionary or a list. In case of a
        dictionary, all table column names must be provided in row;
        additional keys that are not yet column names in the table
        will be discarded. In case of a list, the list elements must
        be in the same order as the table columns.

        Returns
        -------

        n : int, the total number of rows in the table

        """
        if isinstance(row, dict):
            newrow = [row[colname] for colname in self.table.columns]
            self.add_row(newrow)
        if isinstance(row, list):
            self.table.add_row(row)
        return len(self.data)

    def add_column(self, data, name):
        """Append a single column to the current table.

        Parameters
        ----------
        data : list or array-like, data to be filled into the table; required
        to have the same length as the existing table
        name : string, column name

        Returns
        -------
        None

        """
        self.table.add_column(Column(data, name=name))

    def _check_columns(self, colnames):
        """Checks whether all of the elements in colnames exist as
        column names in `self.table`."""

        return all([col in self.column_names for col in colnames])


def mpc_observations(targetid, bib=None):
    """Obtain all available observations of a small body from the Minor
    Planet Center (http://www.minorplanetcenter.net) and provides them in
    the form of an Astropy table.

    Parameters
    ----------
    targetid : str, mandatory
        target identifier
    bib : SBPy Bibliography instance, optional, default None
        Bibliography instance that will be populated

    Returns
    -------
    Astropy Table

    Examples
    --------
    >>> from sbpy.data import mpc_observations
    >>> obs = mpc_observations('ceres')  # doctest: +SKIP

    not yet implemented

    """


def sb_search(field, bib=None):
    """Use the Skybot service (http://vo.imcce.fr/webservices/skybot/) at
    IMCCE to Identify moving objects potentially present in a registered
    FITS images.

    Parameters
    ----------
    field : string, astropy.io.fits Header object, Primary HDU, or Image HDU
      A FITS image file name, HDU data structure, or header with
      defined WCS

    bib : SBPy Bibliography instance, optional, default None
        Bibliography instance that will be populated

    Returns
    -------
    Astropy Table

    Examples
    --------
    >>> from sbpy.data import sb_search
    >>> objects = sb_search('ceres')  # doctest: +SKIP

    not yet implemented

    """


def image_search(targetid, bib=None):
    """Use the Solar System Object Image Search function of the Canadian
    Astronomy Data Centre
    (http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/ssois/) to identify
    images with a specific small body in them.

    Parameters
    ----------
    targetid : str, mandatory
        target identifier
    bib : SBPy Bibliography instance, optional, default None
        Bibliography instance that will be populated

    Returns
    -------
    Astropy Table

    Examples
    --------
    >>> from sbpy.data import Misc  # doctest: +SKIP
    >>> images = Misc.image_search('ceres')  # doctest: +SKIP

    not yet implemented

    """


def pds_ferret(targetid, bib=None):
    """Use the Small Bodies Data Ferret (http://sbntools.psi.edu/ferret/)
    at the Planetary Data System's Small Bodies Node to query for
    information on a specific small body in the PDS.

    Parameters
    ----------
    targetid : str, mandatory
        target identifier
    bib : SBPy Bibliography instance, optional, default None
        Bibliography instance that will be populated

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
