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
    `dict` (DataClass.from_dict) or `array`-like
    (DataClass.from_array) objects. A few high-level functions for
    table modification are provided; other modifications can be
    applied to the table object (`DataClass.table`) directly.
    """

    def __init__(self, **kwargs):
        """Build data table from `**kwargs`."""
        self.table = Table()
        # self.altkeys = {}  # dictionary for alternative column names
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
        tab : `~DataClass` object

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
             of the sub-list or sub-array on the 0-axis is interpreted becomesdictionary creates a row in the data table. Dictionary
             keys are used as column names. If a list of dicitionaries
             is provided, all dictionaries have to provide them same
             set of keys (and units, if used at all).
        data : list or array, mandatory
            data that will be rearranged in Astropy Table format, one array 
            per column
        names : list, mandatory
            column names, must have n names for n `data` arrays

        Returns
        -------
        Astropy Table

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


# class Orbit(DataClass):
#

        """Class for querying, manipulating, integrating, and fitting orbital elements

#     Every function of this class returns an Astropy Table object; the
#     columns in these tables are not fixed and depend on the function
#     generating the table or the user input.

#     The `Orbit` class also provides interfaces to OpenOrb
#     (https://github.com/oorb/oorb) for orbit fitting and REBOUND
#     (https://github.com/hannorein/rebound) for orbit integrations.

#     """


#     @classmethod
#     def from_horizons(cls, targetid, epoch=None, center='500@10',
#                       bib=None):
#         """Load orbital elements from JPL Horizons
#         (https://ssd.jpl.nasa.gov/horizons.cgi).

#         Parameters
#         ----------
#         targetid : str, mandatory
#             target identifier
#         epoch : astropy Time instance or iterable, optional, default None
#             epoch of elements; if None is provided, current date is used
#         center : str, optional, default '500@10' (Sun)
#             center body of orbital elements
#         bib : SBPy Bibliography instance, optional, default None
#             Bibliography instance that will be populated

#         preliminary implementation

#         Returns
#         -------
#         Astropy Table

#         Examples
#         --------
#         >>> from sbpy.data import Orbit
#         >>> from astropy.time import Time
#         >>> epoch = Time('2018-05-14', scale='utc')
#         >>> orb = Orbit.from_horizons('ceres', epoch)
#         """

#         if epoch is None:
#             epoch = [Time.now()]
#         elif isinstance(epoch, Time):
#             epoch = [epoch]

#         # for now, use CALLHORIZONS for the query; this will be replaced with
#         # a dedicated query
#         el = callhorizons.query(targetid)
#         el.set_discreteepochs([ep.jd for ep in epoch])
#         el.get_elements(center=center)
#         data = [el[field] for field in el.fields]
#         names = el.fields
#         #meta = {'name': 'orbital elements from JPL Horizons'}
#         # table = Table([el[field] for field in el.fields],
#         #               names=el.fields,
#         #               meta={'name': 'orbital elements from JPL Horizons'})
#         # # Astropy units will be integrated in the future

#         if bib is not None:
#             bib['Horizons orbital elements query'] = {'implementation':
#                                                       '1996DPS....28.2504G'}

#         return cls.from_array(data, names)

#     @classmethod
#     def from_mpc(cls, targetid, bib=None):
#         """Load orbital elements from the Minor Planet Center
#         (http://minorplanetcenter.net/).

#         Parameters
#         ----------
#         targetid : str, mandatory
#             target identifier
#         bib : SBPy Bibliography instance, optional, default None
#             Bibliography instance that will be populated

#         Returns
#         -------
#         Astropy Table

#         Examples
#         --------
#         >>> from sbpy.data import Orbit
#         >>> orb = Orbit.from_mpc('ceres')

#         not yet implemented

#         """

#     @classmethod
#     def from_astdys(cls, targetid, bib=None):
#         """Load orbital elements from AstDyS
#         (http://hamilton.dm.unipi.it/astdys/).

#         Parameters
#         ----------
#         targetid : str, mandatory
#             target identifier
#         bib : SBPy Bibliography instance, optional, default None
#             Bibliography instance that will be populated

#         Returns
#         -------
#         Astropy Table

#         Examples
#         --------
#         >>> from sbpy.data import Orbit
#         >>> orb = Orbit.from_mpc('ceres')

#         not yet implemented

#         """

#     @classmethod
#     def from_state(cls, pos, vel):
#         """Convert state vector (positions and velocities) or orbital elements.

#         Parameters
#         ----------
#         pos : `Astropy.coordinates` instance, mandatory
#             positions vector
#         vel : `Astropy.coordinates` instance, mandatory
#             velocity vector

#         Returns
#         -------
#         Astropy Table

#         Examples
#         --------
#         >>> from sbpy.data import Orbit
#         >>> import astropy.coordinates as coords
#         >>> r = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=1, y=0, z=0, unit=u.au))
#         >>> v = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=30, y=0, z=0, unit=u.km / u.s))
#         >>> orb = Orbit.from_state(r, v)

#         not yet implemented

#         """

#     def to_state(self, epoch):
#         """Convert orbital elements to state vector (positions and velocities)

#         Parameters
#         ----------
#         epoch : `astropy.time.Time` object, mandatory
#           The epoch(s) at which to compute state vectors.

#         Returns
#         -------
#         pos : `Astropy.coordinates` instance
#             positions vector
#         vel : `Astropy.coordinates` instance
#             velocity vector

#         Examples
#         --------
#         >>> from astropy.time import Time
#         >>> from sbpy.data import Orbit
#         >>> orb = Orbit.from_mpc('ceres')
#         >>> state = orb.to_state(Time('2015-03-06')

#         not yet implemented

#         """

#     def orbfit(self, eph):
#         """Function that fits an orbit solution to a set of ephemerides using
#         the OpenOrb (https://github.com/oorb/oorb) software which has
#         to be installed locally.

#         Parameters
#         ----------
#         eph : `Astropy.table`, mandatory
#             set of ephemerides with mandatory columns `ra`, `dec`, `epoch` and
#             optional columns `ra_sig`, `dec_sig`, `epoch_sig`

#         additional parameters will be identified in the future

#         Returns
#         -------
#         Astropy Table

#         Examples
#         --------
#         >>> from sbpy.data import Orbit, Ephem
#         >>> eph = Ephem.from_array([ra, dec, ra_sigma, dec_sigma,
#         >>>                         epochs, epochs_sigma],
#         >>>                         names=['ra', 'dec', 'ra_sigma',
#         >>>                                'dec_sigma', 'epochs',
#         >>>                                'epochs_sigma'])
#         >>> orb = Orbit.orbfit(eph)

#         not yet implemented

#         """

#     def integrate(self, time, integrator='IAS15'):
#         """Function that integrates an orbit over a given range of time using
#         the REBOUND (https://github.com/hannorein/rebound) package

#         Parameters
#         ----------
#         time : `Astropy.units` quantity, mandatory
#             Time range over which the orbit will be integrated.
#         integrator : str, option, default 'IAS15'
#             Integrator type to be used for the integration.

#         Returns
#         -------
#         REBOUND simulation object

#         Examples
#         --------
#         >>> from sbpy.data import Orbit
#         >>> orb = Orbit.from...
#         >>> sim = orb.integrate(1000*u.year)

#         not yet implemented

#         """

#     @classmethod
#     def from_rebound(cls, sim):
#         """Obtain orbital elements from REBOUND
#         (https://github.com/hannorein/rebound) simulation instance.

#         Parameters
#         ----------
#         sim : REBOUND simulation instance, mandatory
#             Simulation from which to obtain orbital elements.

#         Returns
#         -------
#         Astropy Table

#         Examples
#         --------
#         >>> from sbpy.data import Orbit
#         >>> orb = Orbit.from...
#         >>> sim = Orbit.integrate(orb, time=1000*u.year)
#         >>> future_orb = Orbit.from_rebound(sim)

#         not yet implemented

#         """

# # class Ephem(DataClass):
# #     """Class for storing and querying ephemerides

# #     The `Ephem` class provides an interface to PyEphem
# #     (http://rhodesmill.org/pyephem/) for ephemeris calculations.

# #     """

# #     @classmethod
# #     def from_horizons(cls, targetid, epoch, observatory, center='500@10',
# #                       bib=None):
# #         """Load orbital elements from JPL Horizons (https://ssd.jpl.nasa.gov/horizons.cgi).

# #         Parameters
# #         ----------
# #         targetid : str, mandatory
# #             Target identifier.
# #         epoch : astropy Time instance or iterable, optional, default None
# #             Epoch of elements; if None is provided, current date is used.
# #         center : str, optional, default '500@10' (Sun)
# #             Center body of orbital elements.
# #         bib : SBPy Bibliography instance, optional, default None
# #             Bibliography instance that will be populated.

# #         preliminary implementation

# #         Returns
# #         -------
# #         Astropy Table

# #         Examples
# #         --------
# #         >>> from sbpy.data import Orbit
# #         >>> from astropy.time import Time
# #         >>> epoch = Time('2018-05-14', scale='utc')
# #         >>> orb = Orbit.from_horizons('ceres', epoch)

# #         """

# #         try:
# #             dummy = epoch[0]
# #         except TypeError:
# #             epoch = [epoch]


# #         # for now, use CALLHORIZONS for the query; this will be replaced with
# #         # a dedicated query
# #         el = callhorizons.query(targetid)
# #         el.set_discreteepochs([ep.jd for ep in epoch])
# #         el.get_ephemerides(observatory)

# #         print(el.dates)

# #         data = [el[field] for field in el.fields]
# #         names = el.fields
# #         #meta = {'name': 'orbital elements from JPL Horizons'}
# #         # table = Table([el[field] for field in el.fields],
# #         #               names=el.fields,
# #         #               meta={'name': 'orbital elements from JPL Horizons'})
# #         # # Astropy units will be integrated in the future

# #         if bib is not None:
# #             bib['Horizons orbital elements query'] = {'implementation':
# #                                                       '1996DPS....28.2504G'}

# #         return cls.from_array(data, names)

# #     @classmethod
# #     def from_mpc(cls, targetid, epoch, observatory='500', bib=None):
# #         """Load ephemerides from the Minor Planet Center (http://minorplanetcenter.net/).

# #         Parameters
# #         ----------
# #         targetid : str, mandatory
# #             target identifier
# #         epochs : astropy Time instance or iterable, optional, default None
# #             epoch of elements; if None is provided, current date is used
# #         observatory : str, optional, default '500' (geocentric)
# #             location of observer
# #         bib : SBPy Bibliography instance, optional, default None
# #             Bibliography instance that will be populated

# #         Returns
# #         -------
# #         Astropy Table

# #         Examples
# #         --------
# #         >>> from sbpy.data import Ephem
# #         >>> from astropy.time import Time
# #         >>> epoch = Time('2018-05-14', scale='utc')
# #         >>> eph = Ephem.from_mpc('ceres', '568', epoch)

# #         not yet implemented

# #         """

# #     def report_to_mpc(bib=None):
# #         """Format as a report to the Minor Planet Center
# #         (http://minorplanetcenter.net/).

# #         Parameters
# #         ----------
# #         bib : SBPy Bibliography instance, optional, default None
# #             Bibliography instance that will be populated

# #         additional parameters will be identified in the future

# #         Returns
# #         -------
# #         str

# #         Examples
# #         --------
# #         >>> from sbpy.data import Ephem
# #         >>> eph = Ephem.from_array...
# #         >>> report = eph.report_to_mpc()

# #         not yet implemented

# #         """

# #     @classmethod
# #     def from_imcce(cls, targetid, epoch, observatory='500', bib=None):
# #         """Load orbital elements from IMCCE (http://vo.imcce.fr/webservices/miriade/).

# #         Parameters
# #         ----------
# #         targetid : str, mandatory
# #             target identifier
# #         epochs : astropy Time instance or iterable, optional, default None
# #             epoch of elements; if None is provided, current date is used
# #         observatory : str, optional, default '500' (geocentric)
# #             location of observer
# #         bib : SBPy Bibliography instance, optional, default None
# #             Bibliography instance that will be populated

# #         Returns
# #         -------
# #         Astropy Table

# #         Examples
# #         --------
# #         >>> from sbpy.data import Ephem
# #         >>> from astropy.time import Time
# #         >>> epoch = Time('2018-05-14', scale='utc')
# #         >>> eph = Ephem.from_imcce('ceres', '568', epoch)

# #         not yet implemented

# #         """

# #     @classmethod
# #     def from_lowell(cls, targetid, epoch, observatory='500', bib=None):
# #         """Load orbital elements from Lowell Observatory (http://asteroid.lowell.edu/).

# #         Parameters
# #         ----------
# #         targetid : str, mandatory
# #             target identifier
# #         epochs : astropy Time instance or iterable, optional, default None
# #             epoch of elements; if None is provided, current date is used
# #         observatory : str, optional, default '500' (geocentric)
# #             location of observer
# #         bib : SBPy Bibliography instance, optional, default None
# #             Bibliography instance that will be populated

# #         Returns
# #         -------
# #         Astropy Table

# #         Examples
# #         --------
# #         >>> from sbpy.data import Ephem
# #         >>> from astropy.time import Time
# #         >>> epoch = Time('2018-05-14', scale='utc')
# #         >>> eph = Ephem.from_lowell('ceres', '568', epoch)

# #         not yet implemented

# #         """

# #     @classmethod
# #     def from_pyephem(cls, orb, location, epoch):
# #         """Function that derives ephemerides based on an `Astropy.table`
# #         containing orbital elements using PyEphem (http://rhodesmill.org/pyephem/).

# #         Parameters
# #         ----------
# #         orb : `Astropy.table`, mandatory
# #             complete set of orbital elements
# #         location : str or dictionary, mandatory
# #             name of location or a dictionary fully describing the location
# #         epoch : `Astropy.time` object

# #         Examples
# #         --------
# #         >>> from sbpy.data import Ephem, Orbit
# #         >>> orb = Orbit.from_...
# #         >>> eph = Ephem.from_pyephem(orb,
# #         >>>                          location={'name':'Flagstaff',
# #         >>>                                    'geolon':35.199167,
# #         >>>                                    'geolat':-111.631111,
# #         >>>                                    'altitude':'2106'},
# #         >>>                          epoch=epoch)

# #         not yet implemented

# #         """

# # class Phys(DataClass):
# #     """Class for storing and querying physical properties"""

# #     @classmethod
# #     def from_horizons(cls, targetid, bib=None):
# #         """Load physical properties from JPL Horizons
# #         (https://ssd.jpl.nasa.gov/horizons.cgi)

# #         Parameters
# #         ----------
# #         targetid : str, mandatory
# #             target identifier
# #         bib : SBPy Bibliography instance, optional, default None
# #             Bibliography instance that will be populated

# #         Returns
# #         -------
# #         Astropy Table

# #         Examples
# #         --------
# #         >>> from sbpy.data import Phys
# #         >>> phys = Phys.from_horizons('ceres'(

# #         not yet implemented

# #         """

# #     @classmethod
# #     def from_lowell(cls, targetid, bib=None):
# #         """Load physical properties from Lowell Observatory
# #         (http://asteroid.lowell.edu/).

# #         The Lowell database will provide a database of physical
# #         properties which is a compilation of a number of different sources.

# #         Parameters
# #         ----------
# #         targetid : str, mandatory
# #             target identifier
# #         bib : SBPy Bibliography instance, optional, default None
# #             Bibliography instance that will be populated

# #         Returns
# #         -------
# #         Astropy Table

# #         Examples
# #         --------
# #         >>> from sbpy.data import Phys
# #         >>> phys = Phys.from_astorb('ceres'(

# #         not yet implemented

# #         """

# #     def derive_absmag(self):
# #         """Derive absolute magnitude from diameter and geometric albedo"""

# #     def derive_diam(self):
# #         """Derive diameter from absolute magnitude and geometric albedo"""

# #     def derive_pv(self):
# #         """Derive geometric albedo from diameter and absolute magnitude"""

# #     def derive_bondalbedo(self):
# #         """Derive Bond albedo from geometric albedo and photometric phase slope"""


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
