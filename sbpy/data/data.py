# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
================
SBPy Data Module
================

created on June 22, 2017
"""

from astropy.table import Table, Column
from astropy.time import Time
import astropy.units as u
import callhorizons

__all__ = ['Orbit', 'Ephem', 'Phys', 'Misc', 'DataClass']


class DataClass():

    def __init__(self, *args, **kwargs):
        """Create Astropy.table from kwargs"""
        self.table = Table()
        for key, val in kwargs.items():
            try:
                unit = val.unit
                val = val.value
            except:
                unit = None
            # check if val is already list-like
            try:
                val[0]
            except TypeError:
                val = [val]

            self.table[key] = Column(val, unit=unit)

    @classmethod
    def from_dict(cls, data):
        """Create orbital elements table from dictionary or list of
        dictionaries
        
        Parameters
        ----------
        data : dictionary or list of dicts, mandatory
            data that will be rearranged in Astropy Table format

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> import astropy.units as u
        >>> orb = Orbit.from_dict({'a': 2.7674*u.au, 
        >>>                        'e': .0756,
        >>>                        'i': 10.59321*u.deg})

        """
        return cls(**data)
    
    @classmethod
    def from_array(cls, data, names):
        """Create orbital elements table from lists or arrays

        Parameters
        ----------
        data : list or array, mandatory
            data that will be rearraned in Astropy Table format, one array per 
            column
        names : list, mandatory
            column names, must have n names for n `data` arrays

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> import astropy.units as u
        >>> from numpy.random import random as r
        >>> orb = Orbit.from_array(data=[r(100)*2*u.au,
        >>>                              r(100),
        >>>                              r(100)*180*u.deg],
        >>>                        names=['a', 'e', 'i'])

        """

        return cls.from_dict(dict(zip(names, data)))
    
    def __getattr__(self, field):
        if field in self.table.columns:
            return self.table[field]
        else:
           raise AttributeError ("field '{:s}' does not exist".format(field))

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
        return self.table
    
            
class Orbit(DataClass):
    """Class for querying, manipulating, integrating, and fitting orbital elements

    Every function of this class returns an Astropy Table object; the
    columns in these tables are not fixed and depend on the function
    generating the table or the user input.

    The `Orbit` class also provides interfaces to `OpenOrb`_ for orbit
    fitting and `REBOUND`_ for orbit integrations.

    .. _OpenOrb: https://github.com/oorb/oorb
    .. _REBOUND: https://github.com/hannorein/rebound

    """



    @classmethod
    def from_horizons(cls, targetid, epoch=None, center='500@10',
                      bib=None):
        """Load orbital elements from `JPL Horizons`_

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epoch : astropy Time instance or iterable, optional, default None
            epoch of elements; if None is provided, current date is used
        center : str, optional, default '500@10' (Sun)
            center body of orbital elements
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        preliminary implementation
        
        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> orb = Orbit.from_horizons('ceres', epoch)

        .. _JPL Horizons: https://ssd.jpl.nasa.gov/horizons.cgi

        """

        if epoch is None:
            epoch = [Time.now()]
        elif isinstance(epoch, Time):
            epoch = [epoch]

        # for now, use CALLHORIZONS for the query; this will be replaced with
        # a dedicated query
        el = callhorizons.query(targetid)
        el.set_discreteepochs([ep.jd for ep in epoch])
        el.get_elements(center=center)
        data = [el[field] for field in el.fields]
        names = el.fields
        #meta = {'name': 'orbital elements from JPL Horizons'}
        # table = Table([el[field] for field in el.fields],
        #               names=el.fields,
        #               meta={'name': 'orbital elements from JPL Horizons'})
        # # Astropy units will be integrated in the future

        if bib is not None:
            bib['Horizons orbital elements query'] = {'implementation':
                                                      '1996DPS....28.2504G'}
            
        return cls.from_array(data, names)


    def from_mpc(self, targetid, bib=None):
        """Load orbital elements from the `Minor Planet Center`_

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
        >>> from sbpy.data import Orbit
        >>> orb = Orbit.from_mpc('ceres')

        not yet implemented

        .. _Minor Planet Center: http://minorplanetcenter.net/

        """

    def from_astdys(self, targetid, bib=None):
        """Load orbital elements from `AstDyS`_

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
        >>> from sbpy.data import Orbit
        >>> orb = Orbit.from_mpc('ceres')

        not yet implemented

        .. _AstDyS: http://hamilton.dm.unipi.it/astdys/ 

        """

    def from_state(self, pos, vel):
        """Convert state vector (positions and velocities) or orbital elements

        Parameters
        ----------
        pos : `Astropy.coordinates` instance, mandatory
            positions vector
        vel : `Astropy.coordinates` instance, mandatory
            velocity vector
        
        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> import astropy.coordinates as coords
        >>> r = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=1, y=0, z=0, unit=u.au))
        >>> v = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=30, y=0, z=0, unit=u.km / u.s))
        >>> orb = Orbit.from_state(r, v)

        not yet implemented

        """

    def to_state(self, pos, vel):
        """Convert orbital elements to state vector (positions and velocities)

        Parameters
        ----------
        obs : `Astropy.table` instance, mandatory
            orbital elements
        
        Returns
        -------
        pos : `Astropy.coordinates` instance
            positions vector
        vel : `Astropy.coordinates` instance
            velocity vector

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> orb = Orbit.from_mpc('ceres')        
        >>> state = Orbit.to_state(orb)      

        not yet implemented

        """

    def orbfit(self, eph):
        """Function that fits an orbit solution to a set of ephemerides using
        the `OpenOrb`_ software which has to be installed locally.

        Parameters
        ----------
        eph : `Astropy.table`, mandatory
            set of ephemerides with mandatory columns `ra`, `dec`, `epoch` and 
            optional columns `ra_sig`, `dec_sig`, `epoch_sig` 
        
        additional parameters will be identified in the future

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Orbit, Ephem
        >>> eph = Ephem.from_array([ra, dec, ra_sigma, dec_sigma, 
        >>>                         epochs, epochs_sigma],
        >>>                        names=['ra', 'dec', 'ra_sigma', 
        >>>                               'dec_sigma', 'epochs', 
        >>>                               'epochs_sigma'])
        >>> orb = Orbit.orbfit(eph)

        not yet implemented

        .. _OpenOrb: https://github.com/oorb/oorb

        """
        
    def integrate(self, orb, time, integrator='IAS15'):
        """Function that integrates an orbit over a given range of time using the `REBOUND`_ package

        Parameters
        ----------
        orb : `Astropy.table`, mandatory
            complete set of orbital elements
        time : `Astropy.units` quantity, mandatory
            time range over which the orbit will be integrated 
        integrator : str, option, default 'IAS15'
            integrator type to be used for the integration

        Returns
        -------
        REBOUND simulation object

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> orb = Orbit.from...
        >>> sim = Orbit.integrate(orb, time=1000*u.year)

        not yet implemented

        .. _REBOUND: https://github.com/hannorein/rebound
        """

    def from_rebound(self, sim):
        """Obtain orbital elements from `REBOUND`_ simulation instance

        Parameters
        ----------
        sim : REBOUND simulation instance, mandatory
            simulation from which to obtain orbital elements

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> orb = Orbit.from...
        >>> sim = Orbit.integrate(orb, time=1000*u.year)
        >>> future_orb = Orbit.from_rebound(sim)

        not yet implemented

        .. _REBOUND: https://github.com/hannorein/rebound

        """
        
class Ephem(DataClass):
    """Class for storing and querying ephemerides
    
    The `Ephem` class provides an interface to `PyEphem`_ for
    ephemeris calculations.
    
    .. _PyEphem: http://rhodesmill.org/pyephem/
    
    """

    @classmethod
    def from_horizons(cls, targetid, epoch, observatory, center='500@10',
                      bib=None):
        """Load orbital elements from `JPL Horizons`_

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epoch : astropy Time instance or iterable, optional, default None
            epoch of elements; if None is provided, current date is used
        center : str, optional, default '500@10' (Sun)
            center body of orbital elements
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        preliminary implementation
        
        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> orb = Orbit.from_horizons('ceres', epoch)

        .. _JPL Horizons: https://ssd.jpl.nasa.gov/horizons.cgi

        """

        try:
            dummy = epoch[0]
        except TypeError:
            epoch = [epoch]

            
        # for now, use CALLHORIZONS for the query; this will be replaced with
        # a dedicated query
        el = callhorizons.query(targetid)
        el.set_discreteepochs([ep.jd for ep in epoch])
        el.get_ephemerides(observatory)

        print(el.dates)

        data = [el[field] for field in el.fields]
        names = el.fields
        #meta = {'name': 'orbital elements from JPL Horizons'}
        # table = Table([el[field] for field in el.fields],
        #               names=el.fields,
        #               meta={'name': 'orbital elements from JPL Horizons'})
        # # Astropy units will be integrated in the future

        if bib is not None:
            bib['Horizons orbital elements query'] = {'implementation':
                                                      '1996DPS....28.2504G'}
            
        return cls.from_array(data, names)

    @classmethod
    def from_mpc(cls, targetid, epoch, observatory='500', bib=None):
        """Load ephemerides from the `Minor Planet Center`_

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default None
            epoch of elements; if None is provided, current date is used
        observatory : str, optional, default '500' (geocentric)
            location of observer
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_mpc('ceres', '568', epoch)

        not yet implemented

        .. _Minor Planet Center: http://minorplanetcenter.net/

        """

    @classmethod
    def report_to_mpc(cls, eph, bib=None):
        """Format ephemerides `Astropy.table` to report to `Minor Planet Center`_ 

        Parameters
        ----------
        eph : `Astropy.table` of ephemerides
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated
        
        additional parameters will be identified in the future
        
        Returns
        -------
        str

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> eph = Ephem.from_array...
        >>> report = Ephem.report_to_mpc(eph)

        not yet implemented

        .. _Minor Planet Center: http://minorplanetcenter.net/

        """

    @classmethod
    def from_imcce(cls, targetid, epoch, observatory='500', bib=None):
        """Load orbital elements from `IMCCE`_
           
        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default None
            epoch of elements; if None is provided, current date is used
        observatory : str, optional, default '500' (geocentric)
            location of observer
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_imcce('ceres', '568', epoch)

        not yet implemented

        .. _IMCCE: http://vo.imcce.fr/webservices/miriade/

        """

    @classmethod
    def from_lowell(cls, targetid, epoch, observatory='500', bib=None):
        """Load orbital elements from `Lowell Observatory`_

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default None
            epoch of elements; if None is provided, current date is used
        observatory : str, optional, default '500' (geocentric)
            location of observer
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_lowell('ceres', '568', epoch)

        not yet implemented

        .. _Lowell Observatory: http://asteroid.lowell.edu/ 

        """

    @classmethod
    def from_pyephem(cls, orb, location, epoch):
        """Function that derives ephemerides based on an `Astropy.table`
        containing orbital elements using `PyEphem`_.
        
        Parameters
        ----------
        orb : `Astropy.table`, mandatory
            complete set of orbital elements
        location : str or dictionary, mandatory
            name of location or a dictionary fully describing the location
        epoch : `Astropy.time` object
            
        Examples
        --------
        >>> from sbpy.data import Ephem, Orbit
        >>> orb = Orbit.from_...
        >>> eph = Ephem.from_pyephem(orb, 
        >>>                          location={'name':'Flagstaff', 
        >>>                                    'geolon':35.199167, 
        >>>                                    'geolat':-111.631111, 
        >>>                                    'altitude':'2106'},
        >>>                          epoch=epoch)

        not yet implemented

        .. _PyEphem: http://rhodesmill.org/pyephem/

        """        

class Phys():
    """Class for storing and querying physical properties"""

    @classmethod
    def from_horizons(cls, targetid, bib=None):
        """Load physical properties from `JPL Horizons`_

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
        >>> from sbpy.data import Phys
        >>> phys = Phys.from_horizons('ceres'(

        not yet implemented

        .. _JPL Horizons: https://ssd.jpl.nasa.gov/horizons.cgi

        """

    @classmethod
    def from_lowell(cls, targetid, bib=None):
        """Load physical properties from `Lowell Observatory`_ 

        The Lowell database will provide a database of physical
        properties which is a compilation of a number of different sources.

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
        >>> from sbpy.data import Phys
        >>> phys = Phys.from_astorb('ceres'(

        not yet implemented

        .. _Lowell Observatory: http://asteroid.lowell.edu/ 

        """

    def derive_absmag(self):
        """Derive absolute magnitude from diameter and geometric albedo"""

    def derive_diam(self):
        """Derive diameter from absolute magnitude and geometric albedo"""
        
    def derive_pv(self):
        """Derive geometric albedo from diameter and absolute magnitude"""

    def derive_bondalbedo(self):
        """Derive Bond albedo from geometric albedo and photometric phase slope"""

    
        

class Misc():
    """Class for obtaining miscellaneous data on small bodies"""

    def mpc_observations(targetid, bib=None):
        """Function that obtains all available observations of a small body from the `Minor Planet Center`_ and provides them in the form of a `sbpy.data.Ephem` Astropy table

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
        >>> from sbpy.data import Misc
        >>> eph = Misc.mpc_observations('ceres')

        not yet implemented

        .. _Minor Planet Center: http://www.minorplanetcenter.net

        """

    def sb_search(filename, bib=None):
        """Function that uses the `Skybot`_ service at IMCCE to identify moving objects potentially present in a registered FITS images 

        Parameters
        ----------
        filename : str, mandatory
            filename of FITS image
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Misc
        >>> eph = Misc.sb_search('ceres')

        not yet implemented

        .. _Skybot: http://vo.imcce.fr/webservices/skybot/

        """
        
    def image_search(targetid, bib=None):
        """Function that uses the Solar System Object Image Search function of the `Canadian Astronomy Data Centre`_ to identify images with a specific small body in them

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
        >>> from sbpy.data import Misc
        >>> eph = Misc.image_search('ceres')

        not yet implemented

        .. _Canadian Astronomy Data Centre: http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/ssois/

        """
        
    def pds_ferret(targetid, bib=None):
        """Function that uses the `Small Bodies Data Ferret`_ at the Planetary Data System's Small Bodies Node to query for all existing information on a specific small body in the PDS

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
        >>> from sbpy.data import Misc
        >>> eph = Misc.pds_ferret('ceres')

        not yet implemented

        .. _Small Bodies Data Ferret: http://sbntools.psi.edu/ferret/

        """

    def asteroid_or_comet(targetid):
        """Checks if an object is an asteroid, or a comet, based on its targetid

        Examples
        --------
        >>> from sbpy.data import Misc
        >>> print(Misc.asteroid_or_comet('2P')
            'comet'
        
        note yet implemented

        """
        
    def altident(targetid, bib=None):
        """Query Lowell database to obtain alternative target names for `targetitd`

        Examples
        --------
        >>> from sbpy.data import Misc
        >>> print(Misc.altident('3552'))
            ['3552', 'Don Quixote', '1983 SA']

        not yet implemented

        """

    def from_packed(targetid):
        """Convert packed designation/number to unpacked"""

    def to_packed(targetid):
        """ Convert unpacked identifier to packed designation/number"""
        
