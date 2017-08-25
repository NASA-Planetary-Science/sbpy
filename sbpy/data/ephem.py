# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=====================
SBPy data.Ephem Class
=====================

created on June 22, 2017
"""

__all__ = ['Ephem']

from astropy.table import Table, Column
from astropy.time import Time
import astropy.units as u
import callhorizons
from .data import DataClass

class Ephem(DataClass):
    """Class for storing and querying ephemerides
    
    The `Ephem` class provides an interface to `PyEphem`_ for
    ephemeris calculations.
    
    .. _PyEphem: http://rhodesmill.org/pyephem/
    
    """

    @classmethod
    def from_horizons(cls, targetid, epoch, observatory, center='500@10',
                      bib=None):
        """Load orbital elements from `JPL Horizons`_.

        Parameters
        ----------
        targetid : str, mandatory
            Target identifier.
        epoch : astropy Time instance or iterable, optional, default None
            Epoch of elements; if None is provided, current date is used.
        center : str, optional, default '500@10' (Sun)
            Center body of orbital elements.
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated.

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

    def report_to_mpc(bib=None):
        """Format as a report to the `Minor Planet Center`_.

        Parameters
        ----------
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
        >>> report = eph.report_to_mpc()

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

