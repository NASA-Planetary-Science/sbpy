# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
SBPy data.Ephem Module
======================

Class for storing and querying ephemerides

created on June 04, 2017
"""


from .core import DataClass

__all__ = ['Ephem']

class Ephem(DataClass):
    """Class for storing and querying ephemerides
    
    The `Ephem` class provides an interface to PyEphem
    (http://rhodesmill.org/pyephem/) for ephemeris calculations.
    
    """

    @classmethod
    def from_horizons(cls, targetid, epoch, observatory, center='500@10',
                      bib=None):
        """Load orbital elements from JPL Horizons (https://ssd.jpl.nasa.gov/horizons.cgi).

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

        """

        from astropy.time import Time
        
        if epoch is None:
            epoch = [Time.now()]
        elif isinstance(epoch, Time):
            epoch = [epoch]
            
        # for now, use CALLHORIZONS for the query; this will be replaced with
        # a dedicated query
        import callhorizons
        el = callhorizons.query(targetid)
        el.set_discreteepochs([ep.jd for ep in epoch])
        el.get_ephemerides(observatory)

        data = [el[field] for field in el.fields]
        names = el.fields
        #meta = {'name': 'orbital elements from JPL Horizons'}
        # table = Table([el[field] for field in el.fields],
        #               names=el.fields,
        #               meta={'name': 'orbital elements from JPL Horizons'})
        # # Astropy units will be integrated in the future

        from .. import bib        
        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem', {'implementation':
                                             '1996DPS....28.2504G'})
        
        return cls.from_array(data, names)

    @classmethod
    def from_mpc(cls, targetid, epoch, observatory='500', bib=None):
        """Load ephemerides from the Minor Planet Center (http://minorplanetcenter.net/).

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
        >>> from sbpy.data import Ephem # doctest: +SKIP
        >>> from astropy.time import Time # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc') # doctest: +SKIP
        >>> eph = Ephem.from_mpc('ceres', '568', epoch) # doctest: +SKIP

        not yet implemented

        """

    def report_to_mpc(bib=None):
        """Format as a report to the Minor Planet Center
        (http://minorplanetcenter.net/).

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
        >>> from sbpy.data import Ephem # doctest: +SKIP
        >>> eph = Ephem.from_array... # doctest: +SKIP
        >>> report = eph.report_to_mpc() # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_imcce(cls, targetid, epoch, observatory='500', bib=None):
        """Load orbital elements from IMCCE (http://vo.imcce.fr/webservices/miriade/).
           
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
        >>> from sbpy.data import Ephem # doctest: +SKIP
        >>> from astropy.time import Time # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc') # doctest: +SKIP
        >>> eph = Ephem.from_imcce('ceres', '568', epoch) # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_lowell(cls, targetid, epoch, observatory='500', bib=None):
        """Load orbital elements from Lowell Observatory (http://asteroid.lowell.edu/).

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
        >>> from sbpy.data import Ephem # doctest: +SKIP
        >>> from astropy.time import Time # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc') # doctest: +SKIP
        >>> eph = Ephem.from_lowell('ceres', '568', epoch) # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_pyephem(cls, orb, location, epoch):
        """Function that derives ephemerides based on an `Astropy.table`
        containing orbital elements using PyEphem (http://rhodesmill.org/pyephem/).
        
        Parameters
        ----------
        orb : `Astropy.table`, mandatory
            complete set of orbital elements
        location : str or dictionary, mandatory
            name of location or a dictionary fully describing the location
        epoch : `Astropy.time` object
            
        Examples
        --------
        >>> from sbpy.data import Ephem, Orbit # doctest: +SKIP
        >>> orb = Orbit.from_... # doctest: +SKIP
        >>> eph = Ephem.from_pyephem(orb,  # doctest: +SKIP
        >>>                          location={'name':'Flagstaff', # doctest: +SKIP
        >>>                                    'geolon':35.199167, # doctest: +SKIP
        >>>                                    'geolat':-111.631111, # doctest: +SKIP
        >>>                                    'altitude':'2106'},# doctest: +SKIP
        >>>                          epoch=epoch) # doctest: +SKIP

        not yet implemented

        """        
