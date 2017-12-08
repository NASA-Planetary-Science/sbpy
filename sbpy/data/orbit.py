# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
SBPy data.Orbit Module
======================

Class for querying, manipulating, integrating, and fitting orbital elements.

created on June 04, 2017
"""


from .core import DataClass

__all__ = ['Orbit']

class Orbit(DataClass):
    """Class for querying, manipulating, integrating, and fitting orbital elements

    Every function of this class returns an Astropy Table object; the
    columns in these tables are not fixed and depend on the function
    generating the table or the user input.

    The `Orbit` class also provides interfaces to OpenOrb
    (https://github.com/oorb/oorb) for orbit fitting and REBOUND
    (https://github.com/hannorein/rebound) for orbit integrations.

    """



    @classmethod
    def from_horizons(cls, targetid, epoch=None, center='500@10',
                      bib=None):
        """Load orbital elements from JPL Horizons
        (https://ssd.jpl.nasa.gov/horizons.cgi).

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
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> from astropy.time import Time # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc') # doctest: +SKIP
        >>> orb = Orbit.from_horizons('Ceres', epoch) # doctest: +SKIP
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

    @classmethod
    def from_mpc(cls, targetid, bib=None):
        """Load orbital elements from the Minor Planet Center
        (http://minorplanetcenter.net/).

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
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> orb = Orbit.from_mpc('ceres') # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_astdys(cls, targetid, bib=None):
        """Load orbital elements from AstDyS
        (http://hamilton.dm.unipi.it/astdys/).

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
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> orb = Orbit.from_mpc('ceres') # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_state(cls, pos, vel):
        """Convert state vector (positions and velocities) or orbital elements.

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
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> import astropy.coordinates as coords # doctest: +SKIP
        >>> r = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=1, y=0, z=0, unit=u.au)) # doctest: +SKIP
        >>> v = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=30, y=0, z=0, unit=u.km / u.s)) # doctest: +SKIP
        >>> orb = Orbit.from_state(r, v) # doctest: +SKIP

        not yet implemented

        """

    def to_state(self, epoch):
        """Convert orbital elements to state vector (positions and velocities)

        Parameters
        ----------
        epoch : `astropy.time.Time` object, mandatory
          The epoch(s) at which to compute state vectors.
        
        Returns
        -------
        pos : `Astropy.coordinates` instance
            positions vector
        vel : `Astropy.coordinates` instance
            velocity vector

        Examples
        --------
        >>> from astropy.time import Time # doctest: +SKIP
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> orb = Orbit.from_mpc('ceres') # doctest: +SKIP
        >>> state = orb.to_state(Time('2015-03-06') # doctest: +SKIP

        not yet implemented

        """

    def orbfit(self, eph):
        """Function that fits an orbit solution to a set of ephemerides using
        the OpenOrb (https://github.com/oorb/oorb) software which has
        to be installed locally.

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
        >>> from sbpy.data import Orbit, Ephem # doctest: +SKIP
        >>> eph = Ephem.from_array([ra, dec, ra_sigma, dec_sigma, # doctest: +SKIP 
        >>>                         epochs, epochs_sigma], # doctest: +SKIP
        >>>                         names=['ra', 'dec', 'ra_sigma', # doctest: +SKIP
        >>>                                'dec_sigma', 'epochs',  # doctest: +SKIP
        >>>                                'epochs_sigma']) # doctest: +SKIP
        >>> orb = Orbit.orbfit(eph) # doctest: +SKIP

        not yet implemented

        """
        
    def integrate(self, time, integrator='IAS15'):
        """Function that integrates an orbit over a given range of time using
        the REBOUND (https://github.com/hannorein/rebound) package

        Parameters
        ----------
        time : `Astropy.units` quantity, mandatory
            Time range over which the orbit will be integrated.
        integrator : str, option, default 'IAS15'
            Integrator type to be used for the integration.

        Returns
        -------
        REBOUND simulation object

        Examples
        --------
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> orb = Orbit.from... # doctest: +SKIP
        >>> sim = orb.integrate(1000*u.year) # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_rebound(cls, sim):
        """Obtain orbital elements from REBOUND
        (https://github.com/hannorein/rebound) simulation instance.

        Parameters
        ----------
        sim : REBOUND simulation instance, mandatory
            Simulation from which to obtain orbital elements.

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> orb = Orbit.from... # doctest: +SKIP
        >>> sim = Orbit.integrate(orb, time=1000*u.year) # doctest: +SKIP
        >>> future_orb = Orbit.from_rebound(sim) # doctest: +SKIP

        not yet implemented

        """
