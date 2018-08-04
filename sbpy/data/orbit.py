# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Orbit Module
======================

Class for querying, manipulating, integrating, and fitting orbital elements.

created on June 04, 2017
"""
from numpy import ndarray
from astropy.time import Time
from astropy.table import vstack
from astroquery.jplhorizons import Horizons

from .. import bib
from .core import DataClass

__all__ = ['Orbit']


class Orbit(DataClass):
    """Class for querying, manipulating, integrating, and fitting orbital
    elements"""

    @classmethod
    def from_horizons(cls, targetids, id_type='smallbody',
                      epochs=None, center='500@10',
                      **kwargs):
        """Load target orbital elements from
        `JPL Horizons <https://ssd.jpl.nasa.gov/horizons.cgi>`_ using
        `astroquery.jplhorizons.HorizonsClass.elements`

        Parameters
        ----------
        targetids : str or iterable of str
            Target identifier, i.e., a number, name, designation, or JPL
            Horizons record number, for one or more targets.
        id_type : str, optional
            The nature of ``targetids`` provided; possible values are
            ``'smallbody'`` (asteroid or comet), ``'majorbody'`` (planet or
            satellite), ``'designation'`` (asteroid or comet designation),
            ``'name'`` (asteroid or comet name), ``'asteroid_name'``,
            ``'comet_name'``, ``'id'`` (Horizons id).
            Default: ``'smallbody'``
        epochs : `~astropy.time.Time` object or iterable thereof, or dictionary, optional
            Epochs of elements to be queried; a list, tuple or
            `~numpy.ndarray` of `~astropy.time.Time` objects or Julian
            Dates as floats should be used for a number of discrete
            epochs; a dictionary including keywords ``start``,
            ``step``, and ``stop`` can be used to generate a range of
            epochs (see
            `~astroquery.jplhorizons.HorizonsClass.Horizons.elements`
            for details); if ``None`` is provided, current date and
            time are used. Default: ``None``
        center : str, optional, default ``'500@10'`` (center of the Sun)
            Elements will be provided relative to this position.
        **kwargs : optional
            Arguments that will be provided to
            `astroquery.jplhorizons.HorizonsClass.elements`.

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Orbit.from_horizons('Ceres', epochs=epoch)
        """

        # modify epoch input to make it work with astroquery.jplhorizons
        # maybe this stuff should really go into that module....
        if epochs is None:
            epochs = [Time.now().jd]
        elif isinstance(epochs, Time):
            epochs = [Time(epochs).jd]
        elif isinstance(epochs, dict):
            for key, val in epochs.items():
                if isinstance(val, Time):
                    epochs[key] = str(val.utc)

       # if targetids is a list, run separate Horizons queries and append
        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        # append elements table for each targetid
        all_elem = None
        for targetid in targetids:

            # load elements using astroquery.jplhorizons
            obj = Horizons(id=targetid, id_type=id_type,
                           location=center, epochs=epochs)
            elem = obj.elements(**kwargs)

            # workaround for current version of astroquery to make
            # column units compatible with astropy.table.QTable
            # should really change '---' units to None in
            # astroquery.jplhorizons.__init__.py
            for column_name in elem.columns:
                if elem[column_name].unit == '---':
                    elem[column_name].unit = None

            if all_elem is None:
                all_elem = elem
            else:
                all_elem = vstack([all_elem, elem])

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Orbit', {'data service':
                                             '1996DPS....28.2504G'})

        return cls.from_table(all_elem)

    @classmethod
    def from_mpc(cls, targetid):
        """Load orbital elements from the Minor Planet Center
        (http://minorplanetcenter.net/).

        Parameters
        ----------
        targetid : str, mandatory
            target identifier

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> orb = Orbit.from_mpc('ceres') # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_astdys(cls, targetid):
        """Load orbital elements from AstDyS
        (http://hamilton.dm.unipi.it/astdys/).

        Parameters
        ----------
        targetid : str, mandatory
            target identifier

        Returns
        -------
        `~Orbit` object

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
        `~Orbit` object

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
        epoch : `~astropy.time.Time` object, mandatory
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
        `~Orbit` object

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
        `~Orbit` object

        Examples
        --------
        >>> from sbpy.data import Orbit # doctest: +SKIP
        >>> orb = Orbit.from... # doctest: +SKIP
        >>> sim = Orbit.integrate(orb, time=1000*u.year) # doctest: +SKIP
        >>> future_orb = Orbit.from_rebound(sim) # doctest: +SKIP

        not yet implemented

        """
