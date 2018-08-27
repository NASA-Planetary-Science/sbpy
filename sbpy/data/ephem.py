# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Ephem Module
======================

Class for storing and querying ephemerides

created on June 04, 2017
"""
from numpy import ndarray
from astropy.time import Time
from astropy.table import vstack
from astroquery.jplhorizons import Horizons

from .. import bib
from .core import DataClass

__all__ = ['Ephem']


class Ephem(DataClass):
    """Class for querying, manipulating, and calculating ephemerides"""

    @classmethod
    def from_horizons(cls, targetids, id_type='smallbody',
                      epochs=None, location='500', **kwargs):
        """Load target ephemerides from
        `JPL Horizons <https://ssd.jpl.nasa.gov/horizons.cgi>`_ using
        `astroquery.jplhorizons.HorizonsClass.ephemerides`

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
            `~astroquery.jplhorizons.HorizonsClass.Horizons.ephemerides`
            for details); if ``None`` is provided, current date and
            time are used. Default: ``None``
        location : str, optional, default ``'500'`` (geocentric)
            Location of the observer.
        **kwargs : optional
            Arguments that will be provided to
            `astroquery.jplhorizons.HorizonsClass.ephemerides`.

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_horizons('ceres', epochs=epoch)
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
        elif isinstance(epochs, (list, tuple, ndarray)):
            new_epochs = [None] * len(epochs)
            for i in range(len(epochs)):
                if isinstance(epochs[i], Time):
                    new_epochs[i] = epochs[i].jd
                else:
                    new_epochs[i] = epochs[i]
            epochs = new_epochs

        # if targetids is a list, run separate Horizons queries and append
        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        # append ephemerides table for each targetid
        all_eph = None
        for targetid in targetids:

            # load ephemerides using astroquery.jplhorizons
            obj = Horizons(id=targetid, id_type=id_type,
                           location=location, epochs=epochs)
            eph = obj.ephemerides(**kwargs)

            # workaround for current version of astroquery to make
            # column units compatible with astropy.table.QTable
            # should really change '---' units to None in
            # astroquery.jplhorizons.__init__.py
            for column_name in eph.columns:
                if eph[column_name].unit == '---':
                    eph[column_name].unit = None

            if all_eph is None:
                all_eph = eph
            else:
                all_eph = vstack([all_eph, eph])

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem', {'data service':
                                             '1996DPS....28.2504G'})

        return cls.from_table(all_eph)

    @classmethod
    def from_mpc(cls, targetid, epoch, observatory='500'):
        """
        Load ephemerides from the
        `Minor Planet Center <http://minorplanetcenter.net>`_.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default ``None``
            epoch of elements; if ``None`` is provided, current date is used
        observatory : str, optional, default ``'500'`` (geocentric)
            location of observer

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> from astropy.time import Time  # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc')  # doctest: +SKIP
        >>> eph = Ephem.from_mpc('ceres', '568', epoch)  # doctest: +SKIP

        not yet implemented

        """

    def report_to_mpc():
        """
        Format ephemerides as a report to the
        `Minor Planet Center <http://minorplanetcenter.net>`_.

        Returns
        -------
        list of strings

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> eph = Ephem.from_array...  # doctest: +SKIP
        >>> report = eph.report_to_mpc()  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_imcce(cls, targetid, epoch, observatory='500'):
        """
        Load orbital elements from
        `IMCCE <http://vo.imcce.fr/webservices/miriade/>`_.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default ``None``
            epoch of elements; if ``None`` is provided, current date is used
        observatory : str, optional, default '500' (geocentric)
            location of observer

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> from astropy.time import Time  # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc')  # doctest: +SKIP
        >>> eph = Ephem.from_imcce('ceres', '568', epoch)  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_lowell(cls, targetid, epoch, observatory='500'):
        """
        Load orbital elements from
        Lowell Observatory's `astorb <http://asteroid.lowell.edu/>`_.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default ``None``
            epoch of elements; if ``None`` is provided, current date is used
        observatory : str, optional, default '500' (geocentric)
            location of observer

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> from astropy.time import Time  # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc')  # doctest: +SKIP
        >>> eph = Ephem.from_lowell('ceres', '568', epoch)  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_pyephem(cls, orb, location, epoch):
        """
        Derives ephemerides using
        `PyEphem <http://rhodesmill.org/pyephem/>`_.

        Parameters
        ----------
        orb : ``~Orbit`` object, mandatory
            complete set of orbital elements
        location : str or dictionary, mandatory
            name of location or a dictionary fully describing the location
        epoch : astropy ``Time`` object

        Examples
        --------
        >>> from sbpy.data import Ephem, Orbit  # doctest: +SKIP
        >>> orb = Orbit.from_...  # doctest: +SKIP
        >>> eph = Ephem.from_pyephem(orb,  # doctest: +SKIP
        >>>                          location={'name': 'Flagstaff',  # doctest: +SKIP
        >>>                                    'geolon': 35.199167,  # doctest: +SKIP
        >>>                                    'geolat': -111.631111,  # doctest: +SKIP
        >>>                                    'altitude': '2106'},  # doctest: +SKIP
        >>>                          epoch=epoch)  # doctest: +SKIP

        not yet implemented

        """
