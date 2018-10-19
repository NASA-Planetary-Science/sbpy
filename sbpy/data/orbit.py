# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Orbit Module
======================

Class for querying, manipulating, integrating, and fitting orbital elements.

created on June 04, 2017
"""
import os
from numpy import array, ndarray, double, arange
from astropy.time import Time
from astropy.table import vstack
from astroquery.jplhorizons import Horizons
import astropy.units as u

from .. import bib
from . import conf, DataClass

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
                    val.format = 'iso'
                    val.out_subfmt = 'date_hm'
                    epochs[key] = "'"+val.value+"'"
                else:
                    epochs[key] = "'"+epochs[key]+"'"
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
            bib.register('sbpy.data.Orbit.from_horizons',
                         {'data service': '1996DPS....28.2504G'})

        return cls.from_table(all_elem)

    @classmethod
    def from_mpc(cls, targetid):
        """Load orbital elements from the Minor Planet Center
        (http: // minorplanetcenter.net/).

        Parameters
        ----------
        targetid: str, mandatory
            target identifier

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        >> > from sbpy.data import Orbit  # doctest: +SKIP
        >> > orb = Orbit.from_mpc('ceres')  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_astdys(cls, targetid):
        """Load orbital elements from AstDyS
        (http: // hamilton.dm.unipi.it/astdys/).

        Parameters
        ----------
        targetid: str, mandatory
            target identifier

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        >> > from sbpy.data import Orbit  # doctest: +SKIP
        >> > orb = Orbit.from_mpc('ceres')  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_state(cls, pos, vel):
        """Convert state vector(positions and velocities) or orbital elements.

        Parameters
        ----------
        pos: `Astropy.coordinates` instance, mandatory
            positions vector
        vel: `Astropy.coordinates` instance, mandatory
            velocity vector

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        >> > from sbpy.data import Orbit  # doctest: +SKIP
        >> > import astropy.coordinates as coords  # doctest: +SKIP
        # doctest: +SKIP
        >> > r = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=1, y=0, z=0, unit=u.au))
        # doctest: +SKIP
        >> > v = coords.HeliocentricTrueEcliptic(coords.CartesianRepresentation(x=30, y=0, z=0, unit=u.km / u.s))
        >> > orb = Orbit.from_state(r, v)  # doctest: +SKIP

        not yet implemented

        """

    def to_state(self, epoch):
        """Convert orbital elements to state vector(positions and velocities)

        Parameters
        ----------
        epoch: `~astropy.time.Time` object, mandatory
          The epoch(s) at which to compute state vectors.

        Returns
        -------
        pos: `Astropy.coordinates` instance
            positions vector
        vel: `Astropy.coordinates` instance
            velocity vector

        Examples
        --------
        >> > from astropy.time import Time  # doctest: +SKIP
        >> > from sbpy.data import Orbit  # doctest: +SKIP
        >> > orb = Orbit.from_mpc('ceres')  # doctest: +SKIP
        >> > state = orb.to_state(Time('2015-03-06')  # doctest: +SKIP

        not yet implemented

        """

    def orbfit(self, eph):
        """Function that fits an orbit solution to a set of ephemerides using
        the OpenOrb(https: // github.com/oorb/oorb) software which has
        to be installed locally.

        Parameters
        - ---------
        eph: `Astropy.table`, mandatory
            set of ephemerides with mandatory columns `ra`, `dec`, `epoch` and
            optional columns `ra_sig`, `dec_sig`, `epoch_sig`

        additional parameters will be identified in the future

        Returns
        - ------ `~Orbit` object

        Examples
        - -------
        >> > from sbpy.data import Orbit, Ephem  # doctest: +SKIP
        >> > eph=Ephem.from_array([ra, dec, ra_sigma, dec_sigma,  # doctest: +SKIP
        >> >                         epochs, epochs_sigma],  # doctest: +SKIP
        >> >                         names=['ra', 'dec', 'ra_sigma',  # doctest: +SKIP
        >> >                                'dec_sigma', 'epochs',  # doctest: +SKIP
        >> >                                'epochs_sigma'])  # doctest: +SKIP
        >> > orb=Orbit.orbfit(eph)  # doctest: +SKIP

        not yet implemented

        """

    def integrate(self, time, integrator='IAS15'):
        """Function that integrates an orbit over a given range of time using
        the REBOUND(https: // github.com/hannorein/rebound) package

        Parameters
        - ---------
        time: `Astropy.units` quantity, mandatory
            Time range over which the orbit will be integrated.
        integrator: str, option, default 'IAS15'
            Integrator type to be used for the integration.

        Returns
        - ------
        REBOUND simulation object

        Examples
        - -------
        >> > from sbpy.data import Orbit  # doctest: +SKIP
        >> > orb=Orbit.from...  # doctest: +SKIP
        >> > sim=orb.integrate(1000*u.year)  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_rebound(cls, sim):
        """Obtain orbital elements from REBOUND
        (https: // github.com/hannorein/rebound) simulation instance.

        Parameters
        - ---------
        sim: REBOUND simulation instance, mandatory
            Simulation from which to obtain orbital elements.

        Returns
        - ------ `~Orbit` object

        Examples
        - -------
        >> > from sbpy.data import Orbit  # doctest: +SKIP
        >> > orb=Orbit.from...  # doctest: +SKIP
        >> > sim=Orbit.integrate(orb, time=1000*u.year)  # doctest: +SKIP
        >> > future_orb=Orbit.from_rebound(sim)  # doctest: +SKIP

        not yet implemented

        """

    # functions using pyoorb

    def _to_oo(self, timescale='UTC'):
        """Converts this orbit object to a openorb-compatible orbit array"""

        # identify orbit type based on available table columns
        orbittype = None
        for testtype in ['KEP', 'COM', 'CART']:
            try:
                self._translate_columns(
                    conf.oorb_orbit_fields[testtype][1:6])
                orbittype = testtype
                break
            except KeyError:
                pass

        if orbittype is None:
            raise ValueError(
                'orbit type cannot be determined from elements')

        # assemble orbit array for oorb_ephemeris
        if orbittype == 'COM':
            # cometary orbit: id q e i node argperi t_p otype epoch t H G
            orbits = array(array([arange(0, len(self.table), 1),
                                  self['q'].to('au').value,
                                  self['e'].data,
                                  self['i'].to('radian').value,
                                  self['Omega'].to('radian').value,
                                  self['w'].to('radian').value,
                                  (self['Tp_jd'].to('d').value -
                                   2400000.5),
                                  [conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  (self['epoch'].to('d').value
                                   - 2400000.5),
                                  [conf.oorb_timeScales[timescale]] *
                                  len(self.table),
                                  self['H'].value,
                                  self['G'].data]).transpose(),
                           dtype=double, order='F')
        elif orbittype == 'KEP':
            # keplerian orbit: id a e i node argperi M otype epoch ttype H G
            orbits = array(array([arange(0, len(self.table), 1),
                                  self['a'].to('au').value,
                                  self['e'].data,
                                  self['incl'].to('radian').value,
                                  self['Omega'].to('radian').value,
                                  self['w'].to('radian').value,
                                  self['M'].to('radian').value,
                                  [conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  (self['epoch'].to('d').value
                                   - 2400000.5),
                                  [conf.oorb_timeScales[timescale]] *
                                  len(self.table),
                                  self['H'].value,
                                  self['G'].data]).transpose(),
                           dtype=double, order='F')
        elif orbittype == 'CART':
            # cartesian orbit: id x y z dx dy dz otype epoch ttype H G
            orbits = array(array([arange(0, len(self.table), 1),
                                  self['x'].to('au').value,
                                  self['y'].to('au').value,
                                  self['z'].to('au').value,
                                  self['dx'].to('au/d').value,
                                  self['dy'].to('au/d').value,
                                  self['dz'].to('au/d').value,
                                  [conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  self['datetime_jd'].to('d').value
                                  - 2400000.5,
                                  [conf.oorb_timeScales[timescale]] *
                                  len(self.table),
                                  self['H'].data,
                                  self['G'].data]).transpose(),
                           dtype=double, order='F')

        return orbits

    def oo_transform(self, orbittype, timescale='UTC', ephfile='de430'):
        """Uses pyoorb to transform this orbit object to a different
        orbit type definition.

        Parameters
        ----------
        orbittype : str
            Orbit definition to be transformed to; available orbit
            definitions are ``KEP`` (Keplerian elements), ``CART``
            (cartesian elements), ``COM`` (cometary elements).
        timescale : str
            Timescale to be used in the transformation; the following
            values are allowed: ``'UTC'``, ``'UT1'``, ``'TT'``,
            ``'TAI'``. Default: ``'UTC'``
        ephfile : str, optional
            Planet and Lunar ephemeris file version as provided by JPL
            to be used in the propagation. Default: ``'de430'``

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        Obtain the current state vector (cartesian definition, ``CART``) for
        asteroid Ceres.

        >>> from sbpy.data import Orbit
        >>> ceres = Orbit.from_horizons('Ceres')
        >>> statevec = ceres.oo_transform('CART')  # doctest: +SKIP
        >>> print(statevec.table)  # doctest: +SKIP
           id           x                   y           ... epoch_scale  H    G
                        AU                  AU          ...             mag
        ------- ------------------ -------------------- ... ----------- ---- ----
        1 Ceres -2.525066589355839 -0.34964875372044524 ...         UTC 3.34 0.12

        """
        import pyoorb

        # initialize pyoorb
        ephfile = os.path.join(os.getenv('OORB_DATA'), ephfile+'.dat')
        pyoorb.pyoorb.oorb_init(ephfile)

        oo_orbits, err = pyoorb.pyoorb.oorb_element_transformation(
            in_orbits=self._to_oo(timescale),
            in_element_type={'CART': 1, 'COM': 2, 'KEP': 3,
                             'DEL': 4, 'EQX': 5}[orbittype])

        if err != 0:
            RuntimeError('pyoorb failed with error code {:d}'.format(err))

        # reorder data in Orbit object
        columns = conf.oorb_orbit_fields[orbittype]

        orbits = self.from_array(oo_orbits.transpose(), names=columns)

        # apply units
        for i, col in enumerate(orbits.column_names):
            orbits[col].unit = conf.oorb_orbit_units[orbittype][i]

        # replace id column with actual target names from original orbits
        orbits.table.replace_column('id', self['targetname'])

        # replace orbtype and epoch_scale columns
        orbits.table.replace_column('orbtype',
                                    [orbittype] * len(orbits.table))
        orbits.table.replace_column('epoch_scale',
                                    [timescale] * len(orbits.table))

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem.from_oo',
                         {'method': '2009M&PS...44.1853G',
                          'implementation': 'https://github.com/oorb/oorb'})

        return orbits

    def oo_propagate(self, epoch, timescale='UTC',
                     dynmodel='N', ephfile='de430'):
        """Uses pyoorb to propagate this `~Orbit` object.

        Parameters
        ----------
        epoch : `~astropy.time.Time` object or float
            Epoch to which the orbit will be propagated to. A float
            value will be interpreted as Julian date.
        timescale : str, optional
            Timescale to be used in the propagation; the following
            values are allowed: ``'UTC'``, ``'UT1'``, ``'TT'``,
            ``'TAI'``. Default: ``'UTC'``
        dynmodel : str, optional
            The dynamical model to be used in the propagation: ``'N'``
            for n-body simulation or ``'2'`` for a 2-body
            simulation. Default: ``'N'``
        ephfile : str, optional
            Planet and Lunar ephemeris file version as provided by JPL
            to be used in the propagation. Default: ``'de430'``

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        Propagate the orbit of Ceres 100 days into the future.
        >>> from sbpy.data import Orbit
        >>> from astropy.time import Time
        >>> epoch = Time.now().jd + 100
        >>> ceres = Orbit.from_horizons('Ceres')
        >>> future_ceres = ceres.oo_propagate(epoch)  # doctest: +SKIP
        >>> print(future_ceres.table)  # doctest: +SKIP
           id           a                  e          ... epoch_scale  H    G  
                        AU                            ...             mag      
        ------- ----------------- ------------------- ... ----------- ---- ----
        1 Ceres 2.767911178119476 0.07574650026062148 ...         UTC 3.34 0.12
        """

        import pyoorb

        # initialize pyoorb
        ephfile = os.path.join(os.getenv('OORB_DATA'), ephfile+'.dat')
        pyoorb.pyoorb.oorb_init(ephfile)

        # identify orbit type based on available table columns
        orbittype = None
        for testtype in ['KEP', 'COM', 'CART']:
            try:
                self._translate_columns(
                    conf.oorb_orbit_fields[testtype][1:6])
                orbittype = testtype
                break
            except KeyError:
                pass

        if orbittype is None:
            raise ValueError(
                'orbit type cannot be determined from elements')

        if isinstance(epoch, Time):
            ooepoch = [epoch.jd-2400000.5, conf.oorb_timeScales[timescale]]
        else:
            ooepoch = [epoch-2400000.5, conf.oorb_timeScales[timescale]]

        oo_orbits, err = pyoorb.pyoorb.oorb_propagation(
            in_orbits=self._to_oo(timescale),
            in_epoch=ooepoch,
            in_dynmodel=dynmodel)

        if err != 0:
            RuntimeError('pyoorb failed with error code {:d}'.format(err))

        # reorder data in Orbit object
        columns = conf.oorb_orbit_fields[orbittype]

        orbits = self.from_array(oo_orbits.transpose(), names=columns)

        # apply units
        for i, col in enumerate(orbits.column_names):
            orbits[col].unit = conf.oorb_orbit_units[orbittype][i]

        # replace id column with actual target names from original orbits
        orbits.table.replace_column('id', self.table['targetname'])

        # replace orbtype and epoch_scale columns
        orbits.table.replace_column('orbtype',
                                    [orbittype] * len(orbits.table))
        orbits.table.replace_column('epoch_scale',
                                    [timescale] * len(orbits.table))

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem.from_oo',
                         {'method': '2009M&PS...44.1853G',
                          'implementation': 'https://github.com/oorb/oorb'})

        return orbits
