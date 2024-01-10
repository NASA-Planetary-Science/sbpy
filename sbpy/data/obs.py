# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
sbpy data.Obs Module
====================

Class for storing and querying observations

created on July 3, 2019
"""

from astropy.time import Time

try:
    from astroquery.mpc import MPC
except ImportError:
    pass

from astropy.table import vstack, hstack

from .ephem import Ephem
from .core import QueryError
from ..bib import cite
from ..utils.decorators import requires

__all__ = ['Obs']


class Obs(Ephem):
    """Class for querying, storing, and manipulating observations """

    @classmethod
    @requires("astroquery")
    @cite({'data source': 'https://minorplanetcenter.net/db_search',
           'software: astroquery': '2019AJ....157...98G'})
    def from_mpc(cls, targetid, id_type=None, **kwargs):
        """Load available observations for a target from the
        `Minor Planet Center <https://minorplanetcenter.net>`_ using
        `~astroquery.mpc.MPCClass.get_observations`.

        Parameters
        ----------
        targetid : str
            Target identifier, resolvable by the Minor Planet
            Ephemeris Service.

        id_type : str, optional
            ``'asteroid number'``, ``'asteroid designation'``,
            ``'comet number'``, ``'comet designation'``, or ``None``.
            ``None`` attempts to automatically identify the target id type
            using `~sbpy.data.Names`. Default: ``None``

        **kwargs : optional
            Additional keyword arguments are passed to
            `~astroquery.mpc.MPC.query_object`

        Returns
        -------
        `~Obs` object
            The resulting object will be populated with the same fields
            defined in `~astroquery.mpc.MPCClass.get_observations`.


        Examples
        --------
        >>> from sbpy.data import Obs
        >>> obs = Obs.from_mpc('12893') # doctest: +REMOTE_DATA
        >>> orb[:3]  # doctest: +SKIP
        number   desig   discovery note1 ...         DEC         mag band observatory
                                         ...         deg         mag
        ------ --------- --------- ----- ... ------------------- --- ---- -----------
         12893 1998 QS55        --    -- ...  -15.78888888888889 0.0   --         413
         12893 1998 QS55        --    -- ... -15.788944444444445 0.0   --         413
         12893 1998 QS55         *     4 ...   5.526472222222222 0.0   --         809
        """

        from ..data import Names

        if id_type is None:
            id_type = Names.asteroid_or_comet(targetid)
            if id_type == 'asteroid':
                ident = Names.parse_asteroid(targetid)
            elif id_type == 'comet':
                ident = Names.parse_comet(targetid)
            if 'number' in ident:
                id_type += ' number'
            elif 'designation' in ident:
                id_type += ' designation'

        try:
            results = MPC.get_observations(targetid, id_type=id_type,
                                           **kwargs)
        except (RuntimeError, ValueError) as e:
            raise QueryError(
                ('Error raised by '
                 'astroquery.mpc.MPCClass.get_observations: '
                 '{}').format(e))

        results['epoch'] = Time(results['epoch'].to('d').value,
                                scale='utc', format='jd')

        return cls.from_table(results)

    @requires("astroquery")
    @cite({'software: astroquery': '2019AJ....157...98G'})
    def supplement(self, service='jplhorizons', id_field='targetname',
                   epoch_field='epoch', location='500',
                   modify_fieldnames='obs', **kwargs):
        """Supplement observational data with ephemerides
        queried from the selected service.

        Parameters
        ----------
        service : str, optional
            Service from which to acquire data: ``'jplhorizons'``,
            ``'mpc'``, or ``'miriade'``, corresponding to the
            `JPL Horizons system <https://ssd.jpl.nasa.gov/horizons/>`_
            (using `~sbpy.data.Ephem.from_horizons`),
            the `Minor Planet Center ephemeris service
            <https://minorplanetcenter.net/iau/MPEph/MPEph.html>`_
            (using `~sbpy.data.Ephem.from_mpc`), and
            the `IMCCE Miriade service
            <http://vo.imcce.fr/webservices/miriade/>`_
            (using `~sbpy.data.from_miriade`). Default:
            ``'jplhorizons'``
        id_field : str, optional
            Field name that corresponds to a suitable target identifier in
            this `~sbpy.data.Obs` object. Default: ``'targetname'``
        epoch_field : str, optional
            Field name that corresponds to a suitable epoch identifier in
            this `~sbpy.data.Obs` object. The corresponding column must be
            of type `~astropy.time.Time`. Default: ``'epoch'``
        location : str, optional
            Location of the observer for the data stored in this
            `~sbpy.data.Obs` object. Default: ``'500'`` (geocenter)
        modify_fieldnames : str, optional
            Defines whether field names in this `~sbpy.data.Obs` object
            (``'obs'``) or in the supplemental data to be queried (``'eph'``)
            will be modified by adding a suffix in case of field name
            collisions. Default: ``'obs'``
        **kwargs : optional
            Additional keyword arguments are passed to the corresponding
            ephemerides query service.

        Returns
        -------
        `~Obs` object
            The resulting object will contain all data from this
            `~sbpy.data.Obs` object as well as the queried ephemeris data.

        Notes
        -----
        * Not all available service are equally suited for this kind
          of query: only the JPL Horizons system enables quick queries
          for a large number of epochs. Queries using the other
          services may take a long time depending on the number of
          epochs and targets.


        Examples
        --------
        >>> from sbpy.data import Obs
        >>> obs = Obs.from_mpc('2019 AA', id_type='asteroid designation') # doctest: +SKIP
        >>> data = obs.supplement(id_field='designation') # doctest: +SKIP
        >>> data.field_names # doctest: +SKIP
        <TableColumns names=('number','desig','discovery','note1','note2','epoch','RA_obs','DEC_obs','mag','band','observatory','target','RA','DEC','delta','V','alpha','elong','RAcosD_rate','DEC_rate','delta_rate')>
        """

        try:
            targetids = set(self[id_field])
        except (TypeError, KeyError):
            raise QueryError('cannot use field {} as id_field.'.format(
                id_field))

        all_obs = None
        all_eph = None
        for targetid in targetids:

            if all_obs is None:
                all_obs = self.table[self[id_field] == targetid]
            else:
                all_obs = vstack([all_obs,
                                  self.table[self[id_field] == targetid]])

            if service == 'jplhorizons':
                eph = Ephem.from_horizons(
                    targetid,
                    epochs=self[self[id_field] == targetid][epoch_field],
                    location=location,
                    **kwargs)
                eph.table.remove_column('epoch')
            elif service == 'mpc':
                eph = Ephem.from_mpc(
                    targetid,
                    epochs=self[self[id_field] == targetid][epoch_field],
                    location=location,
                    **kwargs)
                eph.table.remove_column('Date')
            elif service == 'miriade':
                eph = Ephem.from_miriade(
                    targetid,
                    epochs=self[self[id_field] == targetid][epoch_field],
                    location=location,
                    **kwargs)
                eph.table.remove_column('epoch')
            else:
                raise QueryError('service {} not known.'.format(service))

            if all_eph is None:
                all_eph = eph.table
            else:
                all_eph = vstack([all_eph, eph.table])

        # identify field names that both obs and eph have in common
        fieldnames_intersect = set(all_eph.columns).intersection(
            all_obs.columns)
        for fieldname in fieldnames_intersect:
            if modify_fieldnames == 'obs':
                all_obs.rename_column(fieldname, fieldname+'_obs')
            elif modify_fieldnames == 'eph':
                all_eph.rename_column(fieldname, fieldname+'_eph')

        return Obs.from_table(hstack([all_obs, all_eph]),
                              meta=self.meta)
