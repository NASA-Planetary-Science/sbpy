# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
sbpy data.Obs Module
====================

Class for storing and querying observations

created on July 3, 2019
"""

from astroquery.mpc import MPC

from .ephem import Ephem
from ..bib import cite
from .core import QueryError

__all__ = ['Obs']


class Obs(Ephem):
    """Class for querying, storing, and manipulating observations """

    @classmethod
    @cite({'data source':
           'https://minorplanetcenter.net/db_search'})
    def from_mpc(cls, targetid, id_type=None, **kwargs):
        """Load available observations for a target from the
        `Minor Planet Center <http://minorplanetcenter.net>`_ using
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

        return cls.from_table(results)
