# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
sbpy.data
---------

:author: Michael Mommert (mommermiscience@gmail.com)
"""


class Conf():

    # property name alternatives for Orbit, Ephem, Phys
    # default name: [list of alternative names]
    namealts = {'i': ['inc', 'incl'],
                'epoch': ['datetime_jd'],
                'Omega': ['argper'],
                'w': ['longnode']
                }

    # reverse namealts for dict of alternative names pointing to
    # default names
    altnames = {}
    for key, vals in namealts.items():
        for val in vals:
            altnames[val] = key

    # definitions for use of pyoorb in Orbits
    oorb_timeScales = {'UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4}
    oorb_elemType = {'CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5}
    oorb_orbit_fields = {'COM': ['id', 'q', 'e', 'incl', 'Omega',
                                 'w', 'Tp_jd', 'orbtype', 'epoch',
                                 'epoch_scale', 'H', 'G'],
                         'KEP': ['id', 'a', 'e', 'incl', 'Omega', 'w', 'M',
                                 'orbtype', 'epoch', 'epoch_scale', 'H',
                                 'G'],
                         'CART': ['id', 'x', 'y', 'z', 'vx', 'vy', 'vz',
                                  'orbtype', 'epoch', 'epoch_scale', 'H',
                                  'G']}


conf = Conf()

from .core import (DataClass, mpc_observations, sb_search,
                   image_search, pds_ferret)
from .ephem import Ephem
from .orbit import Orbit
from .phys import Phys
from .names import Names

__all__ = ['DataClass', 'Ephem', 'Orbit', 'Phys', 'Names', 'conf', 'Conf',
           'mpc_observations', 'sb_search', 'image_search',
           'pds_ferret']
