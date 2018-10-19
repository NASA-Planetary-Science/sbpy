# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
sbpy.data
---------

:author: Michael Mommert (mommermiscience@gmail.com)
"""


class Conf():

    # alternative field names for DataClass
    # append field description as final item in the list (for documentation)
    fieldnames = [
        ['targetname', 'id', 'Target Identifier'],

        # orbital elements
        ['i', 'inc', 'incl', 'Inclination'],
        ['epoch', 'datetime_jd', 'Date', 'date', 'Epoch'],
        ['Omega', 'longnode', 'Longitude of the Ascending Node'],
        ['w', 'argper', 'Argument of the Periapsis'],

        # ephemerides
        ['r', 'r_hel', 'heldist', 'Heliocentric Distance'],
        ['delta', 'Delta', 'obsdist', 'Distance to the Observer'],
        ['ra', 'RA', 'Right Ascension'],
        ['dec', 'DEC', 'Dec', 'Declination'],
        ['ra_rate', 'RA_rate', 'ra_rates', 'RA_rates', 'dRA',
         'dra', 'RA Rate'],
        ['dec_rate', 'DEC_rate', 'Dec_rate', 'dec_rates', 'DEC_rates',
         'Dec_rates', 'dDec', 'dDEC', 'ddec', 'Dec Rate'],
        ['alpha', 'phaseangle', 'Phase', 'Solar Phase Angle'],
        ['elong', 'solarelong', 'solarelongation', 'elongation',
         'Elongation', 'Solar Elongation'],
        ['V', 'Vmag', 'V-band Magnitude'],
        ['hlon', 'EclLon', 'ecllon', 'HelEclLon',
         'helecllon', 'Heliocentric Ecliptic Longitude'],
        ['hlat', 'EclLat', 'ecllat', 'HelEclLat', 'helecllat',
         'Heliocentric Ecliptic Latitude'],
        ['el', 'EL', 'elevation', 'alt',
         'altitude', 'Elevation'],
        ['lunar_elong', 'elong_moon', 'elongation_moon',
         'lunar_elongation', 'lunarelong', 'Lunar Elongation'],
        ['vx', 'dx', 'dx/dt', 'x Velocity Component'],
        ['vy', 'dy', 'dy/dt',  'y Velocity Component'],
        ['vz', 'dz', 'dz/dt', 'z Velocity Component'],

        # physical properties
        ['d', 'D', 'diam', 'Diameter'],
        ['pv', 'pV', 'p_v', 'p_V', 'V-band Geometric Albedo'],
    ]

    fieldname_idx = {}
    for idx, field in enumerate(fieldnames):
        for alt in field:
            fieldname_idx[alt] = idx

    # property name alternatives for Orbit, Ephem, Phys
    # default name: [list of alternative names]
    # namealts = {
        #     'targetname': ['id'],  # target identifier
        #     # orbital elements
        #     'i': ['inc', 'incl'],  # inclination
        #     'epoch': ['datetime_jd', 'Date', 'date'],  # epoch
        #     'Omega': ['longnode'],  # longitude of the ascending node
        #     'w': ['argper'],  # argument of periapsis
        #     # ephemerides
        #     'r': ['r_hel', 'heldist'],  # heliocentric distance
        #     'delta': ['Delta', 'obsdist'],  # distance to observer
        #     'ra': ['RA'],  # right ascension
        #     'dec': ['DEC', 'Dec'],  # declination
        #     'ra_rate': ['RA_rate', 'ra_rates', 'RA_rates',
        #                 'dRA', 'dra'],  # RA rate
        #     'dec_rate': ['DEC_rate', 'Dec_rate', 'dec_rates', 'DEC_rates',
        #                  'Dec_rates', 'dDec', 'dDEC', 'ddec'],  # DEC rate
        #     'alpha': ['phaseangle', 'Phase'],  # solar phase angle
        #     'elong': ['solarelong', 'solarelongation',
        #               'elongation', 'Elongation'],  # solar elongation
        #     'V': ['Vmag'],  # V-band magnitude
        #     'hlon': ['EclLon', 'ecllon',
        #              'HelEclLon', 'helecllon'],  # heliocentric ecliptic longitude
        #     'hlat': ['EclLat', 'ecllat',
        #              'HelEclLat', 'helecllat'],  # heliocentric ecliptic latitude
        #     'el': ['EL', 'elevation', 'alt', 'altitude'],  # topocentric elevation
        #     'lunar_elong': ['elong_moon', 'elongation_moon',
        #                     'lunar_elongation', 'lunarelong'],  # lunar elongation
        #     'vx': ['dx', 'dx/dt'],  # x velocity component
        #     'vy': ['dy', 'dy/dt'],  # x velocity component
        #     'vz': ['dz', 'dz/dt'],  # x velocity component

        #     # physical properties
        #     'd': ['D', 'diam'],  # diameter
        #     'pv': ['pV', 'p_v', 'p_V'],  # V-band geometric albedo
        # }

        # # reverse namealts for dict of alternative names pointing to
        # # default names
        # altnames = {}
        # for key, vals in namealts.items():
        #     for val in vals:
        #         altnames[val] = key

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
    oorb_orbit_units = {'COM': [None, 'au', None, 'deg', 'deg',
                                'deg', 'd', None, 'd',
                                None, 'mag', None],
                        'KEP': [None, 'au', None, 'deg', 'deg', 'deg', 'deg',
                                None, 'd', None, 'mag', None],
                        'CART': [None, 'au', 'au', 'au', 'au/d', 'au/d',
                                 'au/d', None, 'd', None, 'mag', None]}

    oorb_ephem_fields = ['MJD', 'RA', 'DEC', 'RA_rate', 'DEC_rate', 'alpha',
                         'elong', 'r', 'Delta', 'V', 'pa', 'TopEclLon',
                         'TopEclLat', 'OppTopEclLon', 'OppTopEclLat',
                         'HelEclLon', 'HelEclLat', 'OppHelEclLon',
                         'OppHelEclLat', 'EL', 'ELsun', 'ELmoon',
                         'lunarphase', 'lunarelong', 'x', 'y', 'z',
                         'vx', 'vy', 'vz', 'obsx', 'obsy', 'obsz',
                         'trueanom']
    oorb_ephem_units = ['d', 'deg', 'deg', 'deg/d', 'deg/d', 'deg',
                        'deg', 'au', 'au', 'mag', 'deg', 'deg',
                        'deg', 'deg', 'deg',
                        'deg', 'deg', 'deg',
                        'deg', 'deg', 'deg', 'deg',
                        None, 'deg', 'au', 'au', 'au',
                        'au/d', 'au/d', 'au/d', 'au', 'au', 'au', 'deg']


conf = Conf()

from .core import (DataClass, mpc_observations, sb_search,
                   image_search, pds_ferret)
from .ephem import Ephem
from .orbit import Orbit
from .phys import Phys
from .names import Names, natural_sort_key

__all__ = ['DataClass', 'Ephem', 'Orbit', 'Phys', 'Names', 'conf', 'Conf',
           'mpc_observations', 'sb_search', 'image_search',
           'pds_ferret']
