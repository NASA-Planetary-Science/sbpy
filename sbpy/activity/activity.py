# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
SBPy Activity Module
====================

created on June 23, 2017
"""

__all__ = ['Activity']


class Activity():

    class Haser():
        """Haser model implementation"""

        def __init__(self, Q, gamma):
            """Parameters
            ----------
            Q : `Astropy.units` quantity or iterable, mandatory
                production rate usually in units of `u.molecule / u.s`
            gamma : `Astropy.units` quantity or iterable, mandatory
                scalen length usually in units of `u.km`
            """

            
        def volume_density(self, aperture, eph):
            """calculates coma volume density"""
