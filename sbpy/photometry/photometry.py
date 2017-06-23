# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
SBPy Photometry Module
======================

created on June 23, 2017
"""

__all__ = ['Photometry']


class Photometry():

    def diam2mag(phys, eph, model=None):
        """Function to calculate the apparent bightness of a body from its physical properties and ephemerides"""
