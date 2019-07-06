# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
sbpy data.Obs Module
====================

Class for storing and querying observations

created on July 3, 2019
"""

from .ephem import Ephem

__all__ = ['Obs']


class Obs(Ephem):
    """Class for querying, storing, and manipulating observations """

    pass
