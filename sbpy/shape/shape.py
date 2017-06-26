# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=================
SBPy Shape Module
=================

created on June 26, 2017
"""

__all__ = ['ROLO', 'Kaasalainen']

class ModelClass():

    def __init__(self):
        self.shape = None
    
    def load_obj(self, filename):
        """Load .OBJ shape model"""

    def get_facet(self, facetid):
        """Retrieve information for one specific surface facet"""

    def iof(self, facetid, ephem):
        """Derive I/F for one specific surface facet"""

    def integrated_flux(self, ephem):
        """Derive surface-integrated flux"""

    def lightcurve(self, ephem, time):
        """Derive lightcurve"""
        
        
class ROLO(ModelClass):

    def __init__(self):
        self.properties = None
        # setup model properties



class Kaasalainen(ModelClass):

    def __init__(self):
        self.properties = None
        # setup model properties

    def convexinv(self, ephem):
        """Obtain shape model through lightcurve inversion, optimizing all
        parameters and uses spherical harmonics functions for shape
        representation"""

    def conjgradinv(self, ephem):
        """Obtain shape model through lightcurve inversion, optimizing only
        shape and uses directly facet areas as parameters"""
