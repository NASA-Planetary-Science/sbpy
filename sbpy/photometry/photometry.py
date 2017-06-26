# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
SBPy Photometry Module
======================

created on June 23, 2017
"""

__all__ = ['ModelClass', 'HG', 'HG12', 'HG1G2']

class ModelClass():

    def fit(self, eph):
        """Fit photometric model to photometric data stored in sbpy.data.Ephem
        object
        
        Parameters
        ----------
        eph : `sbpy.data.Ephem` instance, mandatory
            photometric data; must contain `phase` (phase angle) and `mag` 
            (apparent magnitude) columns; `mag_err` optional
        
        Returns
        -------
        fit Chi2
        
        Examples
        --------
        >>> from sbpy.photometry import HG
        >>> from sbpy.data import Misc
        >>> eph = Misc.mpc_observations('Bennu')
        >>> hg = HG()
        >>> chi2 = hg.fit(eph)

        not yet implemented

        """

    def mag(self, phase):
        """Derive magnitude for a given photometric model as a function of 
        phase angle
        
        Parameters
        ----------
        phase : list or array, mandatory
            phase angles 

        Returns
        -------
        sbpy.data.Ephem instance 
        
        Examples
        --------
        >>> from sbpy.photometry import HG1G2
        >>> hg1g2 = HG1G2(12, 0.1, -0.2)
        >>> mag = hg1g2.mag([0, 5, 15, 30 ,60, 90])

        not yet implemented

        """

    def distance_module(self, eph):
        """Account magnitudes for distance module (distance from observer, 
        distance to the Sun); return modified magnitudes

        Parameters
        ----------
        eph : list or array, mandatory
            phase angles 

        Returns
        -------
        sbpy.data.Ephem instance
        
        Examples
        --------
        TBD

        not yet implemented

        """
        

class HG(ModelClass):
    """HG photometric phase model (Bowell XXX)"""
    def __init__(self, **kwargs):
        self.H = None
        self.G = None

class HG12(ModelClass):
    """HG12 photometric phase model (Muinonen et al. 2010)"""
    def __init__(self, **kwargs):
        self.H = None
        self.G12 = None

class HG1G2(ModelClass):
    """HG1G2 photometric phase model (Muinonen et al. 2010)"""
    def __init__(self, **kwargs):
        self.H = None
        self.G1 = None
        self.G2 = None
        


# class Photometry():

#     def diam2mag(phys, eph, model=None):
#         """Function to calculate the apparent bightness of a body from its physical properties and ephemerides"""
