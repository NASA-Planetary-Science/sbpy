# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
SBPy Photometry Module
======================

created on June 23, 2017
"""

__all__ = ['DiskIntegratedModelClass', 'HG', 'HG12', 'HG1G2', 'DiskFunctionModel', 'LommelSeeliger', 'Lambert', 'PhaseFunctionModel', 'ROLOPhase', 'ResolvedPhotometricModelClass', 'ROLO']

class DiskIntegratedModelClass():

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
        >>> from sbpy.photometry import HG # doctest: +SKIP
        >>> from sbpy.data import Misc # doctest: +SKIP
        >>> eph = Misc.mpc_observations('Bennu') # doctest: +SKIP
        >>> hg = HG() # doctest: +SKIP
        >>> chi2 = hg.fit(eph) # doctest: +SKIP

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
        >>> hg1g2 = HG1G2(12, 0.1, -0.2) # doctest: +SKIP
        >>> mag = hg1g2.mag([0, 5, 15, 30 ,60, 90]) # doctest: +SKIP

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


class HG(DiskIntegratedModelClass):
    """HG photometric phase model (Bowell XXX)"""
    def __init__(self, **kwargs):
        self.H = None
        self.G = None

class HG12(DiskIntegratedModelClass):
    """HG12 photometric phase model (Muinonen et al. 2010)"""
    def __init__(self, **kwargs):
        self.H = None
        self.G12 = None

class HG1G2(DiskIntegratedModelClass):
    """HG1G2 photometric phase model (Muinonen et al. 2010)"""
    def __init__(self, **kwargs):
        self.H = None
        self.G1 = None
        self.G2 = None


class DiskFunctionModel(object):
    """Base class for disk-function model"""
    pass


class LommelSeeliger(DiskFunctionModel):
    """Lommel-Seeliger model class"""
    pass


class Lambert(DiskFunctionModel):
    """Lambert model class"""
    pass


class PhaseFunctionModel(object):
    """Base class for phase function model"""
    pass


class ROLOPhase(PhaseFunctionModel):
    """ROLO phase function model class"""
    A0 = None
    A1 = None
    C0 = None
    C1 = None
    C2 = None
    C3 = None
    C4 = None


class ResolvedPhotometricModelClass(object):
    """Base class for disk-resolved photometric model"""
    # composite model as the product of a disk function and a phase function
    pass

class ROLO(ResolvedPhotometricModelClass):
    """ROLO disk-resolved photometric model"""
    pass


# class Photometry():

#     def diam2mag(phys, eph, model=None):
#         """Function to calculate the apparent bightness of a body from its physical properties and ephemerides"""
