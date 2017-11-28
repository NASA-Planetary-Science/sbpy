# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Photometry Module

created on June 23, 2017
"""

from astropy.modeling import Fittable1DModel, Fittable2DModel, Parameter
import numpy as np

__all__ = ['DiskIntegratedModelClass', 'HG', 'HG12', 'HG1G2',
           'DiskFunctionModel', 'LommelSeeliger', 'Lambert',
           'LunarLambert', 'PhaseFunctionModel', 'ROLOPhase',
           'ResolvedPhotometricModelClass', 'ROLO']


class DiskIntegratedModelClass(Fittable1DModel):
    """Base class for disk-integrated phase function model"""

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


class DiskFunctionModel(Fittable2DModel):
    """Base class for disk-function model"""
    pass


class PhaseFunctionModel(Fittable1DModel):
    """Base class for phase function model"""
    pass


class HG(DiskIntegratedModelClass):
    """IAU HG photometric phase model (Bowell XXX)"""
    H = Parameter(default=3.2, description='H parameter')
    G = Parameter(default=0.28, description='G parameter')

 HG12(DiskIntegratedModelClass):
    """IAU HG12 photometric phase model (Muinonen et al. 2010)"""
    H = Parameter(default=3.2, description='H parameter')
    G = Parameter(default=0.2, description='G12 parameter')


class HG1G2(DiskIntegratedModelClass):
    """IAU HG1G2 photometric phase model (Muinonen et al. 2010)"""
    H = Parameter(default=3.2, description='H parameter')
    G1 = Parameter(default=0.2, description='G1 parameter')
    G2 = Parameter(default=0.2, description='G2 parameter')


class LommelSeeliger(DiskFunctionModel):
    """Lommel-Seeliger model class"""

    @staticmethod
    def evaluate(i, e):
        mu0 = np.cos(i)
        mu = np.cos(e)
        return mu0/(mu0+mu)


class Lambert(DiskFunctionModel):
    """Lambert model class"""

    @staticmethod
    def evaluate(i, e):
        return np.cos(i)


class LunarLambert(DiskFunctionModel):
    """Lunar-Lambert model, or McEwen model class"""
    L = Parameter(default=0.2, description='Partition parameter')

    @staticmethod
    def evaluate(i, e, L):
        mu0 = np.cos(i)
        mu = np.cos(e)
        return (1-L)*mu0/(mu0+mu) + L*mu0

    @staticmethod
    def fit_deriv(i, e, L):
        return -mu0/(mu0+mu) + mu0


class ROLOPhase(PhaseFunctionModel):
    """ROLO phase function model class"""
    A0 = Parameter(default=0.1, description='ROLO A0 parameter')
    A1 = Parameter(default=0.1, description='ROLO A1 parameter')
    C0 = Parameter(default=0.1, description='ROLO C0 parameter')
    C1 = Parameter(default=0.1, description='ROLO C1 parameter')
    C2 = Parameter(default=0.1, description='ROLO C2 parameter')
    C3 = Parameter(default=0.1, description='ROLO C3 parameter')
    C4 = Parameter(default=0.1, description='ROLO C4 parameter')


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
