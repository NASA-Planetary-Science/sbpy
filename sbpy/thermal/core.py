# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Thermal Module

created on June 27, 2017
"""

__all__ = ['ThermalClass', 'STM', 'FRM', 'NEATM']

__doctest_requires__ = {
    "ThermalClass.flux": "astroquery"
}


class ThermalClass():

    def flux(phys, eph, lam):
        """Model flux density for a given wavelength `lam`, or a list/array thereof

        Parameters
        ----------
        phys : `sbpy.data.Phys` instance, mandatory
            provide physical properties
        eph : `sbpy.data.Ephem` instance, mandatory
            provide object ephemerides
        lam : `astropy.units` quantity or list-like, mandatory
            wavelength or list thereof

        Examples
        --------
        >>> from astropy.time import Time
        >>> from astropy import units as u
        >>> from sbpy.thermal import STM
        >>> from sbpy.data import Ephem, Phys
        >>> epoch = Time('2019-03-12 12:30:00', scale='utc')
        >>> eph = Ephem.from_horizons('2015 HW', location='568', epochs=epoch) # doctest: +REMOTE_DATA
        >>> phys = PhysProp('diam'=0.3*u.km, 'pv'=0.3) # doctest: +SKIP
        >>> lam = np.arange(1, 20, 5)*u.micron # doctest: +SKIP
        >>> flux = STM.flux(phys, eph, lam) # doctest: +SKIP

        not yet implemented

        """

    def fit(self, eph):
        """Fit thermal model to observations stored in `sbpy.data.Ephem` instance

        Parameters
        ----------
        eph : `sbpy.data.Ephem` instance, mandatory
            provide object ephemerides and flux measurements

        Examples
        --------
        >>> from sbpy.thermal import STM
        >>> stmfit = STM.fit(eph) # doctest: +SKIP

        not yet implemented

        """


class STM(ThermalClass):
    pass


class FRM(ThermalClass):
    pass


class NEATM(ThermalClass):
    def __init__(self):
        from .. import bib
        bib.register('sbpy.thermal.NEATM', {'method': '1998Icar..131..291H'})
