# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
========================
SBPy Spectroscopy Module
========================

created on June 23, 2017
"""

__all__ = ['Spectrum', 'SpectralModel']



class SpectralModel():
    """Range of spectral models"""
    def haser():
        """Haser model
        
        should allow direct creation of a `sbpy.actvity.Haser` instance"""
        pass

    def emission_lines():
        """Emission lines"""
        pass

    def reflectance():
        """Reflectance spectrum (asteroids)"""



class Spectrum():

    def __init__(self, flux, dispersionaxis, unit):
        self.flux = flux
        self.axis = dispersionaxis
        self.unit = unit


    @classmethod
    def read(cls, filename, columns='auto'):
        """Read spectrum from file
        
        Parameters
        ----------
        filename : str, mandatory
            data file name
        columns : str or list-like, optional, default: 'auto'
            file format, `auto` will try to recognize the format 
            automatically

        Returns
        -------
        `Spectrum` instance

        Examples
        --------
        >>> spec = Spectrum.read('2012_XY.dat')

        not yet implemented

        """

        
    def write(self, filename, columns='all'):
        """Write spectrum to file

        Parameters
        ----------
        filename : str, mandatory
            data file name
        columns : str or list-like, optional: default: 'all'
            file format; `all` will write all fields to the file

        Examples
        --------
        >>> spec = Spectrum.read('2012_XY.dat')
        >>> spec.write('2012_XY.dat.bak')

        not yet implemented

        """
            
    def convert_units(self, **kwargs):
        """Convert Spectrum units as provided by user

        Examples
        --------
        >>> spec.convert_units(flux_unit=u.K)
        >>> spec.convert_units(dispersion_unit=u.km/u.s)
        
        not yet implemented

        """


    def baseline(self, subtract=False):
        """fit baseline to `Spectrum` instance

        Parameters
        ----------
        subtract : bool, optional, default=False
            if `True`, subtract the baseline

        Returns
        -------
        float

        Examples
        --------
        >>> baseline = spec.baseline()
        >>> spec.baseline(subtract=True)

        not yet implemented

        """

    def slope(self, subtract=False):
        """fit slope to `Spectrum` instance

        Parameters
        ----------
        subtract : bool, optional, default=False
            if `True`, subtract the slope

        Returns
        -------
        float

        Examples
        --------
        >>> slope = spec.slope()
        >>> spec.slope(subtract=True)

        not yet implemented

        """

    def fit(self, spec):
        """Fit `SpectralModel` to different model types

        Parameters
        ----------
        spec : str, mandatory
            `SpectralModel` instance to fit

        Examples
        --------
        >>> spec_model = SpectralModel(type='Haser', molecule='H2O')
        >>> spec.fit(spec_model)
        >>> print(spec.fit_info)

        not yet implemented 
        
        """

    def plot(self):
        """Plot spectrum 
       
        Returns
        -------
        `matplotlib.pyplot` instance

        not yet implemented
        """
