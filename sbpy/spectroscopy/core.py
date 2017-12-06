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
        >>> spec = Spectrum.read('2012_XY.dat') # doctest: +SKIP

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
        >>> spec = Spectrum.read('2012_XY.dat') # doctest: +SKIP
        >>> spec.write('2012_XY.dat.bak') # doctest: +SKIP

        not yet implemented

        """
            
    def convert_units(self, **kwargs):
        """Convert Spectrum units as provided by user

        Examples
        --------
        >>> spec.convert_units(flux_unit=u.K) # doctest: +SKIP
        >>> spec.convert_units(dispersion_unit=u.km/u.s) # doctest: +SKIP
        
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
        >>> baseline = spec.baseline() # doctest: +SKIP
        >>> spec.baseline(subtract=True) # doctest: +SKIP

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
        >>> slope = spec.slope() # doctest: +SKIP
        >>> spec.slope(subtract=True) # doctest: +SKIP

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
        >>> spec_model = SpectralModel(type='Haser', molecule='H2O') # doctest: +SKIP
        >>> spec.fit(spec_model) # doctest: +SKIP
        >>> print(spec.fit_info) # doctest: +SKIP

        not yet implemented 
        
        """

    def plot(self):
        """Plot spectrum 
       
        Returns
        -------
        `matplotlib.pyplot` instance

        not yet implemented
        """
