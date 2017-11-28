# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
SBPy Activity Module
====================

created on June 23, 2017
"""

import astropy.units as u

__all__ = ['Species', 'Aperture', 'ActivityClass', 'Haser',
           'Vectorial', 'Syndynes']


class Species():

    def scale_length(species, eph):
        """Return scale length of a specific species or list thereof at a
        given heliocentric distance (provided by eph)
        """

        # implement list treatment

        data = {'H2O': 2.4e4, #Cochran and Schleicher 1993
                'OH': 1.6e5 #Cochran and Schleicher 1993
        }
        
        try:
            out = data[species]
        except KeyError:
            raise KeyError(("no scale length available for species "
                            "'{:s}'").format(species))
        return out*u.km/eph.rh**2


    def fluorescence_length(species, eph):
        """Return fluorescence band efficiencies of a specific species and
        transition or list thereof at a given heliocentric distance
        (provided by eph)

        """

        # implement list treatment

        data = {'OH 0-0': 2e-15 #Schleicher and A'Hearn 1988
        }
        
        try:
            out = data[species]
        except KeyError:
            raise KeyError(("no scale length available for species "
                            "'{:s}'").format(species))
        return out* u.erg / u.s / u.molecule * u.au**2/eph.rh**2



class Aperture():

    def circular(r):
        """Defines a circular aperture"""

    def rectangular(a, b):
        """Defines a rectangular aperture"""

    def gaussian(r):
        """Defines a Gaussian aperture"""


        

class ActivityClass():
    """General model implementation"""
    
    def __init__(self, Q):
        self. Q = Q

    def volume_density(self, aperture, eph):
        """Calculate coma volume density in aperture

        Parameters
        ----------
        aperture : `sbpy.activity.Aperture` instance, mandatory
            aperture
        eph : `sbpy.data.Ephem` instance, mandatory
            ephemerides at epoch, required: heliocentric distance `rh`

        Returns
        -------
        integer

        Examples
        --------
        TBD

        not yet implemented

        """



    def column_density(self, aperture, eph):
        """Calculate coma column density in aperture

        Parameters
        ----------
        aperture : `sbpy.activity.Aperture` instance, mandatory
            aperture
        eph : `sbpy.data.Ephem` instance, mandatory
            ephemerides at epoch, required: heliocentric distance `rh`

        Returns
        -------
        integer

        Examples
        --------
        TBD

        not yet implemented

        """
        
    def total_number(self, aperture, eph):
        """Calculate total number of molecules in aperture 

        Parameters
        ----------
        aperture : `sbpy.activity.Aperture` instance, mandatory
            aperture
        eph : `sbpy.data.Ephem` instance, mandatory
            ephemerides at epoch, required: heliocentric distance `rh`

        Returns
        -------
        integer

        Examples
        --------
        TBD

        not yet implemented

        """
        

class Haser(ActivityClass):
    """Haser model implementation"""
    
    def __init__(self, Q, gamma):
        """Parameters
        ----------
        Q : `Astropy.units` quantity or iterable, mandatory
            production rate usually in units of `u.molecule / u.s`
        gamma : `Astropy.units` quantity or iterable, mandatory
            scale length usually in units of `u.km`

        Returns
        -------
        Haser instance

        Examples
        --------
        TBD

        not yet implemented

        """
        
        
        


class Vectorial(ActivityClass):
    """Vectorial model implementation"""
    
    def __init__(self, Q, species):
        """Parameters
        ----------
        Q : `Astropy.units` quantity or iterable, mandatory
            production rate usually in units of `u.molecule / u.s`
        species : dictionary or list of dictionares, mandatory
            defines gas velocity, lifetimes, disassociative lifetimes 

        Returns
        -------
        Vectorial instance

        Examples
        --------
        TBD

        not yet implemented

        """



class Syndynes():
    """Syndynes and Synchrones"""

    def __init__(orb, date):
        self.orb = orb
        self.date = date

    def plot_syndynes(self, betas, ages, location='500'):
        """Parameters
        ----------
        betas : array, mandatory
            beta values
        ages : array, mandatory
            synchrone ages
        location : str, optional, default: '500' (geocentric)
            observer location MPC code

        Returns
        -------
        `matplotlib.pyplot` instance

        Examples
        --------
        TBD

        not yet implemented

        """
        
