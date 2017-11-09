# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
SBPy Activity Module
====================

created on June 23, 2017
"""

__all__ = ['Species', 'Aperture', 'CircularAperture', 'RectangularAperture',
           'GaussianAperture', 'ActivityClass', 'Haser',
           'Vectorial', 'Syndynes']

from abc import ABC, abstractmethod
import astropy.units as u

def photo_length_scale(species, rdot=None, eph=None, source=None):
    """Photodissociation length scale for a gas species.

    The length scale returned is scaled to 1 au.

    Parameters
    ----------
    species : string
      The species to look up.
    rdot : astropy Quantity, optional
      The helocentric radial velocity.  One of `rdot` or `eph` is
      required for some species.
    eph : sbpy Ephem, optional
      The target ephemeris.  Must include heliocentric radial
      velocity.  One of `rdot` or `eph` is required for some species.
    source : string, optional
      Use values from this source.

    Returns
    -------
    length_scale : astropy Quantity
      The length scale quoted at 1 au.

    Example
    -------
    >>> from sbpy.activity import photo_length_scale
    >>> gamma = photo_length_scale('OH', rdot=23 * u.km / u.s)

    not yet implemented

    """
    
    data = {'H2O': 2.4e4, #Cochran and Schleicher 1993
            'OH': 1.6e5 #Cochran and Schleicher 1993
    }
   
    try:
        out = data[species]
    except KeyError:
        raise KeyError(("no scale length available for species "
                        "'{:s}'").format(species))
    return out*u.km/eph.rh**2

def fluorescence_band_strength(species):
    """Fluorescence band efficiencies of a specific species and
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

class Aperture(ABC):
    """Abstract base class for photometric apertures."""
    pass

class CircularAperture(Aperture):
    """A circular aperture projected at the distance of the target.

    Parameters
    ----------
    radius : astropy Quantity
      Angular or projected linear radius for the aperture.

    """

    def __init__(self, radius):
        assert isinstance(radius, u.Quantity)
        assert radius.unit.is_equivalent((u.radian, u.meter))
        self.radius = radius

class RectangularAperture(Aperture):
    """A rectangular aperture projected at the distance of the target.
    
    Parameters
    ----------
    shape : astropy Quantity
      A two-element `Quantity` of angular or projected linear size for
      the aperture.

    """

    def __init__(self, shape):
        assert isinstance(shape, u.Quantity)
        assert shape.unit.is_equivalent((u.radian, u.meter))
        self.shape = shape

class GaussianAperture(Aperture):
    """A Gaussian-shaped aperture, typically used for radio observations.

    Parameters
    ----------
    sigma : astropy Quantity, optional
      The width of the Gaussian beam (square-root of the variance) as
      an angular or projected size.  One of `sigma` or `fwhm` is
      required.
    fwhm : astropy Quantity, optional
      The full-width at half-maximum of the Gaussian beam as an
      angular or projected size.  One of `sigma` or `fwhm` is
      required.

    """

    def __init__(self, sigma=None, fwhm=None):
        assert (sigma is not None) or (fwhm is not None), "One of `sigma` or `fwhm` must be defined."
        # units are tested in self.sigma
        if sigma is not None:
            self.sigma = sigma
        else:
            self.fwhm = fwhm

    @property
    def sigma(self):
        return self._sigma

    @sigma.setter
    def sigma(self, s):
        assert s.is_equivalent((u.radian, u.meter))
        self._sigma = s

    @property
    def fwhm(self):
        return self._sigma * 2.3548200450309493

    @fwhm.setter
    def fwhm(self, f):
        self.sigma = f / 2.3548200450309493

class Activity(ABC):
    """Abstract base class for gas coma models."""
    
    def __init__(self, Q):
        assert Q.unit.is_equivalent((u.s**-1, u.mol / u.s))
        self.Q = Q

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
        
