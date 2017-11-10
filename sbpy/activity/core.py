# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=========================
SBPy Activity Core Module
=========================

Core module functions and classes, especially for handling coma
geometries.


Functions
---------


Classes
-------
Aperture            - Abstract base class for observation apertures.
CircularAperture    - Circular photometric aperture for observing
                      a coma model.
RectangularAperture - Rectangular photometric aperture for observing a
                      coma model.
GaussianAperture    - Gaussian photometric aperture for observing a coma
                      model.


created on June 23, 2017

"""

__all__ = [
    'CircularAperture',
    'RectangularAperture',
    'GaussianAperture',
]

from abc import ABC, abstractmethod
import astropy.units as u
from astropy.units.quantity import Quantity

class Aperture(ABC):
    """Abstract base class for photometric apertures."""

    @abstractmethod
    def coma_equivalent_radius(self):
        """Equivalent circular aperture radius for 1/ρ coma.

        Returns
        -------
        rap : `~Quantity`

        """

        pass

class CircularAperture(Aperture):
    """A circular aperture projected at the distance of the target.

    Parameters
    ----------
    radius : astropy Quantity
      Angular or projected linear radius for the aperture.

    """

    def __init__(self, radius):
        assert isinstance(radius, Quantity)
        assert radius.unit.is_equivalent((u.radian, u.meter))
        self.radius = radius

    def coma_equivalent_radius(self):
        return self.radius

class RectangularAperture(Aperture):
    """A rectangular aperture projected at the distance of the target.
    
    Parameters
    ----------
    shape : astropy Quantity
      A two-element `Quantity` of angular or projected linear size for
      the aperture.

    """

    def __init__(self, shape):
        assert isinstance(shape, Quantity)
        assert shape.unit.is_equivalent((u.radian, u.meter))
        assert len(shape) == 2
        self.shape = shape

class GaussianAperture(Aperture):
    """A Gaussian-shaped aperture, typically used for radio observations.

    Parameters
    ----------
    sigma : astropy Quantity, optional
      The width of the Gaussian beam (square-root of the variance) as
      an angular or projected size.
    fwhm : astropy Quantity, optional
      The full-width at half-maximum of the Gaussian beam as an
      angular or projected size.

    Notes
    -----
    One of `sigma` or `fwhm` is required.

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
        assert isinstance(s, Quantity)
        assert s.unit.is_equivalent((u.radian, u.meter))
        self._sigma = s

    @property
    def fwhm(self):
        return self._sigma * 2.3548200450309493

    @fwhm.setter
    def fwhm(self, f):
        self.sigma = f / 2.3548200450309493

