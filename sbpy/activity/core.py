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
import numpy as np
import astropy.units as u

class Aperture(ABC):
    """Abstract base class for photometric apertures."""

    @abstractmethod
    def coma_equivalent_radius(self):
        """Circular aperture radius that yields same flux for a 1/ρ coma.

        Returns
        -------
        rap : `~astropy.units.Quantity`

        """

        pass

class CircularAperture(Aperture):
    """A circular aperture projected at the distance of the target.

    Parameters
    ----------
    radius : `~astropy.units.Quantity`
      Angular or projected linear radius for the aperture.

    """

    def __init__(self, radius):
        assert isinstance(radius, u.Quantity)
        assert radius.unit.is_equivalent((u.radian, u.meter))
        self.radius = radius

    def coma_equivalent_radius(self):
        return self.radius

class RectangularAperture(Aperture):
    """A rectangular aperture projected at the distance of the target.
    
    Parameters
    ----------
    shape : `~astropy.units.Quantity`
      A two-element `~astropy.units.Quantity` of angular or projected
      linear size for the aperture.

    """

    def __init__(self, shape):
        assert isinstance(shape, u.Quantity)
        assert shape.unit.is_equivalent((u.radian, u.meter))
        assert len(shape) == 2
        self.shape = shape

    def coma_equivalent_radius(self):
        # Int_0^θ Int_0^r 1/r * r * dr dθ
        # limits on r depend on θ; eq. of a line: x * sec(θ)
        # --> Int_0^θ Int_0^x*sec(θ) dr dθ
        # --> Int_0^θ x*sec(θ) dθ
        # --> x * log(tan(θ) + sec(θ))

        # First, integrate the 1/rho distribution over the first
        # octant of the rectangle in polar coordinates.  The azimuthal
        # limits are 0 to arctan(y / x).

        # th = np.arctan(y / x)
        # I = (x / 2) * log(tan(th) + sec(th))
        # sec(th) = cos(th)**-1
        # cos(arctan(y / x)) = 1 / sqrt(1 + (y / x)**2)
        # sec(th) = sqrt(1 + (y / x)**2)
        # I1 = x / 2 * np.log(y / x + np.sqrt(1 + (y / x)**2))

        # Then, integrate the second octant: th = 0 to arctan(y / x)
        # I2 = y / 2 * np.log(x / y + np.sqrt(1 + (x / y)**2))

        # The two octants correspond to 1/4 the full area:
        # I = 4 * (I1 + I2)

        # For the circular aperture, the integral is 2 pi rho.
        # rho = I / 2 / np.pi

        # implement the above, moving constants around
        x, y = self.shape
        I1 = x * np.log(y / x + np.sqrt(1 + (y / x)**2))
        I2 = y * np.log(x / y + np.sqrt(1 + (x / y)**2))
        return (I1 + I2) / np.pi

class GaussianAperture(Aperture):
    """A Gaussian-shaped aperture, typically used for radio observations.

    Parameters
    ----------
    sigma : `~astropy.units.Quantity`, optional
      The width of the Gaussian beam (square-root of the variance) as
      an angular or projected size.
    fwhm : `~astropy.units.Quantity`, optional
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
        assert isinstance(s, u.Quantity)
        assert s.unit.is_equivalent((u.radian, u.meter))
        self._sigma = s

    @property
    def fwhm(self):
        return self._sigma * 2.3548200450309493

    @fwhm.setter
    def fwhm(self, f):
        self.sigma = f / 2.3548200450309493

    def coma_equivalent_aperture(self):
        # This beam is normalized to 1.0 at the center, is that OK?
        return np.sqrt(np.pi / 2) * self.sigma
