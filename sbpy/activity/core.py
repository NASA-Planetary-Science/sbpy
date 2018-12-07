# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Activity Core Module

Core module functions and classes, especially for handling coma
geometries.

created on June 23, 2017
"""

__all__ = [
    'rho_as_angle',
    'rho_as_length',
    'CircularAperture',
    'AnnularAperture',
    'RectangularAperture',
    'GaussianAperture',
]

from abc import ABC, abstractmethod
import numpy as np
import astropy.units as u


def rho_as_angle(rho, eph):
    """Projected linear distance to angular distance.

    Parameters
    ----------
    rho : `~astropy.units.Quantity`
      Projected distance in units of length.

    eph : dictionary-like or `~sbpy.data.Ephem`
      Ephemerides; requires geocentric distance as `delta`.

    Returns
    -------
    rho_l : `~astropy.units.Quantity`

    """

    if rho.unit.is_equivalent(u.m):
        rho_a = np.arctan(rho / eph['delta'].to(u.m))
    elif rho.unit.is_equivalent(u.rad):
        rho_a = rho
    else:
        raise ValueError('rho must have units of length or angle')

    return rho_a


def rho_as_length(rho, eph):
    """Angular distance to projected linear distance.

    Parameters
    ----------
    rho : `~astropy.units.Quantity`
      Projected distance in units of angle.

    eph : dictionary-like or `~sbpy.data.Ephem`
      Ephemerides; requires geocentric distance as `delta`.

    Returns
    -------
    rho_l : `~astropy.units.Quantity`

    """

    if rho.unit.is_equivalent(u.rad):
        rho_l = eph['delta'].to(u.m) * np.tan(rho)
    elif rho.unit.is_equivalent(u.m):
        rho_l = rho
    else:
        raise ValueError('rho must have units of length or angle.')

    return rho_l


class Aperture(ABC):
    """Abstract base class for photometric apertures.

    The shape of the aperture must be passed as the first argument of
    `__init__`, or else `as_length` and `as_angle` must be overridden.

    """

    def __init__(self, dim):
        if not dim.unit.is_equivalent((u.radian, u.meter)):
            raise ValueError(
                'aperture must be defined with angles or lengths.')

        self.dim = dim

    @abstractmethod
    def __str__(self):
        return "Abstract aperture of size {}".format(self.dim)

    def as_angle(self, eph):
        """This aperture in units of angle.

        Parameters
        ----------
        eph : dictionary-like or `~sbpy.data.Ephem`, optional
          Ephemerides at epoch; requires geocentric distance as
          `delta` keyword.  Ignored if the aperture is already in
          units of angle.

        Returns
        -------
        aper

        """

        dim = rho_as_angle(self.dim, eph)
        return type(self)(dim)

    def as_length(self, eph):
        """This aperture in units of length.

        Parameters
        ----------
        eph : dictionary-like or `~sbpy.data.Ephem`, optional
          Ephemerides at epoch; requires geocentric distance as
          `delta` keyword.  Ignored if the aperture is already in
          units of length.

        Returns
        -------
        aper

        """

        dim = rho_as_length(self.dim, eph)
        return type(self)(dim)

    def _convert_unit(self, rho, eph):
        """Make units match those of self."""
        if not self.dim.unit.is_equivalent(rho.unit):
            if rho.unit.is_equivalent(u.m):
                x = rho_as_angle(rho, eph)
            else:
                x = rho_as_length(rho, eph)
        else:
            x = rho

        return x

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
        super().__init__(radius)

    def __str__(self):
        return "Circular aperture, radius {}".format(self.dim)

    @property
    def radius(self):
        return self.dim

    def coma_equivalent_radius(self):
        return self.radius


class AnnularAperture(Aperture):
    """Annular aperture projected at the distance of the target.

    Parameters
    ----------
    shape : `~astropy.units.Quantity`
      A two-element `~astropy.units.Quantity` of angular or projected
      linear size for the inner and outer radius of the aperture.

    """

    def __init__(self, shape):
        if len(shape) != 2:
            raise ValueError('shape must be 2-elements')
        super().__init__(shape)

    def __str__(self):
        return "Annular aperture, radii {0[0].value:}–{0[1]:}".format(self.dim)

    @property
    def shape(self):
        return self.dim

    def coma_equivalent_radius(self):
        return max(self.dim) - min(self.dim)


class RectangularAperture(Aperture):
    """Rectangular aperture projected at the distance of the target.

    Parameters
    ----------
    shape : `~astropy.units.Quantity`
      A two-element `~astropy.units.Quantity` of angular or projected
      linear size for the width and height of the aperture.  The order
      is not significant.

    """

    def __init__(self, shape):
        if len(shape) != 2:
            raise ValueError('shape must be 2-elements')
        super().__init__(shape)

    def __str__(self):
        return "Rectangular aperture, dimensions {0[0].value:}×{0[1]:}".format(self.dim)

    @property
    def shape(self):
        return self.dim

    def coma_equivalent_radius(self):
        # Int_0^θ Int_0^r 1/r * r * dr dθ
        # limits on r depend on θ; eq. of a line: x * sec(θ)
        # --> Int_0^θ Int_0^x*sec(θ) dr dθ
        # --> Int_0^θ x*sec(θ) dθ
        # --> x * log(tan(θ) + sec(θ))

        # First, integrate the 1/rho distribution over the first
        # "octant" of the rectangle in polar coordinates.  The
        # azimuthal limits are 0 to arctan(y / x).  Also, x, and y are
        # the full rectangle dimensions, so they must be halved.

        # th = np.arctan(y / x)
        # I = (x / 2) * log(tan(th) + sec(th))
        # sec(th) = cos(th)**-1
        # cos(arctan(y / x)) = 1 / sqrt(1 + (y / x)**2)
        # sec(th) = sqrt(1 + (y / x)**2)
        # I1 = x / 2 * np.log(y / x + np.sqrt(1 + (y / x)**2))

        # Then, integrate the second "octant": th = 0 to arctan(x / y)
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
    """A Gaussian-shaped aperture, e.g., for radio observations.

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
        if (sigma is None) and (fwhm is None):
            raise ValueError('One of `sigma` or `fwhm` must be defined')

        if sigma is not None:
            super().__init__(sigma)
        else:
            super().__init__(fwhm / 2.3548200450309493)

    def __str__(self):
        return "Gaussian aperture, 1-σ width {}".format(self.dim)

    @property
    def sigma(self):
        return self.dim

    @property
    def fwhm(self):
        return self.dim * 2.3548200450309493

    def __call__(self, rho, eph=None):
        """Evaluate the aperture.

        Parameters
        ----------
        rho : `~astropy.units.Quantity`
          Position to evaluate, in units of length or angle.
        eph : dictionary-like or `~sbpy.data.Ephem`, optional
          Ephemerides at epoch; requires geocentric distance as
          `delta` keyword if the aperture's units and `rho`'s units do
          not match.

        """

        x = self._convert_unit(rho, eph)

        # normalize to 1.0 at the center?
        return np.exp(-x**2 / self.sigma**2 / 2)

    def coma_equivalent_radius(self):
        # This beam is normalized to 1.0 at the center, is that correct?
        return np.sqrt(np.pi / 2) * self.sigma
