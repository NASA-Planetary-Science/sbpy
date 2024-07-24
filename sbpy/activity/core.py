# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Activity Core Module

Core module functions and classes, especially for handling coma
geometries.

created on June 23, 2017
"""


__all__ = [
    'Aperture',
    'CircularAperture',
    'AnnularAperture',
    'RectangularAperture',
    'GaussianAperture',
]


from abc import ABC, abstractmethod

import numpy as np
import astropy.units as u

from .. import data as sbd
from .. import units as sbu


class Aperture(ABC):
    """
    Abstract base class for photometric apertures.

    Notes
    -----
    The shape of the aperture must be passed as the first argument of
    `__init__`, or else `as_length` and `as_angle` must be overridden.

    """

    def __init__(self, dim):
        if not dim.unit.is_equivalent((u.radian, u.meter)):
            raise u.UnitTypeError(
                'aperture must be defined with angles or lengths.')

        self.dim = dim

    def __str__(self):
        """Description of the aperture."""
        # assumes preferred format for __repr__
        return repr(self)[1:-1].replace('Aperture:', ' aperture,')

    @abstractmethod
    def __repr__(self):
        """Preferred format <ShapedAperture: size>"""

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
    def as_angle(self, eph):
        """This aperture in units of angle.


        Parameters
        ----------
        eph : dictionary-like, `~sbpy.data.Ephem`, or `~astropy.units.Quantity`
            The observer-target distance (``delta``).


        Returns
        -------
        aper

        """

        dim = self.dim.to('arcsec', sbu.projected_size(eph))
        return type(self)(dim)

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
    def as_length(self, eph):
        """This aperture in units of length.


        Parameters
        ----------
        eph : dictionary-like, `~sbpy.data.Ephem`, or `~astropy.units.Quantity`
            The observer-target distance (``delta``).


        Returns
        -------
        aper

        """

        dim = self.dim.to('km', sbu.projected_size(eph))
        return type(self)(dim)

    @abstractmethod
    def coma_equivalent_radius(self):
        """Circular aperture radius that yields same flux for a 1/ρ coma.


        Returns
        -------
        rap : `~astropy.units.Quantity`

        """


class CircularAperture(Aperture):
    """Circular aperture projected at the distance of the target.


    Parameters
    ----------
    radius : `~astropy.units.Quantity`
        Angular or projected linear radius for the aperture.

    """

    def __init__(self, radius):
        super().__init__(radius)

    def __repr__(self):
        return '<CircularAperture: radius {}>'.format(self.dim)

    @property
    def radius(self):
        """Aperture radius."""
        return self.dim

    def coma_equivalent_radius(self):
        return self.radius
    coma_equivalent_radius.__doc__ = Aperture.coma_equivalent_radius.__doc__

    @classmethod
    def from_coma_equivalent(cls, aperture):
        """Initialize based on coma equivalent radius.


        Parameters
        ----------
        aperture : `Aperture` or `~sbpy.units.Quantity`
            Another aperture or a radius.

        """

        if isinstance(aperture, Aperture):
            radius = aperture.coma_equivalent_radius()
        else:
            radius = u.Quantity(aperture)

        return cls(radius)


class AnnularAperture(Aperture):
    """Annular aperture projected at the distance of the target.


    Parameters
    ----------
    shape : `~astropy.units.Quantity`
        A two-element `~astropy.units.Quantity` of angular or
        projected linear size for the inner and outer radius of the
        aperture.

    """

    def __init__(self, shape):
        if len(shape) != 2:
            raise ValueError('shape must be 2-elements')
        super().__init__(shape)

    def __repr__(self):
        return ('<AnnularAperture: radii {0[0].value:}–{0[1]:}>'
                .format(self.dim))

    @property
    def shape(self):
        """Annulus inner and outer radii."""
        return self.dim

    def coma_equivalent_radius(self):
        return max(self.dim) - min(self.dim)
    coma_equivalent_radius.__doc__ = Aperture.coma_equivalent_radius.__doc__


class RectangularAperture(Aperture):
    """Rectangular aperture projected at the distance of the target.

    Parameters
    ----------
    shape : `~astropy.units.Quantity`
        A two-element `~astropy.units.Quantity` of angular or
        projected linear size for the width and height of the
        aperture.  The order is not significant.

    """

    def __init__(self, shape):
        if len(shape) != 2:
            raise ValueError('shape must be 2-elements')
        super().__init__(shape)

    def __repr__(self):
        return ("<RectangularAperture: dimensions {0[0].value:}×{0[1]:}>"
                .format(self.dim))

    @property
    def shape(self):
        """Rectangle dimensions."""
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
    coma_equivalent_radius.__doc__ = Aperture.coma_equivalent_radius.__doc__


class GaussianAperture(Aperture):
    """Gaussian-shaped aperture, e.g., for radio observations.

    The aperture is normalized to 1.0 at the center.

    Parameters
    ----------
    sigma : `~astropy.units.Quantity`, optional
        The width of the Gaussian beam (square-root of the variance)
        as an angular or projected size.

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

    def __repr__(self):
        return "<GaussianAperture: 1-σ width {}>".format(self.dim)

    @property
    def sigma(self):
        """Beam Gaussian width."""
        return self.dim

    @property
    def fwhm(self):
        """Beam full-width at half-maximum."""
        return self.dim * 2.3548200450309493

    @sbd.dataclass_input
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
    def __call__(self, rho, eph: sbd.Ephem=None):
        """Evaluate the aperture.


        Parameters
        ----------
        rho : `~astropy.units.Quantity`
            Position to evaluate, in units of length or angle.

        eph : dictionary-like, `~sbpy.data.Ephem`, or `~astropy.units.Quantity`, optional
            The observer-target distance (``delta``).  Use ``eph`` to
            convert between angles and lengths, as needed.

        """

        if eph is not None:
            equiv = sbu.projected_size(eph)
        else:
            equiv = []
        x = rho.to(self.dim.unit, equiv)

        # normalize to 1.0 at the center
        return np.exp(-x**2 / self.sigma**2 / 2)

    def coma_equivalent_radius(self):
        # This beam is normalized to 1.0 at the center.
        return np.sqrt(np.pi / 2) * self.sigma
    coma_equivalent_radius.__doc__ = Aperture.coma_equivalent_radius.__doc__
