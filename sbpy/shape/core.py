# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Asteroid and comet nucleus shape models
"""

from typing import Callable
import astropy.units as u

from ..data.ephem import Ephem
from ..data.decorators import dataclass_input


class Shape:
    """Model asteroid or comet nucleus shape."""


class Sphere(Shape):
    """A spherical object.


    Parameters
    ----------


    """

    def __init__(self, radius: u.physical.length):
        self.radius = radius

    def to_faceted_model(self):
        pass

    @dataclass_input
    def integrate_over_surface(
        self, func: Callable, eph: Ephem, *args, **kwargs
    ) -> u.Quantity:
        """Integrate the function over the surface.

        The integration is over the observed area::




        Parameters
        ----------
        func : callable
            The function to integrate: ``func(*args, i, e, phi, **kwargs)``,
            where :math:`i` is the angle of incidence, :math:`e` is the the
            angle of emittance, and :math:`phi` is the phase angle.

        eph : `Ephem`
            The observing geometry as an ephemeris (or equivalent) object:
                - rh
                - delta
                - phase

        *args, **kwargs
            Additional arguments passed to the function.


        Returns
        -------
        total : `~sbpy.units.Quantity`

        """


# __all__ = ["ModelClass", "Kaasalainen", "Lightcurve"]


# class ModelClass:

#     def __init__(self):
#         self.shape = None

#     def load_obj(self, filename):
#         """Load .OBJ shape model"""

#     def get_facet(self, facetid):
#         """Retrieve information for one specific surface facet"""

#     def iof(self, facetid, eph):
#         """Derive I/F for one specific surface facet"""

#     def integrated_flux(self, eph):
#         """Derive surface-integrated flux"""

#     def lightcurve(self, eph, time):
#         """Derive lightcurve"""


# class Kaasalainen(ModelClass):

#     def __init__(self):
#         self.properties = None
#         # setup model properties

#     def convexinv(self, eph):
#         """Obtain shape model through lightcurve inversion, optimizing all
#         parameters and uses spherical harmonics functions for shape
#         representation"""

#     def conjgradinv(self, eph):
#         """Obtain shape model through lightcurve inversion, optimizing only
#         shape and uses directly facet areas as parameters"""


# class Lightcurve:

#     def __init__(self, eph):
#         self.eph = eph
#         self.fouriercoeff = None
#         self.period = None
#         self.pole = (0, 90)  # ecliptic coordinates

#     def axis_ratio(self):
#         """Derive axis ratio from lightcurve amplitude"""

#     def derive_period(self, method="lomb-scargle"):
#         """Derive lightcurve period using different methods"""

#     def fit_fourier(self, order):
#         """Fit Fourier coefficients to lightcurve"""

#     def fit_pole(self):
#         """Fit pole orientation"""

#     def fit(self):
#         """Fit period, pole orientation, and Fourier coefficients at the same time"""

#     def simulate(self):
#         """Simulate a lightcurve from period, Fourier coefficients, pole orientation"""
