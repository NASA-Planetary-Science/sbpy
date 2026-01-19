# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Asteroid and comet nucleus shape models
"""

__all__ = ["Shape"]

import astropy.units as u

from ..calib import Sun
from ..data.ephem import Ephem
from ..spectroscopy.sources import SinglePointSpectrumError


class Shape:
    """Model asteroid or comet nucleus shape."""

    @staticmethod
    def _incident_sunlight(
        wfb: u.Quantity | str, eph: Ephem, unit: u.Unit, interpolate: bool
    ) -> u.Quantity:
        """Calculate incident sunlight for a normal surface."""

        sun = Sun.from_default()
        try:
            S = sun.observe(wfb, unit=unit, interpolate=interpolate)
        except SinglePointSpectrumError:
            S = sun(wfb, unit=unit)
        S /= eph["rh"][0].to_value("au") ** 2

        return S


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
