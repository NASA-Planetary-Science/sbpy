# Licensed under a 3-clause BSD style license - see LICENSE.rst

from abc import abstractmethod

import numpy as np
from numpy import pi
from astropy.coordinates import Angle
import astropy.constants as const
import astropy.units as u

from .surface import Surface
from ..spectroscopy.sources import BlackbodySource
from ..units.typing import SpectralQuantity, SpectralRadianceQuantity


class SurfaceThermalEmission(Surface):
    """Abstract base class for surface thermal emission.

    The surface is assumed to be illuminated by sunlight.


    Parameters
    ----------
    phys : Phys
        Surface physical parameters:
        - albedo
        - emissivity
        - beaming parameter

    """

    @abstractmethod
    def T(self, i: Angle, rh: u.physical.length) -> u.Quantity:
        """Surface temperature given incidence angle and heliocentric distance.


        Parameters
        ----------
        i : Angle
            Angle from normal of incident light.

        rh : `~astropy.units.Quantity`
            Heliocentric distance.


        Returns
        -------
        T : `~astropy.units.Quantity`
            Surface temperature.

        """

    def radiance(self, wave_freq: SpectralQuantity, n, rs, ro) -> u.Quantity:
        """Observed radiance from a surface.


        Parameters
        ----------
        wave_freq : `~astropy.units.Quantity`
            Wavelength or frequency at which to evaluate the emission.

        n : `numpy.ndarray`
            Surface normal vector.

        rs : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `~astropy.units.Quantity`
            Radial vector from the surface to the observer.


        Returns
        -------
        radiance : `~astropy.units.Quantity`
            Observed radiance.

        """

        n_hat = np.linalg.norm(n.value)
        rs_hat = np.linalg.norm(rs.value)
        ro_hat = np.linalg.norm(ro.value)

        i = np.arccos(np.dot(n_hat, rs_hat))
        e = np.arccos(np.dot(n_hat, ro_hat))
        phi = np.arccos(np.dot(rs_hat, ro_hat))

        rh = np.sqrt(np.dot(rs, rs))
        T = self.T(i, rh)

        epsilon = self.phys["emissivity"]
        BB = BlackbodySource(T)(wave_freq) / pi

        return epsilon * BB * self.emittance(e, phi)


class InstantaneousThermalEquilibrium(SurfaceThermalEmission):
    """Abstract base class for a surface in LTE with sunlight."""

    def T(self, i: Angle, rh: u.physical.length) -> u.Quantity:
        sun = const.L / (4 * pi * rh**2)
        epsilon = self.phys["emissivity"]
        eta = self.phys["eta"]
        T = (self.absorptance(i) * sun / epsilon / eta / const.sigma_sb) ** (1 / 4)
        return T
