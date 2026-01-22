# Licensed under a 3-clause BSD style license - see LICENSE.rst

import abc
import numpy as np
from numpy import pi
import astropy.units as u
import astropy.constants as const
from astropy.modeling.models import BlackBody
from .surface import Surface


class ThermalEmission(abc.ABC):
    """Abstract base class for calculating surface thermal emission."""

    def __init__(self, I, albedo, epsilon, eta, surface: Surface):
        self.I = I
        self.albedo = albedo
        self.epsilon = epsilon
        self.eta = eta
        self.surface: Surface = surface
        self._update_T0 = True

    @property
    def I(self) -> u.Quantity["W/m2"]:
        return self._I

    @I.setter
    @u.quantity_input
    def I(self, i: u.Quantity["W/m2"]):
        self._I = i
        self._update_T0 = True

    @property
    def albedo(self) -> u.Quantity[""]:
        return self._albedo

    @albedo.setter
    @u.quantity_input
    def albedo(self, a: u.Quantity[""]):
        self._albedo = a
        self._update_T0 = True

    @property
    def epsilon(self) -> u.Quantity[""]:
        return self._epsilon

    @epsilon.setter
    @u.quantity_input
    def epsilon(self, e: u.Quantity[""]):
        self._epsilon = e
        self._update_T0 = True

    @property
    def eta(self) -> u.Quantity[""]:
        return self._eta

    @eta.setter
    @u.quantity_input
    def eta(self, e: u.Quantity[""]):
        self._eta = e
        self._update_T0 = True

    @property
    def surface(self) -> Surface:
        return self._surface

    @surface.setter
    @u.quantity_input
    def surface(self, s: Surface):
        self._surface = s
        self._update_T0 = True

    @abc.abstractmethod
    def T0(self) -> u.Quantity[u.K]:
        """Temperature of a surface normal to the incident radiation."""
        pass

    @property
    @abc.abstractmethod
    def T(self, i: u.Quantity["angle"]) -> u.Quantity[u.K]:
        """Temperature of a surface at angle ``i`` with incident radiation."""
        pass

    def emission(self, wave_freq, i, e, phi, unit=None):
        """Observed spectral radiance of a surface element."""

        if unit is None:
            unit = u.W / u.m**2 / u.sr / wave_freq.unit

        T = self.T(i)

        if u.isclose(T, 0 * u.K):
            return 0 * unit

        B = BlackBody(T)(wave_freq)
        # emissivity = 1; it was already accounted for in the temperature
        I_e = B * self.surface.emittance(1, e, phi)

        return I_e.to(unit, u.spectral_density(wave_freq))


class InstantaneousEquilibrium(ThermalEmission):
    @property
    def T0(self) -> u.Quantity[u.K]:
        """Temperature of a surface normal to the incident radiation."""
        if self._update_T0:
            I_a = self.I * self.surface.absorptance(1 - self.albedo, 0 * u.deg)
            self._T0 = (
                (I_a / (pi * self.epsilon * self.eta * const.sigma_sb)) ** (1 / 4)
            ).decompose()
            self._update_T0 = False

        return self._T0

    @u.quantity_input
    def T(self, i: u.Quantity["angle"]) -> u.Quantity[u.K]:
        # emissivity = 1; it was already accounted when T0 was calculated
        return self.T0 * self.surface.absorptance(1, i) ** (1 / 4)
