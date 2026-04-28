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

    @staticmethod
    def _planck_wave(wave, T):
        """Planck function for spectral radiance in per unit wavelength.


        Parameters
        ----------
        wave : float or ndarray
            Wavelength in meters.

        T : float or ndarray
            Temperature in Kelvin.


        Returns
        -------
        B : float or ndarray
            W / (m3 sr).

        """

        # 2 h c2
        c1 = 1.1910429723971884e-16  # J m2 / s

        # h c / k_B
        c2 = 0.014387768775039337  # K m

        a = np.exp(c2 / (wave * T))
        B = c1 / (wave**5 * (a - 1))

        return B

    @staticmethod
    def _planck_freq(freq, T):
        """Planck function for spectral radiance in per unit frequency.


        Parameters
        ----------
        freq : float or ndarray
            Frequency in Hz.

        T : float or ndarray
            Temperature in Kelvin.


        Returns
        -------
        B : float or ndarray
            W / (m2 Hz sr) for frequencies, W / (m3 sr) for wavelengths.

        """

        # 2 h / c2
        c1 = 1.4744994647625417e-50  # J s3 / m2

        # h / k_B
        c2 = 4.799243073366221e-11  # K s

        a = np.exp(c2 * freq / T)
        B = c1 * freq**3 / (a - 1)

        return B

    # @u.quantity_input
    def _emission(
        self,
        wave_freq: u.Quantity["length"] | u.Quantity["frequency"],
        i: u.Quantity["angle"],
        e: u.Quantity["angle"],
        phi: u.Quantity["angle"],
        # unit: u.Unit | None = None,
    ):
        """Same as emission, but without quantities."""

        T = self.T(i).value

        # emissivity = 1; it was already accounted for in the temperature
        e = self.surface.emittance(1, e, phi).value

        if wave_freq.unit.is_equivalent(u.m):
            B = self._planck_wave(wave_freq.si.value, T)
        elif wave_freq.unit.is_equivalent(u.Hz):
            B = self._planck_freq(wave_freq.si.value, T)
        I_e = B * e

        return I_e

    # @u.quantity_input
    def emission(
        self,
        wave_freq: u.Quantity["length"] | u.Quantity["frequency"],
        i: u.Quantity["angle"],
        e: u.Quantity["angle"],
        phi: u.Quantity["angle"],
        unit: u.Unit | None = None,
    ):
        """Observed spectral radiance of a surface element."""

        unit = u.W / u.m**2 / u.sr / wave_freq.unit if unit is None else u.Unit(unit)

        I_e = self._emission(wave_freq, i, e, phi)

        return I_e.to(unit, u.spectral_density(wave_freq))


class InstantaneousEquilibrium(ThermalEmission):
    @property
    def T0(self) -> u.Quantity[u.K]:
        """Temperature of a surface normal to the incident radiation."""
        if self._update_T0:
            self._a0 = self.surface.absorptance(1, 0 * u.deg)
            I_a = self.I * (1 - self.albedo) * self._a0
            self._T0 = (
                (I_a / (pi * self.epsilon * self.eta * const.sigma_sb)) ** (1 / 4)
            ).decompose()
            self._update_T0 = False

        return self._T0

    # @u.quantity_input
    def T(self, i: u.Quantity["angle"]) -> u.Quantity[u.K]:
        # emissivity = 1; it was already accounted when T0 was calculated
        return self.T0 * (self.surface.absorptance(1, i) / self._a0) ** (1 / 4)
