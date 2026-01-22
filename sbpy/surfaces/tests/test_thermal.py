# Licensed under a 3-clause BSD style license - see LICENSE.rst

from numpy import pi
from astropy import units as u
from astropy.modeling.models import BlackBody

from ..thermal import InstantaneousEquilibrium, ThermalEmission
from ..lambertian import LambertianSurface


class ConstantTemperature(ThermalEmission):
    def __init__(self, T0):
        super().__init__(0 * u.W / u.m**2, 0, 1, 1, LambertianSurface())
        self._T0 = T0

    @property
    def T0(self):
        return self._T0

    def T(self, i):
        return self._T0


class TestThermalEmission:
    def test_emission(self):
        bb200 = ConstantTemperature(200 * u.K)
        wave = 10 * u.um
        L = bb200.emission(wave, 0 * u.deg, 0 * u.deg, 0 * u.deg)
        L200 = BlackBody(200 * u.K)(wave).to(L.unit, u.spectral_density(wave))
        assert u.isclose(L, L200)

        L306060 = bb200.emission(wave, 30 * u.deg, 60 * u.deg, 60 * u.deg)
        assert u.isclose(L306060, L200 / 2)

        bb0 = ConstantTemperature(0 * u.K)
        wave = 10 * u.um
        L = bb0.emission(wave, 0 * u.deg, 0 * u.deg, 0 * u.deg)
        assert u.isclose(L, 0 * L.unit)


class TestInstantaneousEquilibrium:
    def test_T0(self):
        S = 1361.16646541 * u.W / u.m**2
        thermal = InstantaneousEquilibrium(S, 0, 1, 1, LambertianSurface())
        T = (S.value / (pi * 5.6703744191844314e-08)) ** (1 / 4) * u.K

        assert u.isclose(thermal.T0, T)

        # only albedo and eta affect temperature (emissivity is for emission)
        thermal.albedo = 0.1
        T = ((1 - 0.1) * S.value / (pi * 5.6703744191844314e-08)) ** (1 / 4) * u.K
        assert u.isclose(thermal.T0, T)

        thermal.eta = 1.5
        T = ((1 - 0.1) * S.value / (pi * 1.5 * 5.6703744191844314e-08)) ** (1 / 4) * u.K
        assert u.isclose(thermal.T0, T)

    def test_temperature(self):
        S = 1361.16646541 * u.W / u.m**2
        thermal = InstantaneousEquilibrium(S, 0, 1, 1, LambertianSurface())
        assert thermal.T(0 * u.deg) == thermal.T0

        T60 = (S.value / (2 * pi * 5.6703744191844314e-08)) ** (1 / 4) * u.K
        assert u.isclose(thermal.T(60 * u.deg), T60)
