# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from astropy.modeling.models import Scale
import astropy.constants as const
import synphot
from ..sources import BlackbodySource, SpectralSource, SynphotRequired, Reddening
from ..coma import Scattered, Thermal
from ..core import SpectralGradient
from ...activity.dust import Afrho, Efrho
from ...calib import Sun
from ...units import hundred_nm


class TestScattered:
    def test_init(self):
        s = Scattered({"rh": 1 * u.au, "delta": 1 * u.au}, "W/(m2 um)")
        assert s.rh == 1 * u.au
        assert s.delta == 1 * u.au
        assert s.unit == u.Unit("W/(m2 um)")

    def test_evaluate(self):
        cross_section = 1 * u.km**2
        S = 0 * u.percent / hundred_nm

        w = [2.0, 2.2] * u.um
        eph = {"rh": 1 * u.au, "delta": 1 * u.au}
        unit = u.Unit("W/(m2 um)")

        # expected value
        sun = Sun.from_default()
        expected = (
            cross_section
            * sun.observe(w, unit=unit)
            / np.pi
            / eph["rh"].to_value("au") ** 2
            / eph["delta"] ** 2
        ).to(unit)

        # test multiple wavelengths and no reddening
        s = Scattered(eph, unit, cross_section=cross_section, S=S)
        assert u.allclose(s(w), expected)

        # test single wavelength and add reddening
        w = 2 * u.um
        s.S = 1 * u.percent / hundred_nm
        R = 1 + ((w - 0.55 * u.um) * s.S).to_value("")
        expected = (
            cross_section
            * sun(w, unit=unit)
            * R
            / np.pi
            / eph["rh"].to_value("au") ** 2
            / eph["delta"] ** 2
        ).to(unit)
        assert u.allclose(s(w), expected)

    def test_afrho(self):
        cross_section = 1 * u.km**2
        S = 0 * u.percent / hundred_nm
        w = 2.0 * u.um
        eph = {"rh": 1 * u.au, "delta": 1 * u.au}
        rap = 1 * u.arcsec
        unit = u.Unit("W/(m2 um)")

        s = Scattered(eph, unit, cross_section=cross_section, S=S)
        a = s.afrho(w, rap)
        b = Afrho.from_cross_section(cross_section, 1, rap, eph)
        assert u.isclose(a, b)


class TestThermal:
    def test_init(self):
        s = Thermal({"rh": 1 * u.au, "delta": 1 * u.au}, "W/(m2 um)")
        assert s.rh == 1 * u.au
        assert s.delta == 1 * u.au
        assert s.unit == u.Unit("W/(m2 um)")

    def test_evaluate(self):
        cross_section = 1 * u.km**2
        Tscale = 1.0

        w = [20.0, 22.0] * u.um
        eph = {"rh": 1 * u.au, "delta": 1 * u.au}
        unit = u.Unit("W/(m2 um)")

        # expected value
        B = BlackbodySource(278 / np.sqrt(eph["rh"].to_value("au")))
        expected = (cross_section * B.observe(w, unit=unit) / eph["delta"] ** 2).to(
            unit
        )

        # test multiple wavelengths and no reddening
        t = Thermal(eph, unit, cross_section=cross_section, Tscale=Tscale)
        assert u.allclose(t(w), expected)

        # test single wavelength and hotter temperatures
        w = 2 * u.um
        Tscale = 1.1
        eph = {"rh": 0.8 * u.au, "delta": 1 * u.au}
        t = Thermal(eph, unit, cross_section=cross_section, Tscale=Tscale)
        B = BlackbodySource(Tscale * 278 / np.sqrt(eph["rh"].to_value("au")))
        expected = (cross_section * B(w, unit=unit) / eph["delta"] ** 2).to(unit)
        assert u.allclose(t(w), expected)

    def test_efrho(self):
        cross_section = 1 * u.km**2
        Tscale = 1.1
        w = 20.0 * u.um
        eph = {"rh": 1 * u.au, "delta": 1 * u.au}
        rap = 1 * u.arcsec
        unit = u.Unit("W/(m2 um)")

        t = Thermal(eph, unit, cross_section=cross_section, Tscale=Tscale)
        a = t.efrho(w, rap)
        b = Efrho.from_cross_section(cross_section, 1, rap, eph)
        assert u.isclose(a, b)
