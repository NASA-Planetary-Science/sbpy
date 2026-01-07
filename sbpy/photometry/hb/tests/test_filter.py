# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import astropy.units as u

from ..filter import HBFilterSet, gamma, gamma_prime


class TestHBFilterSet:
    def test_str(self):
        assert str(HBFilterSet.UC) == "UC"
        assert str(HBFilterSet.COplus) == "CO+"
        assert str(HBFilterSet.H2Oplus) == "H2O+"

    def test_comparisons(self):
        assert HBFilterSet.CN == HBFilterSet.CN
        assert HBFilterSet.C3 < HBFilterSet.BC
        assert HBFilterSet.GC > HBFilterSet.C2

    def test_designation(self):
        assert HBFilterSet.OH.designation == "3090/62"

    def test_wavelength(self):
        assert HBFilterSet.NH.wavelength == 3361 * u.AA

    def test_widths(self):
        widths = HBFilterSet.RC.widths
        assert widths[80] == 53 * u.AA
        assert widths[50] == 58 * u.AA
        assert widths[10] == 71 * u.AA
        assert widths[1] == 92 * u.AA

    def test_fluxd0(self):
        assert HBFilterSet.H2Oplus.fluxd0 == 1.380e-9 * u.erg / u.cm**2 / u.s / u.AA

    def test_solar_color(self):
        assert HBFilterSet.C3.solar_color == 0.497 * u.mag

    def test_gamma(self):
        assert HBFilterSet.C2.gamma("C2") == 5.433e-3 / u.AA
        assert HBFilterSet.NH.gamma("NH") == 1.907e-2 / u.AA
        assert HBFilterSet.NH.gamma("C3") == 1.433e-5 / u.AA
        assert HBFilterSet.NH.gamma("CN") == 0 / u.AA

        with pytest.raises(ValueError):
            HBFilterSet.BC.gamma("C2")

    def test_gamma_prime(self):
        assert HBFilterSet.OH.gamma_prime("OH") == 0.98
        assert HBFilterSet.COplus.gamma_prime("CO+") == 0.99
        assert HBFilterSet.COplus.gamma_prime("C3") == 4.607e-4 * 0.99 / 1.549e-2
        assert HBFilterSet.COplus.gamma_prime("C2") == 0

        with pytest.raises(ValueError):
            HBFilterSet.GC.gamma_prime("C3")


def test_gamma():
    assert gamma("C3", "C3") == 3.352e-3 / u.AA
    assert gamma("CN", "CN") == 1.812e-2 / u.AA
    assert gamma("CN", "C3") == 1.427e-3 / u.AA
    assert gamma("CN", "NH") == 0 / u.AA


def test_gamma_prime():
    assert gamma_prime("H2O+", "H2O+") == 1.00
    assert gamma_prime("C2", "C2") == 0.66
    assert gamma_prime("C2", "C3") == 0
    assert gamma_prime("C2", "OH") == 0
