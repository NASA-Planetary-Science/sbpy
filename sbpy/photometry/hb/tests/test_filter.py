# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u

from ..filter import HBFilterSet


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
