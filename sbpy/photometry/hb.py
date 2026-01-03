# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
HB cometary narrow-band filter set reduction tools.
"""

from enum import Enum
from typing import TypeVar
import astropy.units as u


HBFilterSetType = TypeVar("HBFilterSetType", bound="HBFilterSet")


class HBFilterSet(Enum):
    """Representative filter characteristics (Farnham & Schleicher 2000).


    Examples
    --------

    Get the wavelength of the BC filter:

    >>> from sbpy.photometry.hb import HBFilterSet
    >>> HBFilterSet.BC.wavelength
    <Quantity 445.3 nm>
    """

    OH = 0
    NH = 1
    UC = 2
    CN = 3
    C3 = 4
    COplus = 5
    BC = 6
    C2 = 7
    GC = 8
    H2Oplus = 9
    RC = 10

    def __str__(self) -> str:
        """Filter ID as a string."""
        ions = {"COplus": "CO+", "H2Oplus": "H2O+"}
        return ions.get(self.name, self.name)

    def __lt__(self, filter: HBFilterSetType) -> bool:
        """Compare wavelengths."""
        return self.wavelength < filter.wavelength

    def __gt__(self, filter: HBFilterSetType) -> bool:
        """Compare wavelengths."""
        return self.wavelength > filter.wavelength

    @property
    def designation(self) -> str:
        """Filter designation."""

        designations = [
            "3090/62",
            "3362/58",
            "3448/84",
            "3870/62",
            "4062/62",
            "4266/64",
            "4450/67",
            "5141/118",
            "5260/56",
            "7020/170",
            "7128/58",
        ]

        return designations[self.value]

    @property
    def wavelength(self) -> u.Quantity["nm"]:
        """Center wavelength."""

        wavelengths = [
            309.7,
            336.1,
            344.9,
            386.9,
            406.3,
            426.6,
            445.3,
            513.5,
            525.9,
            702.8,
            713.3,
        ]
        return u.Quantity(wavelengths[self.value], u.nm)

    @property
    def widths(self) -> dict[int, u.Quantity["nm"]]:
        """Power point width.

        Returns
        -------
        widths : dict
            Measured full-width power points at 80%, 50%, 10%, and 1% of the
            peak transmission.

        """

        power_percent = [80, 50, 10, 1]

        widths = [
            [5.2, 5.8, 6.8, 8.7],
            [4.7, 5.4, 6.4, 8.1],
            [7.2, 7.9, 9.3, 11.6],
            [5.0, 5.6, 6.5, 8.2],
            [4.3, 5.8, 6.9, 8.4],
            [5.8, 6.4, 7.4, 9.0],
            [5.5, 6.1, 7.1, 8.6],
            [10.9, 11.9, 14.0, 17.1],
            [5.2, 5.6, 6.5, 7.9],
            [14.8, 16.4, 19.3, 23.9],
            [5.3, 5.8, 7.1, 9.2],
        ]

        return {
            power_percent[i]: u.Quantity(widths[self.value][i], u.nm) for i in range(4)
        }
