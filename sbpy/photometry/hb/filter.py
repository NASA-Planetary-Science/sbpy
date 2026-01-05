# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
HB cometary narrow-band filter set.

Filter set characteristics and data reduction tools.

"""

from enum import Enum
from typing import TypeVar

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable
from astropy.utils.data import get_pkg_data_path


HBFilterSetType = TypeVar("HBFilterSetType", bound="HBFilterSet")

"""HB filter set characteristics and photometric properties."""
filter_data = {
    row["sbpy name"]: row
    for row in QTable(ascii.read(get_pkg_data_path("data", "filter.ecsv")))
}


class HBFilterSet(Enum):
    """Representative filter characteristics and photometric properties.


    Examples
    --------

    Get the wavelength of the BC filter:

    >>> from sbpy.photometry.hb import HBFilterSet
    >>> HBFilterSet.BC.wavelength
    <Quantity 445.3 nm>


    References
    ----------

    Farnham, Schleicher, & A'Hearn 2000. The HB Narrowband Comet Filters:
    Standard Stars and Calibrations.  Icarus 147, 180–204.

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
        return filter_data[self.name]["name"]

    def __lt__(self, filter: HBFilterSetType) -> bool:
        """Compare wavelengths."""
        return self.wavelength < filter.wavelength

    def __gt__(self, filter: HBFilterSetType) -> bool:
        """Compare wavelengths."""
        return self.wavelength > filter.wavelength

    @property
    def designation(self) -> str:
        """Filter designation."""

        return filter_data[self.name]["designation"]

    @property
    def wavelength(self) -> u.Quantity["nm"]:
        """Center wavelength."""

        return filter_data[self.name]["cwave"]

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

        return {power_percent[i]: filter_data[self.name]["width"][i] for i in range(4)}

    @property
    def fluxd0(self) -> u.Quantity["W / (m2 um)"]:
        """Spectral flux density of a 0 magnitude star."""
        return filter_data[self.name]["Flam0"]

    @property
    def XXmBCsun(self) -> u.mag:
        """Solar color: XX-BC."""
        return filter_data[self.name]["XXmBCsun"]
