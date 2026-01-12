# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
HB cometary narrow-band filter set.

Filter set characteristics and data reduction tools.

"""

__all__ = ["filter_data", "HBFilterSet", "gamma", "gamma_prime"]

from enum import Enum
from typing import TypeVar

import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable
from astropy.utils.data import get_pkg_data_path

HBFilterSetType = TypeVar("HBFilterSetType", bound="HBFilterSet")

"""HB filter set characteristics and photometric properties."""
filter_data = {
    row["name"]: row
    for row in QTable(ascii.read(get_pkg_data_path("data", "filter.ecsv")))
}


class HBFilterSet(Enum):
    """Representative filter characteristics and photometric properties.


    Examples
    --------

    Get the wavelength of the BC filter:

    >>> from sbpy.photometry.hb import HBFilterSet
    >>> HBFilterSet.BC.wavelength
    <Quantity 4453. Angstrom>


    References
    ----------

    Farnham, Schleicher, & A'Hearn 2000. The HB Narrowband Comet Filters:
    Standard Stars and Calibrations.  Icarus 147, 180–204.

    """

    OH = "OH"
    NH = "NH"
    UC = "UC"
    CN = "CN"
    C3 = "C3"
    COplus = "CO+"
    BC = "BC"
    C2 = "C2"
    GC = "GC"
    H2Oplus = "H2O+"
    RC = "RC"

    def __str__(self) -> str:
        """Filter ID as a string."""
        return self.value

    def __lt__(self, filter: HBFilterSetType) -> bool:
        """Compare wavelengths."""
        return self.wavelength < filter.wavelength

    def __gt__(self, filter: HBFilterSetType) -> bool:
        """Compare wavelengths."""
        return self.wavelength > filter.wavelength

    @property
    def designation(self) -> str:
        """Filter designation."""

        return filter_data[self.value]["designation"]

    @property
    def wavelength(self) -> u.Quantity[u.AA]:
        """Center wavelength."""

        return filter_data[self.value]["cwave"]

    @property
    def widths(self) -> dict[int, u.Quantity[u.AA]]:
        """Power point width.

        Returns
        -------
        widths : dict
            Measured full-width power points at 80%, 50%, 10%, and 1% of the
            peak transmission.

        """

        power_percent = [80, 50, 10, 1]

        return {power_percent[i]: filter_data[self.value]["width"][i] for i in range(4)}

    @property
    def fluxd0(self) -> u.Quantity[u.Unit("erg / (cm2 s AA)")]:
        """Spectral flux density of a 0 magnitude star."""
        return filter_data[self.value]["Flam0"]

    @property
    def solar_color(self) -> u.Magnitude:
        """Solar color: [filter]-BC."""
        return u.Magnitude(filter_data[self.value]["XXmBCsun"])

    def gamma(self, molecule: str) -> u.Quantity[1 / u.AA]:
        """Fraction of molecule present in this filter, normalized by filter equivalent width.

        Tables VI and VII of Farnham, Schleicher, & A'Hearn 2000.


        Parameters
        ----------
        molecule : str
            The molecule name: OH, NH, CN, C3, CO+, C2, or H2O+.


        Returns
        -------
        gamma : |Quantity|

        """

        if self.value in ["UC", "BC", "GC", "RC"]:
            raise ValueError("gamma is not defined for the continuum filters")
        elif self.value == molecule:
            # Table VI of Farnham, Schleicher, & A'Hearn 2000
            g = filter_data[self.value]["gamma_XX_XX"]
        elif molecule == "C3":
            # Table VII of Farnham, Schleicher, & A'Hearn 2000
            table_vii = {
                "NH": 1.433e-5 / u.AA,
                "CN": 1.427e-3 / u.AA,
                "CO+": 4.607e-4 / u.AA,
            }
            g = table_vii.get(self.value, 0 / u.AA)
        else:
            g = 0 / u.AA

        return g.to("1/AA")

    def gamma_prime(self, molecule: str) -> u.Quantity[u.dimensionless_unscaled]:
        """Fraction of molecule present in this filter.

        Tables VI and VII of Farnham, Schleicher, & A'Hearn 2000.


        Parameters
        ----------
        molecule : str
            The molecule name: OH, NH, CN, C3, CO+, C2, or H2O+.


        Returns
        -------
        gamma_prime : |Quantity|

        """

        if self.value in ["UC", "BC", "GC", "RC"]:
            raise ValueError("gamma is not defined for the continuum filters")
        elif self.value == molecule:
            # Table VI of Farnham, Schleicher, & A'Hearn 2000
            gp = filter_data[self.value]["gamma_prime_XX_XX"]
        elif molecule == "C3":
            data = filter_data[self.value]
            equivalent_width = data["gamma_prime_XX_XX"] / data["gamma_XX_XX"]
            gp = self.gamma(molecule) * equivalent_width
        else:
            gp = 0

        return u.Quantity(gp, u.dimensionless_unscaled)


def gamma(filter: str, molecule: str):
    """Fraction of molecule present filter, normalized by filter equivalent width.

    Tables VI and VII of Farnham, Schleicher, & A'Hearn 2000.


    Parameters
    ----------
    filter : str
        The filter name: OH, NH, CN, C3, CO+, C2, or H2O+.

    molecule : str
        The molecule name: OH, NH, CN, C3, CO+, C2, or H2O+.


    Returns
    -------
    gamma : |Quantity|


    Examples
    --------

    >>> from sbpy.photometry import hb
    >>> hb.gamma("CN", "CN")
    <Quantity 0.01812 1 / Angstrom>

    Using this function is equivalent to using :func:`HBFilterSet.gamma`:

    >>> from sbpy.photometry.hb import HBFilterSet
    >>> HBFilterSet.CN.gamma("CN")
    <Quantity 0.01812 1 / Angstrom>

    """

    return HBFilterSet(filter).gamma(molecule)


def gamma_prime(filter: str, molecule: str) -> u.Quantity[u.dimensionless_unscaled]:
    """Fraction of molecule present in filter.

    Tables VI and VII of Farnham, Schleicher, & A'Hearn 2000.


    Parameters
    ----------
    filter : str
        The filter name: OH, NH, CN, C3, CO+, C2, or H2O+.

    molecule : str
        The molecule name: OH, NH, CN, C3, CO+, C2, or H2O+.


    Returns
    -------
    gamma_prime : |Quantity|


    Examples
    --------

    >>> from sbpy.photometry import hb
    >>> hb.gamma_prime("OH", "OH")
    <Quantity 0.98>

    Using this function is equivalent to using :func:`HBFilterSet.gamma_prime`:

    >>> from sbpy.photometry.hb import HBFilterSet
    >>> HBFilterSet.OH.gamma_prime("OH")
    <Quantity 0.98>

    """

    return HBFilterSet(filter).gamma_prime(molecule)
