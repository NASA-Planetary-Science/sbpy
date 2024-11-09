# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy unit and quantity typing."""

from typing import Union
import astropy.units as u

UnitLike = Union[str, u.Unit]
SpectralQuantity = Union[u.Quantity["length"], u.Quantity["frequency"]]
SpectralFluxDensityQuantity = Union[
    u.Quantity["spectral_flux_density"], u.Quantity["spectral_flux_density_wav"]
]
SpectralRadianceQuantity = Union[
    u.Quantity["surface_brightness"], u.Quantity["surface_brightness_wav"]
]
