# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy unit and quantity typing."""

from typing import Union
from packaging.version import Version

import astropy.units as u
from astropy import __version__ as astropy_version

type UnitLike = Union[str, u.Unit]
type SpectralQuantity = Union[
    u.Quantity[u.physical.length], u.Quantity[u.physical.frequency]
]
type SpectralFluxDensityQuantity = Union[
    u.Quantity[u.physical.spectral_flux_density],
    u.Quantity[u.physical.spectral_flux_density_wav],
]

if Version(astropy_version) < Version("6.1"):
    type SpectralRadianceQuantity = Union[
        u.Quantity[u.physical.spectral_flux_density / u.sr],
        u.Quantity[u.physical.spectral_flux_density_wav / u.sr],
    ]
else:
    type SpectralRadianceQuantity = Union[
        u.Quantity[u.physical.surface_brightness],
        u.Quantity[u.physical.surface_brightness_wav],
    ]
