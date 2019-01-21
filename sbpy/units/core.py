# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy units core module

"""

__all__ = [
    'VegaFluxDensity',
    'VegaMag',
    'hundred_nm'
]

import astropy.units as u

VegaFluxDensity = u.def_unit(
    'VegaFluxDensity', doc='Represents the dust-free flux density of Vega')
VegaMag = u.MagUnit(VegaFluxDensity)
VegaMag.__doc__ = 'Vega-based magnitude'

# for spectral gradients
hundred_nm = u.def_unit('100 nm', represents=100 * u.nm)

u.add_enabled_units((VegaFluxDensity, VegaMag, hundred_nm))
