# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""===============
SBPy Sun Module
===============

Solar spectra and observations thereof.  Requires the `synphot` package.

Users may use any of the built-in spectra, or provide their own.  The package allows any spectrum to be input (solar, stellar, or otherwise).

A few spectra are included in the `sbpy` distribution, others are downloaded and cached as needed.


Classes
-------
Sun - Solar spectrum.


Context Managers
----------------
default_sun - Get/set the default solar spectrum for sbpy.

"""

try:
    import synphot
except ImportError as e:
    from warnings import warn
    from astropy.utils.exceptions import AstropyWarning
    warn(AstropyWarning('synphot is required for the sun module.'))

from .core import *
from . import sources
