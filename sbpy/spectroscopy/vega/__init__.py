# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""================
SBPy Vega Module
================

Vega spectra for conversion to the Vega magnitude system.  Requires
the `synphot` package.

Users may use any of the built-in spectra, or provide their own.

Classes
-------
Bohlin2014 - Vega spectrum from Bohlin (2014).
UserSun    - User-provided Vega spectrum.


References
----------

Bohlin, R. C. 2014, Astronomical Journal, 147, 127.

"""

try:
    from .core import *
except ImportError as e:
    from warnings import warn
    from astropy.utils.exceptions import AstropyWarning
    warn(AstropyWarning('synphot is required for the sun module.'))
