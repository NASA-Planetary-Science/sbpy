# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
================
SBPy Vega Module
================

Vega spectra for conversion to the Vega magnitude system.  Requires
the `synphot` package.

Users may use the built-in spectrum, or provide their own.


Classes
-------
Vega       - Class for Vega standards.

Context Managers
----------------
default_vega - Get/set the default Vega spectrum for sbpy.

References
----------
Bohlin, R. C. 2014, Astronomical Journal, 147, 127.

"""

from .core import *
from . import sources
