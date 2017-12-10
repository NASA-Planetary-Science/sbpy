# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""================
SBPy Vega Module
================

Vega spectra for conversion to the Vega magnitude system.  Requires
the `synphot` package.

Users may use any of the built-in spectra, or provide their own.

Classes
-------
Vega       - Class for Vega standards.
Bohlin2014 - Vega spectrum from Bohlin (2014).


References
----------

Bohlin, R. C. 2014, Astronomical Journal, 147, 127.

"""

from .core import *
from . import sources
