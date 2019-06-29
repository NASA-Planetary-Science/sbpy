# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy Vega Module

Vega spectra for conversion to the Vega magnitude system.  Most
functionality requires the `synphot` package.

Users may use the built-in spectrum, or provide their own.


Classes
-------
Vega       - Class for Vega standards.

"""

from .core import *
from . import sources
