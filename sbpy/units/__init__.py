# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Common planetary astronomy units.

"""

from .core import *


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from astropy.units.utils import generate_unit_summary as _generate_unit_summary
if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())
