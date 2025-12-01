# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Common planetary astronomy units."""

import astropy
from packaging.version import Version as Version

from .core import *


###########################################################################
# DOCSTRING
# This generates a docstring for this module that describes all of the
# standard units defined here.
if Version(astropy.__version__) < Version("7.2.0"):
    from astropy.units.utils import generate_unit_summary as _generate_unit_summary
else:
    from astropy.units.docgen import generate_unit_summary as _generate_unit_summary

if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())
