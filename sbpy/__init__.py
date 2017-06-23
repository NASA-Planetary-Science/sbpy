# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    from .bib import *
    from .data import *
    from .activity import *
    from .photometry import *
    from .imageanalysis import *
    from .spectroscopy import *
    from .obsutil import *
