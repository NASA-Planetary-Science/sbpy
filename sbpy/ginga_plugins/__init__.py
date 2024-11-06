# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from warnings import warn

from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning

try:
    from ginga.misc.Bunch import Bunch
except ImportError:
    warn(AstropyWarning(
        'ginga is not present: sbpy.ginga_plugins will not run.'
    ))

    Bunch = None

warn(
    "sbpy.ginga_plugins is deprecated.  Use the `sbpy_ginga` module instead.",
    AstropyDeprecationWarning,
    stacklevel=2,
)

# path to these plugins
p_path = os.path.split(__file__)[0]


def setup_cometaryenhancements():
    spec = Bunch(path=os.path.join(p_path, 'CometaryEnhancements.py'),
                 module='CometaryEnhancements', klass='CometaryEnhancements',
                 workspace='dialogs')
    return spec
