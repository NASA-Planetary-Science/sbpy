# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ginga.misc.Bunch import Bunch

# path to these plugins
import os.path
p_path = os.path.split(__file__)[0]


def setup_cometaryenhancements():
    spec = Bunch(path=os.path.join(p_path, 'CometaryEnhancements.py'),
                 module='CometaryEnhancements', klass='CometaryEnhancements',
                 workspace='dialogs')
    return spec
