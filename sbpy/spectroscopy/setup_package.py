# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import

import os


def get_package_data():

    paths_test = [os.path.join('data', 'HCN.data'),
                  os.path.join('data', 'HCN.csv')]

    return {'astroquery.jplspec.tests': paths_test,}
