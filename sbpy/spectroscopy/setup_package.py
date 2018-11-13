# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import

import os


def get_package_data():

    paths_test = [os.path.join('data', 'HCN.csv'),
                  os.path.join('data', 'CH3OH.csv'),
                  os.path.join('data', 'CO.csv')]

    return {'sbpy.spectroscopy.tests': paths_test, }
