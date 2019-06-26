# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os


def get_package_data():
    return {'sbpy.units.tests': [os.path.join('data',
                                 'hi05070405_9000036-avg-spec.txt')]}
