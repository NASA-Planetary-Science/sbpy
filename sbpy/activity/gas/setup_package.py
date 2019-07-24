# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os


def get_package_data():
    return {'sbpy.activity.gas': [os.path.join('data', 'schleicher88.txt')]}
