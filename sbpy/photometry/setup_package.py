# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os


def get_package_data():
    return {'sbpy.photometry': [os.path.join('data', '*.fits'),
                                os.path.join('data', '*.txt')]}
