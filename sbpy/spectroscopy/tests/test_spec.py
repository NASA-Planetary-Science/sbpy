import os

import astropy.units as u
from ..spectroscopy import nophotod_prodrate
from astropy.tests.helper import remote_data
from astropy.table import Table
import numpy as np


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


class MockResponseSpec(object):

    def __init__(self, filename):
        self.filename = data_path(filename)

    @property
    def text(self):
        with open(self.filename) as f:
            return f.read()
