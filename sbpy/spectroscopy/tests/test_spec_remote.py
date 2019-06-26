import os

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
from astropy.table import Table
from astroquery.lamda import Lamda
from astroquery.jplspec import JPLSpec
from ...activity.gas import Haser, photo_timescale
from ...data import Ephem
from .. import Spectrum


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)
