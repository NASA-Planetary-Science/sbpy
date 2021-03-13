# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
sbpy.data.dimensions
--------------------

DataClass field dimensions, units, object types.

"""

import astropy.units as u
from astropy.time import Time
from astropy.table import Column


class FieldDimension:
    """DataClass field dimensions.


    Parameters
    ----------
    dimension : string
        Description of the physical dimension, e.g., 'angle', 'inverse time',
        'time-area', 'length per time'.  See Notes for format rules.  This
        string should match a value in ``sbpy.data.Conf.fieldnames_info``.

    unit : `~astropy.unit.Unit` or tuple
        Example units of the dimension.

    obj : Object or tuple, optional
        Expected object type, default is `~astropy.unit.Quantity`.

    equivalences : list, optional
        List of astropy unit equivalencies, e.g., `~astro.unit.spectral()`.


    Notes
    -----
    Dimension format rules:
        * Use 'inverse dimension' rather than '1 / dimension'.
        * Use 'dimension per dimension' rather than 'dimension / dimension'.
        * Use 'dimesion-dimension' rather than 'dimension*dimension'.

    After adding a new dimension, edit test_verify_fields and
    test_verify_fields_error in test/test_dataclass.py, and
    DataClass.verify_fields().

    Do not add dimensions that are not quantities (e.g., target name).

    """

    def __init__(self, dimension, unit, obj=u.Quantity, equivalencies=None):
        self.dimension = dimension
        self.unit = unit
        self.object = obj
        self.equivalencies = equivalencies

    def __str__(self):
        """Format as string."""
        return self.dimension


dimensionless = FieldDimension(
    '', u.dimensionless_unscaled, obj=(u.Quantity, Column))
angle = FieldDimension('angle', u.radian)
angle_per_time = FieldDimension('angle per time', u.radian / u.second)
energy = FieldDimension('energy', u.Joule)
frequency = FieldDimension('frequency', u.Hertz, equivalencies=u.spectral())
inverse_area = FieldDimension('inverse area', u.meter**-2)
inverse_time = FieldDimension('inverse time', u.second**-1)
length = FieldDimension('length', u.meter)
magnitude = FieldDimension('magnitude', u.mag)
magnitude_per_solid_angle = FieldDimension(
    'magnitude per solid angle', u.mag / u.steradian)
percent = FieldDimension('percent', u.percent)
solid_angle = FieldDimension('solid angle', u.steradian)
temperature = FieldDimension(
    'temperature', u.Kelvin, equivalencies=u.temperature())
time = FieldDimension('time', u.second)
time_area = FieldDimension('time-area', u.second * u.meter**2)
length_per_time = FieldDimension('length per time', u.meter / u.second)
time_object = FieldDimension('`~astropy.time.Time`', None, obj=Time)
