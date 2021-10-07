"""This script will take sbpy.data.Conf.fieldnames and turn it into a
   human-readable table
"""
import sys
from io import TextIOWrapper, BytesIO
from numpy import array
from sbpy.data import Conf
from astropy.table import Table
from astropy.io import ascii

out = """
.. _field name list:

sbpy Field Names
================

The following table lists field names that are recognized by `sbpy`
when accessing `~sbpy.data.DataClass` objects, i.e.,
`~sbpy.data.Ephem`, `~sbpy.data.Orbit`, or `~sbpy.data.Phys`
objects. Each row of the following table represents one property; for
each property it lists its description, acceptable field names,
provenance (which `~sbpy.data.DataClass` class should be used to store
this object so that `sbpy` uses it properly), and its physical
dimension (if any).

How do I use this Table?
------------------------

As an example, imagine you are interested in storing an object's right
ascension into a `~sbpy.data.DataClass` object. The field names table
tells you that you should name the field either ``ra`` or ``RA``, that
you should use either a `~sbpy.data.Ephem` or `~sbpy.data.Obs` object
to store the data in, and that the field data should be expressed as
angles. Based on this information, we can create a `~sbpy.data.Obs`
object (presuming that the data were derived from observations):

    >>> from sbpy.data import Obs
    >>> import astropy.units as u
    >>> obs = Obs.from_dict({'ra': [12.345, 12.346, 12.347]*u.deg})
    >>> obs['ra']  # doctest: +SKIP
    <Quantity [12.345, 12.346, 12.347] deg>

Since RA requires an angle as dimension, we use degrees, but we might
as well use radians - `sbpy` will convert the units where necessary.
RA has an alternative field name (``'RA'``), we can now use that name,
too, in order to retrieve the data:

    >>> obs['RA']  # doctest: +SKIP
    <Quantity [12.345, 12.346, 12.347] deg>


The field name list is always up to date, but it might not be
complete. If you think an important alternative name is missing,
please suggest it by opening an issue. However, keep in mind that each
alternative field name has to be **unique** and **unambiguous**. The
source list is located as ``sbpy.data.Conf.fieldnames`` in
``sbpy/data/__init__.py``.

Special Case: Epoch
-------------------

Please note that epochs generally have to be provided as
`~astropy.time.Time` objects. The advantage of using such objects is
that they can be readily transformed into a wide range of formats
(e.g., ISO, Julian Date, etc.) and time scales (e.g., UTC, TT, TDB,
etc.) Hence, `sbpy` requires that all fields referring to a point in
time be provided as `~astropy.time.Time` objects.

Field Name List
---------------

"""

# build table
data = []
for p in Conf.fieldnames_info:
    data.append(['**'+p['description']+'**',
                 ', '.join(['``'+str(f)+'``' for f in p['fieldnames']]),
                 ', '.join([{'orbit': '`~sbpy.data.Orbit`',
                             'ephem': '`~sbpy.data.Ephem`',
                             'obs': '`~sbpy.data.Obs`',
                             'phys': '`~sbpy.data.Phys`'}[
                                 m.replace(',', '')] for m in p['provenance']]),
                 str(p['dimension'])])
data = Table(array(data), names=('Description',
                                 'Field Names',
                                 'Provenance',
                                 'Dimension'))

# redirect stdout
sys.stdout = TextIOWrapper(BytesIO(), sys.stdout.encoding)

# ascii.write will write to stdout
datarst = ascii.write(data, format='rst')

# read formatted table data
sys.stdout.seek(0)
rsttab = sys.stdout.read()

# add formatted table data to out
out += rsttab

# write fieldnames.rst
with open('sbpy/data/fieldnames.rst', 'w') as outf:
    outf.write(out)

sys.stdout.close()
