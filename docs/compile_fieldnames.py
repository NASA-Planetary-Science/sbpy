"""This script will take sbpy.data.conf.fieldnames and turn it into a
   human-readable table
"""
import sys
from io import TextIOWrapper, BytesIO
from numpy import array
from sbpy.data import conf
from astropy.table import Table
from astropy.io import ascii

out = """
.. _propertynames:

List of Property Names
======================

The following table lists names for a wide range of properties that
are accepted and recognized by `sbpy`. Only these names should be used
in `~sbpy.data.DataClass` objects, i.e., `~sbpy.data.Ephem`,
`~sbpy.data.Orbit`, `~sbpy.data.Obs`, and `~sbpy.data.Phys` objects.

Names listed in the "Alternatives" column of the table can be used
synonymously. As an example, heliocentric distance can be addressed as
``'r'`` or ``'heldist'``:

    >>> from sbpy.data import Ephem
    >>> ceres = Ephem.from_horizons('Ceres')
    >>> print(ceres['r']) # doctest: +IGNORE_OUTPUT
    [2.69866993] AU
    >>> print(ceres['heldist']) # doctest: +IGNORE_OUTPUT
    [2.69866993] AU

The list of alternative field names is always up to date, but not
complete.

Developer Information
---------------------

The source list for this table is located in
``sbpy/sbpy/data/__init__.py`` as ``sbpy.data.conf.fieldnames``. If
you think an important alternative is missing, please suggest it by
opening an issue. However, keep in mind that each alternative field
name has to be *unique* and *unambiguous*.


List of Alternative Field Names
-------------------------------

"""

# build table
data = []
for parameter in conf.fieldnames:
    data.append(['**'+parameter[-1]+'**',
                 ', '.join(['``'+str(p)+'``' for p in parameter[:-1]])])
data = Table(array(data), names=('Description',
                                 'Alternative Names'))

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
with open('fieldnames.rst', 'w') as outf:
    outf.write(out)

sys.stdout.close()
