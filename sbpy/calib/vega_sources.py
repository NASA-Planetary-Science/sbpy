# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy Vega Sources

Built-in Vega data.

Spectra parameters are passed to `Vega.from_file()`.  'filename' is a
URL or name of a file distributed with sbpy and located in
sbpy/calib/data (see also sbpy/calib/setup_package.py).  After adding
spectra or photometry here, update ``__init__.py`` docstring,
and docs/calib.rst.

"""


class VegaSpectra:
    Bohlin2014 = {
        'filename': 'alpha_lyr_stis_008-edit.fits',
        'description': 'Dust-free template spectrum of Bohlin 2014',
        'bibcode': '2014AJ....147..127B'
    }


"""Willmer table generated with source files downloaded from ApJS:

https://iopscience.iop.org/0067-0049/236/2/47/suppdata/apjsaabfdft3_ascii.txt
https://iopscience.iop.org/0067-0049/236/2/47/suppdata/apjsaabfdft4_ascii.txt

Commented out header lines.

from astropy.io import ascii
from astropy.table import join
phot = ascii.read('apjsaabfdft4_ascii.txt')
phot['(1)'].name = 'name'
phot['(14)'].name = 'fluxd'
phot['(3)'].name = 'lambda pivot'
phot['(4)'].name = 'lambda eff'
phot['fluxd'].unit = 'erg/(s cm2 AA)'
phot['lambda eff'].unit = 'um'
phot['lambda pivot'].unit = 'um'
phot['name'][phot['name'] == 'WFC3_F098m'] = 'WFC3_F098M'
phot.keep_columns(('name', 'fluxd', 'lambda pivot', 'lambda eff'))
for i in range(len(phot)):
    phot[i]['name'] = phot[i]['name'].replace('_', ' ').replace('cfhtls', 'CFHTLS')
phot.meta['comments'] = [
  'Flux density of Vega',
  'Table 4 of Willmer (2018, ApJS 236, 47)'
]
phot.write('vega-photometry-willmer2018.csv', format='ascii.ecsv', delimiter=',')
"""


class VegaPhotometry:
    Willmer2018 = {
        'filename': 'vega-photometry-willmer2018.csv',
        'description': 'Willmer (2018) flux densities',
        'bibcode': '2018ApJS..236...47W'
    }
