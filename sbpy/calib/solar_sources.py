# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy Solar Spectra

Built-in solar spectra for use with `~sbpy.calib.solar_spectrum`.


These parameters are passed to `Sun.from_file()`.  'filename' is a URL
or name of a file distributed with sbpy and located in sbpy/calib/data
(see also sbpy/calib/setup_package.py).  After adding spectra here,
update ``__init__.py`` docstring, and docs/calib.rst.

"""

sources = [
    'E490_2014',
    'E490_2014LR',
    'Kurucz1993',
    'Castelli1996'
]


class SolarSpectra:
    E490_2014 = {
        'filename': 'e490-00a_2014_hires.csv',
        'wave_unit': 'um',
        'flux_unit': 'W/(m2 um)',
        'description': 'E490-00a (2014) reference solar spectrum (Table 3)',
        'bibcode': 'doi:10.1520/E0490'
    }

    E490_2014LR = {
        'filename': 'e490-00a_2014_hires.csv',
        'wave_unit': 'um',
        'flux_unit': 'W/(m2 um)',
        'description': 'E490-00a (2014) low resolution reference solar spectrum (Table 4)',
        'bibcode': 'doi:10.1520/E0490'
    }

    Kurucz1993 = {
        'filename': 'ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/sun_kurucz93.fits',
        'description': 'Kurucz (1993) model, scaled by Colina et al. (1996)',
        'bibcode': '1993KurCD..13.....K'
    }

    Castelli1996 = {
        'filename': 'ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/sun_castelli.fits',
        'description': 'Castelli model, scaled and presented by Colina et al. (1996)',
        'bibcode': '1996AJ....112..307C'
    }


"""Willmer table generated with source files downloaded from ApJS:

https://iopscience.iop.org/0067-0049/236/2/47/suppdata/apjsaabfdft3_ascii.txt
https://iopscience.iop.org/0067-0049/236/2/47/suppdata/apjsaabfdft4_ascii.txt

Commented out header lines.

from astropy.io import ascii
from astropy.table import join
sun = ascii.read('apjsaabfdft3_ascii.txt')
filters = ascii.read('apjsaabfdft4_ascii.txt')
tab = join(sun, filters, keys=('(1)',))
solar_phot = tab[['(1)', '(6)_1', '(4)_2', '(10)_1']]
solar_phot['(1)'].name = 'name'
solar_phot['(6)_1'].name = 'fluxd'
solar_phot['(4)_2'].name = 'lambda eff'
solar_phot['(10)_1'].name = 'lambda pivot'
solar_phot['fluxd'].unit = 'mag(AB)'
solar_phot['lambda eff'].unit = 'um'
solar_phot['lambda pivot'].unit = 'um'
solar_phot['name'][solar_phot['name'] == 'WFC3_F098m'] = 'WFC3_F098M'
for i in range(len(solar_phot)):
    solar_phot[i]['name'] = solar_phot[i]['name'].replace('_', ' ').replace('cfhtls', 'CFHTLS')
solar_phot.meta['comments'] = [
  'Apparent magnitude of the Sun, (AB mag)',
  'Tables 3 and 4 of Willmer (2018, ApJS 236, 47)'
]
solar_phot.write('solar-photometry-willmer2018.csv', format='ascii.ecsv', delimiter=',')
"""


class SolarPhotometry:
    Willmer2018 = {
        'filename': 'solar-photometry-willmer2018.csv',
        'description': 'Willmer (2018) apparent mangitudes',
        'bibcode': '2018ApJS..236...47W'
    }
