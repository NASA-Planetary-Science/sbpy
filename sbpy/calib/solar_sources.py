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
    'Castelli1996',
    'calspec',
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
        'description': ('E490-00a (2014) low resolution reference solar '
                        'spectrum (Table 4)'),
        'bibcode': 'doi:10.1520/E0490'
    }

    Kurucz1993 = {
        'filename': ('https://archive.stsci.edu/hlsps/reference-atlases/cdbs/'
                     'grid/k93models/standards/sun_kurucz93.fits'),
        'description': 'Kurucz (1993) model, scaled by Colina et al. (1996)',
        'bibcode': '1993KurCD..13.....K'
    }

    Castelli1996 = {
        'filename': ('https://archive.stsci.edu/hlsps/reference-atlases/cdbs/'
                     'grid/k93models/standards/sun_castelli.fits'),
        'description': ('Castelli model, scaled and presented by Colina et '
                        'al. (1996)'),
        'bibcode': '1996AJ....112..307C'
    }

    calspec = {
        "filename": ("https://archive.stsci.edu/hlsps/reference-atlases/cdbs/"
                     "current_calspec/sun_mod_001.fits"),
        "description": "R=5000, created by R. Bohlin from Kurucz Special Model",
        "bibcode": "2014PASP..126..711B"
    }


class SolarPhotometry:
    """Built-in solar photometry.

    Format:
        {
            'filename': 'solar-photometry-willmer2018.json',
            'data': { ... },
            'description': 'Willmer (2018) apparent magnitudes',
            'bibcode': '2018ApJS..236...47W'
        }

    Only one of 'filename' or 'data' required.  The file must be in
    the calib/data directory.

    Willmer data generated with source files downloaded from ApJS:

    https://iopscience.iop.org/0067-0049/236/2/47/suppdata/apjsaabfdft3_ascii.txt
    https://iopscience.iop.org/0067-0049/236/2/47/suppdata/apjsaabfdft4_ascii.txt

    Commented out header lines.

    from astropy.io import ascii
    from astropy.table import join
    import json
    sun = ascii.read('apjsaabfdft3_ascii.txt')
    filters = ascii.read('apjsaabfdft4_ascii.txt')
    tab = join(sun, filters, keys=('(1)',))
    tab['(1)'].name = 'name'
    tab['(7)_1'].name = 'fluxd'
    tab['(4)_2'].name = 'lambda eff'
    tab['(10)_1'].name = 'lambda pivot'
    tab['fluxd'].unit = 'mag(AB)'
    tab['lambda eff'].unit = 'um'
    tab['lambda pivot'].unit = 'um'
    tab['name'][tab['name'] == 'WFC3_F098m'] = 'WFC3_F098M'

    phot = {}
    phot['data'] = {}
    for row in tab:
        name = row['name'].replace('_', ' ').replace('cfhtls', 'CFHTLS')
        phot['data'][name] = {
            'fluxd': [10**(-0.4 * (21.1 + row['fluxd'])), 'erg/(s cm2 AA)'],
            'lambda pivot': [row['lambda pivot'], 'um'],
            'lambda eff': [row['lambda eff'], 'um']
        }

    phot['meta'] = {
      'comments': ['Apparent magnitude of the Sun',
                   'Tables 3 and 4 of Willmer (2018, ApJS 236, 47)']
    }
    with open('solar-photometry-willmer2018.json', 'w') as outf:
        json.dump(phot, outf)

    """
    Willmer2018 = {
        'filename': 'solar-photometry-willmer2018.json',
        'description': 'Willmer (2018) apparent magnitudes',
        'bibcode': '2018ApJS..236...47W'
    }
