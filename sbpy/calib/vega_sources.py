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


class VegaPhotometry:
    """Built-in Vega photometry.

    Format:
        {
            'filename': 'vega-photometry-willmer2018.json',
            'data': { ... },
            'description': 'Willmer (2018) flux densities',
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
    tab = ascii.read('apjsaabfdft4_ascii.txt')
    tab['(1)'].name = 'name'
    tab['(14)'].name = 'fluxd'
    tab['(3)'].name = 'lambda pivot'
    tab['(4)'].name = 'lambda eff'
    tab['fluxd'].unit = 'erg/(s cm2 AA)'
    tab['lambda eff'].unit = 'um'
    tab['lambda pivot'].unit = 'um'
    tab['name'][tab['name'] == 'WFC3_F098m'] = 'WFC3_F098M'

    phot = {}
    phot['data'] = {}
    for row in tab:
        name = row['name'].replace('_', ' ').replace('cfhtls', 'CFHTLS')
        phot['data'][name] = {
            'fluxd': [row['fluxd'], 'erg/(s cm2 AA)'],
            'lambda pivot': [row['lambda pivot'], 'um'],
            'lambda eff': [row['lambda eff'], 'um']
        }

    phot['meta'] = {
      'comments': ['Flux density of Vega)',
                   'Table 4 of Willmer (2018, ApJS 236, 47)']
    }
    with open('vega-photometry-willmer2018.json', 'w') as outf:
        json.dump(phot, outf)

    """
    Willmer2018 = {
        'filename': 'vega-photometry-willmer2018.json',
        'description': 'Willmer (2018) flux densities',
        'bibcode': '2018ApJS..236...47W'
    }
