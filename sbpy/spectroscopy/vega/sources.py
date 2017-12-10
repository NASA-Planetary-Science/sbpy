# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=======================
SBPy Vega Sources Module
=======================

Descriptions of source Vega spectra.

"""

# Parameters passed to Vega.from_file, and 'file' is a URL.  After adding spectra here, update __init__.py docstring.

available = [
    'Bohlin2014'
]

Bohlin2014 = {
    'file': 'ftp://ftp.stsci.edu/cdbs/calspec/alpha_lyr_stis_008.fits',
    'description': 'Vega spectrum.',
    'bibcode': '2014AJ....147..127B'
}
