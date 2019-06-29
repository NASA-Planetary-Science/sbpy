# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
========================
SBPy Vega Sources Module
========================

Descriptions of source Vega spectra.

"""

# Parameters passed to Vega.from_file. 'filename' is a URL.  After
# adding spectra here update __init__.py docstring and
# docs/sbpy/spectroscopy.rst.

available = [
    'Bohlin2014'
]

Bohlin2014 = {
    'filename': 'alpha_lyr_stis_008-edit.fits',
    'description': 'Dust-free template spectrum of Bohlin 2014',
    'bibcode': '2014AJ....147..127B'
}
