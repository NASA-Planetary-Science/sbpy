# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy bandpass Module

"""

__all__ = [
    'bandpass'
]

import os
from astropy.utils.data import get_pkg_data_filename


def bandpass(name):
    """Retrieve bandpass transmission spectrum from sbpy.


    Parameters
    ----------
    name : string
        Name of the bandpass, case insensitive.  See notes for
        available filters.


    Returns
    -------
    bp : `~synphot.SpectralElement`


    Notes
    -----
    Available filters:

    +-------------+---------------------------+
    | Name        | Source                    |
    +=============+===========================+
    | WFC3 F438W  | HST/WFC3 UVIS, v4         |
    +-------------+---------------------------+
    | WFC3 F606W  | HST/WFC3 UVIS, v4         |
    +-------------+---------------------------+
    | Johnson V   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | Cousins R   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | Cousins I   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | SDSS u      | SDSS, dated 2001          |
    +-------------+---------------------------+
    | SDSS g      | SDSS, dated 2001          |
    +-------------+---------------------------+
    | SDSS r      | SDSS, dated 2001          |
    +-------------+---------------------------+
    | SDSS i      | SDSS, dated 2001          |
    +-------------+---------------------------+
    | SDSS z      | SDSS, dated 2001          |
    +-------------+---------------------------+

    WFC3, Johnson V, Cousins R and I filters from [CDBS]_, SDSS filters
    from [SDSS]_.

    References
    ----------
    .. [CDBS] Space Telescope Science Institute.  Calibration Database
       System.  http://www.stsci.edu/hst/observatory/cdbs .

    .. [SDSS] Sloan Digital Sky Survey.  Camera.
       https://www.sdss.org/instruments/camera

    """

    try:
        import synphot
    except ImportError:
        raise ImportError('synphot is required.')

    name2file = {
        'wfc3 f438w': 'wfc3_uvis_f438w_004_syn.fits',
        'wfc3 f606w': 'wfc3_uvis_f606w_004_syn.fits',
        'johnson v': 'johnson_v_004_syn.fits',
        'cousins r': 'cousins_r_004_syn.fits',
        'cousins i': 'cousins_i_004_syn.fits',
        'sdss u': 'sdss-u.fits',
        'sdss g': 'sdss-g.fits',
        'sdss r': 'sdss-r.fits',
        'sdss i': 'sdss-i.fits',
        'sdss z': 'sdss-z.fits',
    }

    fn = get_pkg_data_filename(os.path.join(
        '..', 'photometry', 'data', name2file[name.lower()]))
    bp = synphot.SpectralElement.from_file(fn)
    return bp
