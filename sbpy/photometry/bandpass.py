# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy bandpass Module

"""

__all__ = ["bandpass"]

import os
from astropy.utils.data import get_pkg_data_path
from ..utils.decorators import requires


@requires("synphot")
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
    See references for sources and additional notes.  Available filters:

    +-------------+---------------------------+
    | Name        | Source                    |
    +=============+===========================+
    | 2MASS J     | Cohen et al. 2003         |
    +-------------+---------------------------+
    | 2MASS H     | Cohen et al. 2003         |
    +-------------+---------------------------+
    | 2MASS Ks    | Cohen et al. 2003         |
    +-------------+---------------------------+
    | ATLAS c     | Tonry et al. 2018         |
    +-------------+---------------------------+
    | ATLAS o     | Tonry et al. 2018         |
    +-------------+---------------------------+
    | Cousins R   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | Cousins I   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | Johnson U   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | Johnson B   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | Johnson V   | STScI CDBS, v4            |
    +-------------+---------------------------+
    | LSST u      | LSST GitHub, v1.9         |
    +-------------+---------------------------+
    | LSST g      | LSST GitHub, v1.9         |
    +-------------+---------------------------+
    | LSST r      | LSST GitHub, v1.9         |
    +-------------+---------------------------+
    | LSST i      | LSST GitHub, v1.9         |
    +-------------+---------------------------+
    | LSST y      | LSST GitHub, v1.9         |
    +-------------+---------------------------+
    | LSST z      | LSST GitHub, v1.9         |
    +-------------+---------------------------+
    | PS1 g       | Tonry et al. 2012         |
    +-------------+---------------------------+
    | PS1 r       | Tonry et al. 2012         |
    +-------------+---------------------------+
    | PS1 i       | Tonry et al. 2012         |
    +-------------+---------------------------+
    | PS1 w       | Tonry et al. 2012         |
    +-------------+---------------------------+
    | PS1 y       | Tonry et al. 2012         |
    +-------------+---------------------------+
    | PS1 z       | Tonry et al. 2012         |
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
    | WFC3 F438W  | HST/WFC3 UVIS, v4         |
    +-------------+---------------------------+
    | WFC3 F606W  | HST/WFC3 UVIS, v4         |
    +-------------+---------------------------+
    | WISE W1     | Jarrett et al. 2011       |
    +-------------+---------------------------+
    | WISE W2     | Jarrett et al. 2011       |
    +-------------+---------------------------+
    | WISE W3     | Jarrett et al. 2011       |
    +-------------+---------------------------+
    | WISE W4     | Jarrett et al. 2011       |
    +-------------+---------------------------+


    References
    ----------

    .. [CDBS] Space Telescope Science Institute.  HST Calibration Reference
       Data System.  https://hst-crds.stsci.edu/ .

    .. [COH03] Cohen, M. et al. 2003.  Spectral Irradiance Calibration
       in the Infrared.  XIV.  The Absolute Calibration of 2MASS.  AJ
       126, 1090.

    .. [JAR11] Jarrett, T. H. et al. 2011.  The Spitzer-WISE Survey of
       the Ecliptic Poles. ApJ 735, 112.

    .. [LSST23] LSST Throughputs Repository.  Version 1.9 of the baseline
       "total" throughput curves are used, which include: mirrors, lenses,
       filter, detector, and the zenith atmosphere.
       https://github.com/lsst/throughputs  Commit 2273280.

    .. [SDSS] Sloan Digital Sky Survey.  Camera.
       www.sdss.org/instruments/camera .

    .. [TON12] Tonry, J. L. et al. 2012.  The Pan-STARRS1 Photometric
       System.  ApJ 750, 99.

    .. [TON18] Tonry, J. L. et al. 2018.  ATLAS: A High-cadence All-sky
       Survey System.  PASP 130, 064505.

    """

    import synphot

    name2file = {
        "2mass j": "2mass-j-rsr.txt",
        "2mass h": "2mass-h-rsr.txt",
        "2mass ks": "2mass-ks-rsr.txt",
        "atlas c": "atlas-c.txt",
        "atlas o": "atlas-o.txt",
        "cousins r": "cousins_r_004_syn.fits",
        "cousins i": "cousins_i_004_syn.fits",
        "johnson u": "johnson_u_004_syn.fits",
        "johnson b": "johnson_b_004_syn.fits",
        "johnson v": "johnson_v_004_syn.fits",
        "lsst u": "lsst_u.dat",
        "lsst g": "lsst_g.dat",
        "lsst r": "lsst_r.dat",
        "lsst i": "lsst_i.dat",
        "lsst z": "lsst_z.dat",
        "lsst y": "lsst_y.dat",
        "ps1 g": "ps1-gp1.txt",
        "ps1 r": "ps1-rp1.txt",
        "ps1 i": "ps1-ip1.txt",
        "ps1 w": "ps1-wp1.txt",
        "ps1 z": "ps1-zp1.txt",
        "ps1 y": "ps1-yp1.txt",
        "sdss u": "sdss-u.fits",
        "sdss g": "sdss-g.fits",
        "sdss r": "sdss-r.fits",
        "sdss i": "sdss-i.fits",
        "sdss z": "sdss-z.fits",
        "wfc3 f438w": "wfc3_uvis_f438w_004_syn.fits",
        "wfc3 f606w": "wfc3_uvis_f606w_004_syn.fits",
        "wise w1": "WISE-RSR-W1.EE.txt",
        "wise w2": "WISE-RSR-W2.EE.txt",
        "wise w3": "WISE-RSR-W3.EE.txt",
        "wise w4": "WISE-RSR-W4.EE.txt",
    }
    wave_unit = {
        "lsst u": "nm",
        "lsst g": "nm",
        "lsst r": "nm",
        "lsst i": "nm",
        "lsst z": "nm",
        "lsst y": "nm",
        "ps1 g": "nm",
        "ps1 r": "nm",
        "ps1 i": "nm",
        "ps1 w": "nm",
        "ps1 z": "nm",
        "ps1 y": "nm",
    }

    fn = get_pkg_data_path(
        os.path.join("..", "photometry", "data", name2file[name.lower()])
    )

    # wave_unit is deprecated for FITS files in synphot 1.4
    kwargs = {}
    if name.lower() in wave_unit:
        kwargs["wave_unit"] = wave_unit[name.lower()]

    bp = synphot.SpectralElement.from_file(fn, **kwargs)
    return bp
