# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""===============
SBPy Sun Module
===============

Solar spectra and observations thereof.  Requires the `synphot` package.

Users may use any of the built-in spectra, or provide their own.  The package allows any spectrum to be input (solar, stellar, or otherwise).

Classes
-------

Sun          - Solar spectrum.

Castelli1996 - Castelli high resolution (1-Å resolution) solar
               spectrum (Colina et al. 1996), 2001.5 Å to 2.5 μm.

E490_2014    - E490-00a composite solar spectrum (Table 3 of ASTM 2014),
               1195 Å to 1.0 mm.  Spectral resolution ~100 in the UV,
               ~1000 in the optical to near-IR, ~100 in the near- to
               mid-IR, ~1 beyond 10 μm.

E490_2014LR  - E490-00a low-resolution composite solar spectrum (Table
               4 of ASTM 2014), 1195 Å to 1.0 mm.  Spectral resolution
               ~1 to ~10.

Kurucz93     - Kurucz (1993) low resolution solar spectrum, scaled by
               Colina et al. (1996), 90.9 Å to 160 μm.


References
----------

ASTM 2014, E490-00a(2014) Standard Solar Constant and Zero Air Mass
Solar Spectral Irradiance Tables, ASTM International, West
Conshohocken, PA. doi:10.1520/E0490

Kurucz 1993, Smithsonian Astropys. Obs. CD-ROM No. 13.

Colina et al. 1996, Astronomical Journal 112, 307.

"""

try:
    from .core import *
except ImportError as e:
    from warnings import warn
    from astropy.utils.exceptions import AstropyWarning
    warn(AstropyWarning('synphot is required for the sun module.'))
