# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
===============
SBPy Sun Module
===============

Solar spectra and observations thereof.  Requires the `synphot` package.

Users may use any of the built-in spectra, or provide their own.  The
package allows any spectrum to be input (solar, stellar, or
otherwise).

A few spectra are included in the `sbpy` distribution, others are
downloaded and cached as needed::

  E490_2014 - E490-00a (2014) standard spectrum (ASTM 2014).
  E490_2014LR - Low resolution version of E490 (ASTM 2014).
  Kurucz1993 - Kurucz (1993) model scaled by Colina et al. (1996).
  Castelli1996 - Catelli model spectrum from Colina et al. (1996).


Classes
-------
Sun - Solar spectrum.

Context Managers
----------------
default_sun - Get/set the default solar spectrum for sbpy.

References
----------
ASTM 2014, E490-00a(2014) Standard Solar Constant and Zero Air Mass
Solar Spectral Irradiance Tables, ASTM International, West
Conshohocken, PA. doi:10.1520/E0490

Kurucz 1993, Smithsonian Astropys. Obs. CD-ROM No. 13.
 
Colina et al. 1996, Astronomical Journal 112, 307.

"""

from .core import *
from . import sources
