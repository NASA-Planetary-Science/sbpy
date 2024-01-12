# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""sbpy photometric calibration

Most functionality requires the `synphot` package.

Users may use any of the built-in spectra or photometry for the Sun or
Vega, or provide their own.


Solar calibration
-----------------

The E490 spectra are included in the `sbpy` distribution.  The others
are downloaded and cached as needed::

  E490_2014 - E490-00a (2014) standard spectrum (ASTM 2014).
  E490_2014LR - Low resolution version of E490 (ASTM 2014).
  Kurucz1993 - Kurucz (1993) model scaled by Colina et al. (1996).
  Castelli1996 - Catelli model spectrum from Colina et al. (1996).
  calspec - STSCI CALSPEC R~5000 model solar spectrum, Bohlin et al. (2014).

Filter photometry from Willmer (2018).


Vega calibration
----------------

The spectrum of Bohlin (2014) is included with `sbpy`.

Filter photometry from Willmer (2018).


References
----------
ASTM 2014, E490-00a(2014) Standard Solar Constant and Zero Air Mass
Solar Spectral Irradiance Tables, ASTM International, West
Conshohocken, PA. doi:10.1520/E0490

Bohlin, R. C. 2014, Astronomical Journal, 147, 127.

Colina et al. 1996, Astronomical Journal 112, 307.

Kurucz 1993, Smithsonian Astropys. Obs. CD-ROM No. 13.

Willmer 2018, ApJS 236, 47.

"""

from .core import *
from . import solar_sources
from . import vega_sources
