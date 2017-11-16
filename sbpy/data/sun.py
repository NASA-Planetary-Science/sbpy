# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
===============
SBPy Sun Module
===============

This module requires `synphot` and solar spectra available in Space
Telescope Science Institute's `Calibration Reference Data System
<http://www.stsci.edu/hst/observatory/crds/throughput.html>`.  See
`sbpy`'s installation instructions for details.


Classes
-------
solar_fluxd   - Solar spectrum in flux density.

"""

import os
from astropy.utils.exceptions import AstropyWarning
from warnings import warn
from .. import bib


__all__ = [
    'solar_fluxd',
]


# Load solar spectra, if possible
try:
    import synphot as S
    _synphot = True
except ImportError as e:
    warn(AstropyWarning('synphot is required for solar_fluxd'))
    _synphot = False

try:
    path = os.sep.join((os.environ['PYSYN_CDBS'], 'grid', 'k93models',
                       'standards'))
    _cdbs = True
except KeyError as e:
    warn(AstropyWarning('PYSYN_CDBS environment variable is required in order to locate solar spectra.'))
    _cdbs = False

_sun = {}

if _synphot and _cdbs:
    fn = os.sep.join((path, 'sun_kurucz93.fits'))
    not_found = []
    try:
        _sun['K93'] = (S.SourceSpectrum.from_file(fn),
                       ('1993KurCD..13.....K', '1996AJ....112..307C'))
    except FileNotFoundError as e:
        not_found.append(fn)

    fn = os.sep.join((path, 'sun_castelli.fits'))
    try:
        _sun['C96'] = (S.SourceSpectrum.from_file(fn), '1996AJ....112..307C')
    except FileNotFoundError as e:
        not_found.append(fn)

    if len(_sun) == 0:
        warn(AstropyWarning('No solar spectra found, tried: ', not_found))

    del fn, not_found

def solar_fluxd(wave, unit='W / (m2 um)', source='K93'):
    """Model spectrum of the Sun.

    Parameters
    ----------
    wave : `~astropy.units.Quantity`
      The wavelengths of the returned spectrum.

    unit : string, `~astropy.units.Unit`, optional
      Spectral units.

    source : string
      A key indicating the requested source of the spectrum.  See
      References below for keys and references.


    Returns
    -------
    fluxd : `~astropy.units.Quantity`
      Spectrum of the Sun, binned to match `wave`.


    References
    ----------
    [K93] Kurucz (1993 Smithsonian Astropys. Obs. CD-ROM No. 13) solar
    spectrum from UV to mid-IR at low resolution, scaled by Colina,
    Bohlin & Castelli (1996, AJ 112, 307).

    [C94] Colina et al. (1996, AJ 112, 307) solar spectrum at 1-Å
    resolution from UV to near-IR (0.12 to 2.5 μm).

    """
    
    if not (_synphot and _cdbs and len(_sun)):
        warn(AstropyWarning('Requirements for `sbpy.activty.sun` not met.'))
        return None

    try:
        sun, bibcode = _sun[source]
    except KeyError as e:
        raise KeyError('Solar spectrum key {} not found.  {} spectra loaded: {}'.format(source, len(_sun), ', '.join(_sun.keys())))

    bib.register('data.sun.solar_fluxd', bibcode)

    # Method adapted from http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
    specele = S.SpectralElement(S.ConstFlux1D(1))
    obs = S.Observation(sun, specele, binset=wave, force='extrap')
    return obs.sample_binned(wave, unit)
