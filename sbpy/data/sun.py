# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
===============
SBPy Sun Module
===============

This module requires `synphot` and solar spectra available in Space
Telescope Science Institute's `Calibration Reference Data System
<http://www.stsci.edu/hst/observatory/crds/throughput.html>`.  See
`synphot`'s installation instructions for details.


Functions
---------
solar_filt    - Observe a solar spectrum through a filter.
solar_fluxd   - Solar spectrum in flux density.

"""

import os
from warnings import warn
import numpy as np
from astropy.utils.exceptions import AstropyWarning
from .. import bib


__all__ = [
    'solar_fluxd',
    'solar_filt',
]


# Load solar spectra, if possible
try:
    import synphot
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
        # special handling because file includes NaNs
        #_sun['K93'] = (synphot.SourceSpectrum.from_file(fn),
        #               ('1993KurCD..13.....K', '1996AJ....112..307C'))
        header, wavelengths, fluxes = synphot.specio.read_spec(fn)
        i = ~np.isnan(fluxes)
        _sun['K93'] = (
            synphot.SourceSpectrum(
                synphot.Empirical1D, points=wavelengths[i],
                lookup_table=fluxes[i], keep_neg=False,
                meta={'header': header}),
            ('1993KurCD..13.....K', '1996AJ....112..307C'))
    except FileNotFoundError as e:
        not_found.append(fn)

    fn = os.sep.join((path, 'sun_castelli.fits'))
    try:
        _sun['C96'] = (synphot.SourceSpectrum.from_file(fn), '1996AJ....112..307C')
    except FileNotFoundError as e:
        not_found.append(fn)

    if len(_sun) == 0:
        warn(AstropyWarning('No solar spectra found, tried: ', not_found))

    fn = os.sep.join((os.environ['PYSYN_CDBS'], 'calspec',
                      'alpha_lyr_stis_008.fits'))
    try:
        _vega = (synphot.SourceSpectrum.from_file(fn), '2014AJ....147..127B')
    except FileNotFoundError as e:
        warn(AstropyWarning('Vega spectrum not found, tried: ', fn))
        _vega = None
        
    del fn, not_found, path

def _solar_source(source):
    """Retrieve solar spectrum and register in bibliography.

    Parameters
    ----------
    source : string
      Source key (see `solar_fluxd`).

    Returns
    -------
    sun : `~synphot.SourceSpectrum`
      The requested spectrum of the Sun.

    """
    
    if not (_synphot and _cdbs and len(_sun)):
        warn(AstropyWarning('Requirements for `sbpy.data.sun` not met.'))
        return None

    try:
        sun, bibcode = _sun[source]
    except KeyError as e:
        raise KeyError('Solar spectrum key {} not found.  {} spectra loaded: {}'.format(source, len(_sun), ', '.join(_sun.keys())))

    bib.register('data.sun', bibcode)
    bib.register('data.sun', 'synphot')

    return sun

def solar_fluxd(wave, unit='W / (m2 um)', source='K93'):
    """Model spectrum of the Sun.

    Parameters
    ----------
    wave : `~astropy.units.Quantity` or `None`
      The wavelengths of the returned spectrum.  Set to `None` to
      return the original wavelengths of the source spectrum.

    unit : string, `~astropy.units.Unit`, optional
      Spectral units of the output: flux density, 'vegamag', 'ABmag',
      or 'STmag'.

    source : string or `~synphot.SourceSpectrum`
      A key indicating the requested source of the spectrum (see
      References) or a `SourceSpectrum` object.


    Returns
    -------
    wave : `~astropy.units.Unit`, optional
      Returned if the original wavelengths of the spectrum were
      requested (see `wave` input parameter).

    fluxd : `~astropy.units.Quantity`
      Spectrum of the Sun.  If multiple wavelengths are requested, the
      solar spectrum will be binned to match `wave`.  Otherwise, the
      spectrum at that wavelength will be returned.


    Notes
    -----
    Spectral data is from STScI's Calibration Reference Data System.
    See `synphot` installation instructions for details.


    References
    ----------
    [K93] Kurucz (1993 Smithsonian Astropys. Obs. CD-ROM No. 13) solar
    spectrum from UV to mid-IR at low resolution, scaled by Colina,
    Bohlin & Castelli (1996, AJ 112, 307).

    [C96] Colina et al. (1996, AJ 112, 307) solar spectrum at 1-Å
    resolution from UV to near-IR (0.12 to 2.5 μm).

    """

    if not (_synphot and _cdbs and len(_sun)):
        warn(AstropyWarning('Requirements for `sbpy.data.sun` not met.'))
        return None

    sun = _solar_source(source)

    if wave is None:
        return sun.waveset, sun(sun.waveset, unit)

    if np.size(wave) > 1:
        # Method adapted from http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
        specele = synphot.SpectralElement(synphot.ConstFlux1D(1))
        obs = synphot.Observation(sun, specele, binset=wave, force='extrap')
        fluxd = obs.sample_binned(wave, unit)
    else:
        fluxd = sun(wave, unit)

    return fluxd

def solar_filt(bp, unit='W / (m2 um)', source='K93', **kwargs):
    """Solar spectrum observed through a filter.

    Parameters
    ----------
    bp : string or `~synphot.spectrum.SpectralElement`
      The name of a filter, or a transmission spectrum.  See notes for
      acceptable filter names.

    unit : string, `~astropy.units.Unit`, optional
      Spectral units of the output.

    source : string
      A key indicating the requested source of the spectrum.  See
      `solar_fluxd` for details.

    **kwargs
      Additional keyword arguments for
      `~synphot.observation.Observation`.

    Returns
    -------
    wave : `~astropy.units.Quantity`
      Effective wavelength.

    fluxd : `~astropy.units.Quantity`
      The flux density.


    Notes
    -----
    Filter reference data is from STScI's Calibration Reference Data
    System.  See `synphot` installation instructions for details.

    * ``'bessel_j'`` (Bessel *J*)
    * ``'bessel_h'`` (Bessel *H*)
    * ``'bessel_k'`` (Bessel *K*)
    * ``'cousins_r'`` (Cousins *R*)
    * ``'cousins_i'`` (Cousins *I*)
    * ``'johnson_u'`` (Johnson *U*)
    * ``'johnson_b'`` (Johnson *B*)
    * ``'johnson_v'`` (Johnson *V*)
    * ``'johnson_r'`` (Johnson *R*)
    * ``'johnson_i'`` (Johnson *I*)
    * ``'johnson_j'`` (Johnson *J*)
    * ``'johnson_k'`` (Johnson *K*)

    """

    if not (_synphot and _cdbs and len(_sun)):
        warn(AstropyWarning('Requirements for `sbpy.data.sun` not met.'))
        return None

    assert isinstance(bp, (str, synphot.SpectralElement))
    
    if isinstance(bp, str):
        bp = synphot.SpectralElement.from_filter(bp)
    
    sun = _solar_source(source)
    obs = synphot.Observation(sun, bp, **kwargs)

    opts = {}
    if unit == 'vegamag':
        if _vega is None:
            warn(AstropyWarning('Requirements for `sbpy.data.sun` not met.'))
            return None

        opts['vegaspec'] = _vega[0]
        bib.register('sbpy.data.sun.solar_filt', _vega[1])

    return obs.effective_wavelength(), obs.effstim(unit, **opts)
