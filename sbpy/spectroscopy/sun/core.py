# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
SBPy Sun Core Module
====================
"""

try:
    import synphot
except ImportError as e:
    from warnings import warn
    from astropy.utils.exceptions import AstropyWarning
    warn(AstropyWarning('synphot is required for the sun module.'))

from astropy.utils.state import ScienceState
from ..core import SpectralStandard

__all__ = [
    'Sun',
    'default_sun'
]

class Sun(SpectralStandard):
    """Solar spectrum.

    Parameters
    ----------
    wave : `~astropy.units.Quantity`
      The spectral wavelengths.

    fluxd : `~astropy.units.Quantity`
      The solar flux densities, at 1 au.

    description : string, optional
      A brief description of the source spectrum.

    bibcode : string, optional
      Bibliography code for `sbpy.bib.register`.

    meta : dict, optional
      Any additional meta data, passed on to
      `~synphot.SourceSpectrum`.


    Attributes
    ----------
    wave        - Wavelengths of the source spectrum.
    fluxd       - Source spectrum.
    description - Brief description of the source spectrum.
    meta        - Meta data.


    Examples
    --------
    Create solar standard from `synphot.SourceSpectrum`::
      >>> import astropy.constants as const
      >>> import synphot
      >>> source = synphot.SourceSpectrum(synphot.BlackBody1D, temperature=5770)
      >>> source *= (3.14159 * const.R_sun**2 / const.au**2).decompose().value
      >>> sun = Sun(source, description='5770 K blackbody Sun')

    Create solar standard from an array::
      >>> import numpy as np
      >>> import astropy.units as u
      >>> import astropy.constants as const
      >>> from astropy.modeling.blackbody import blackbody_lambda
      >>> wave = np.logspace(-1, 2) * u.um
      >>> fluxd = blackbody_lambda(wave, 5770 * u.K) * np.pi * u.sr
      >>> fluxd *= (const.R_sun**2 / const.au**2).decompose().value
      >>> sun = Sun.from_array(wave, fluxd, description='5770 K blackbody Sun')

    Create solar standard from a file::
      >>> sun = Sun.from_file('filename')        # doctest: +SKIP

    Evaluate solar standard at 1 μm::
      >>> print(sun(1 * u.um))                   # doctest: +SKIP

    """
    pass

class default_sun(ScienceState):
    """The default solar spectrum to use.

    To retrieve the current default::

      >>> sun = default_sun.get()

    To change it::

      >>> from sbpy.spectroscopy.sun import default_sun
      >>> with default_sun.set('E490_2014'):
      ...     # E490_2014 in effect
      ...     pass

    Or, you may use a string::

      >>> with default_sun.set('E490_2014LR'):
      ...     # E490_2014LR in effect
      ...     pass

    """
    _value = 'E490_2014'

    @staticmethod
    def get_sun_from_string(arg):
        """Return a Sun instance from a string."""
        import os
        import sys
        from astropy.utils.data import _is_url
        from . import sources
        
        try:
            parameters = getattr(sources, arg).copy()
            filename = parameters.pop('file')

            if not _is_url(filename):
                path = os.path.dirname(__file__)
                filename = os.sep.join((path, filename))
                
            sun = Sun.from_file(filename, **parameters)
        except AttributeError:
            msg = 'Unknown solar spectrum "{}".  Valid spectra:\n{}'.format(arg, sources.available)
            raise ValueError(msg)

        return sun

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return cls.get_sun_from_string(value)
        elif isinstance(value, Sun):
            return value
        else:
            raise TypeError("default_sun must be a string or Sun instance.")

