# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
SBPy Sun Core Module
====================
"""

from astropy.utils.state import ScienceState
from ..core import SpectralStandard

__all__ = [
    'Sun',
    'default_sun'
]

__doctest_requires__ = {'Sun': 'synphot'}


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

    def __repr__(self):
        if self.description is None:
            return '<Sun>'
        else:
            return '<Sun: {}>'.format(self.description)

    @classmethod
    def from_builtin(cls, name):
        """Solar spectrum from a built-in `sbpy` source.

        Parameters
        ----------
        name : string
          The name of a solar spectrum parameter set in
          `sbpy.spectroscopy.sun.sources`.

        """

        import os
        from astropy.utils.data import _is_url
        from . import sources

        try:
            parameters = getattr(sources, name).copy()

            if not _is_url(parameters['filename']):
                # find in the module's location
                path = os.path.dirname(__file__)
                parameters['filename'] = os.sep.join(
                    (path, parameters['filename']))

            return Sun.from_file(**parameters)
        except AttributeError:
            msg = 'Unknown solar spectrum "{}".  Valid spectra:\n{}'.format(
                name, sources.available)
            raise ValueError(msg)


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

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            return Sun.from_builtin(value)
        elif isinstance(value, Sun):
            return value
        else:
            raise TypeError("default_sun must be a string or Sun instance.")
