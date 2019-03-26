# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
====================
SBPy Sun Core Module
====================
"""

import os
import astropy.units as u
from astropy.utils.state import ScienceState
from astropy.utils.data import get_pkg_data_filename
from ..sources import SpectralStandard
from . import sources as sources

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
        Spectral wavelengths.

    fluxd : `~astropy.units.Quantity`
        Solar spectral flux density at 1 au.

    description : string, optional
        Brief description of the source spectrum.

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
    Get the default solar spectrum:
    >>> sun = Sun.from_default()

    Create solar standard from `synphot.SourceSpectrum`:
    >>> import astropy.constants as const
    >>> import synphot
    >>> source = (synphot.SourceSpectrum(synphot.BlackBody1D, temperature=5770)
    ...           (3.14159 * const.R_sun**2 / const.au**2).decompose())
    >>> sun = Sun(source, description='5770 K blackbody Sun')

    Create solar standard from an array:
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import astropy.constants as const
    >>> from astropy.modeling.blackbody import blackbody_lambda
    >>> wave = np.logspace(-1, 2) * u.um
    >>> fluxd = (blackbody_lambda(wave, 5770 * u.K) * np.pi * u.sr
    ...          (const.R_sun**2 / const.au**2).decompose())
    >>> sun = Sun.from_array(wave, fluxd, description='5770 K blackbody Sun')

    Create solar standard from a file:
    >>> sun = Sun.from_file('filename')        # doctest: +SKIP

    Interpolate to 0.62 μm:
    >>> sun(0.62 * u.um)                       # doctest: +FLOAT_CMP
    <Quantity 1720.5108871 W / (m2 um)>

    Observe as through a spectrometer:
    >>> import numpy as np
    >>> import astropy.units as u
    >>> sun = Sun.from_default()
    >>> wave = np.linspace(1, 2.5) * u.um
    >>> fluxd = sun.observe(wave)              # doctest: +IGNORE_OUTPUT

    Observe as through a filter:
    >>> sun = Sun.from_default()
    >>> sun.observe('johnson_v')               # doctest: +FLOAT_CMP
    <Quantity [1839.93273227] W / (m2 um)>

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

        from astropy.utils.data import _is_url

        try:
            parameters = getattr(sources, name).copy()

            if not _is_url(parameters['filename']):
                # find in the module's location
                parameters['filename'] = get_pkg_data_filename(
                    os.path.join('data', parameters['filename']))

            return Sun.from_file(**parameters)
        except AttributeError:
            msg = 'Unknown solar spectrum "{}".  Valid spectra:\n{}'.format(
                name, sources.available)
            raise ValueError(msg)

    @classmethod
    def from_default(cls):
        """Return the `sbpy` default solar spectrum."""
        return default_sun.get()

    @staticmethod
    def show_builtin():
        """List built-in solar spectra."""
        from astropy.table import Table
        rows = []
        for name in sources.available:
            source = getattr(sources, name)
            rows.append((name, source['description']))
        Table(rows=rows, names=('name', 'description')).pprint(
            max_lines=-1, max_width=-1)


class default_sun(ScienceState):
    """Get/set the `sbpy` default solar spectrum.

    To retrieve the current default:

    >>> sun = default_sun.get()

    To change it:

    >>> from sbpy.spectroscopy import default_sun
    >>> with default_sun.set('E490_2014'):
    ...     # E490_2014 in effect
    ...     pass

    Or, you may use a string:

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
