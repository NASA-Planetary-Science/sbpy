# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=====================
SBPy Vega Core Module
=====================
"""

__all__ = [
    'Vega'
]

import os
import astropy.units as u
from astropy.utils.state import ScienceState
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from ...spectroscopy.sources import SpectralStandard
from . import sources

try:
    import synphot
except ImportError:
    synphot = None

__doctest_requires__ = {'Vega': 'synphot'}


class Vega(SpectralStandard):
    """Vega spectrum.

    Parameters
    ----------
    wave : `~astropy.units.Quantity`
        Spectral wavelengths.

    fluxd : `~astropy.units.Quantity`
        Spectral flux density.

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
    Get the default Vega spectrum:
    >>> vega = Vega.from_default()               # doctest: +IGNORE_OUTPUT

    Create Vega from a file:
    >>> vega = Vega.from_file('filename')        # doctest: +SKIP

    Evaluate Vega at 1 μm (interpolation):
    >>> print(vega(1 * u.um))                    # doctest: +FLOAT_CMP
    6.326492514857613e-09 W / (m2 um)

    Observe Vega through as if through a spectrometer:
    >>> wave = np.linspace(0.4, 0.6) * u.um
    >>> spec = vega.observe(wave)

    Observe Vega through a filter:
    >>> V = sbpy.util.get_bandpass('Johnson V')
    >>> fluxd = vega.observe(V)

    User provided calibration:
    >>> from sbpy.calib import vega_fluxd
    >>> vega = Vega(None)
    >>> with vega_fluxd.set({'V': 3674 * u.Jy}):
    >>>     print(vega.observe('V'))
    3674.0 Jy

    """

    def __repr__(self):
        if self.description is None:
            return '<Vega>'
        else:
            return '<Vega: {}>'.format(self.description)

    @classmethod
    def from_builtin(cls, name):
        """Vega spectrum from a built-in `sbpy` source.

        Parameters
        ----------
        name : string
          The name of a Vega spectrum parameter set in
          `sbpy.spectroscopy.vega.sources`.

        """

        from astropy.utils.data import _is_url

        try:
            parameters = getattr(sources, name).copy()

            if not _is_url(parameters['filename']):
                # find in the module's location
                parameters['filename'] = get_pkg_data_filename(
                    os.path.join('data', parameters['filename']))

            vega = Vega.from_file(**parameters)
        except AttributeError:
            msg = 'Unknown Vega spectrum "{}".  Valid spectra:\n{}'.format(
                name, sources.available)
            raise ValueError(msg)

        return vega

    @classmethod
    def from_default(cls):
        """``Vega`` object with the `sbpy` default Vega spectrum.

        The spectrum will be ``None`` if `synphot` is not available.

        """
        from .. import vega_spectrum
        if synphot:
            vega = vega_spectrum.get()
        else:
            vega = cls(None)
        return vega

    @staticmethod
    def show_builtin():
        """List built-in Vega spectra."""
        rows = []
        for name in sources.available:
            source = getattr(sources, name)
            rows.append((name, source['description']))
        Table(rows=rows, names=('name', 'description')).pprint(
            max_lines=-1, max_width=-1)

    def filt(self, bp, unit='W / (m2 um)', **kwargs):
        """Observe Vega through a single filter.

        Uses the `sbpy` calibration system.  If `bp` has been set in
        `~sbpy.calib.vega_fluxd`, then those values will be returned.
        Otherwise, the Vega spectrum defined by this object will be
        filtered with the provided bandpass.


        Parameters
        ----------
        bp: string or `~synphot.SpectralElement`
            The name of a filter, or a transmission spectrum as a
            `~synphot.SpectralElement`.  See
            `sbpy.spectroscopy.sources.SpectralSource` for built-in
            filter names.

        unit: string, `~astropy.units.Unit`, optional
            Spectral flux density units of the output.

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        lambda_eff: `~astropy.units.Quantity`
            Effective wavelength.  ``None`` if it cannot be calculated.

        lambda_pivot: `~astropy.units.Quantity`
            Pivot wavelength.  ``None`` if it cannot be calculated.

        fluxd: `~astropy.units.Quantity`
            Spectral flux density.

        """

        from ...calib import vega_fluxd

        source_fluxd = vega_fluxd.get()
        return super().filt(bp, unit=unit, source_fluxd=source_fluxd,
                            **kwargs)
