# Licensed under a 3-clause BSD style license - see LICENSE.rst

from typing import Optional

import numpy as np
import astropy.units as u

from ..calib import Sun
from ..units.typing import SpectralRadianceQuantity, UnitLike
from ..spectroscopy.sources import SpectralSource, SinglePointSpectrumError
from .surface import Surface

__all__ = ["ScatteredLight", "ScatteredSunlight"]


class ScatteredLight:
    """Scatter light off a surface.


    Examples
    --------

    Scatter sunlight from a Lambertian surface.

    >>> import astropy.units as u
    >>> from sbpy.surfaces import ScatteredLight, LambertianSurface
    >>> from sbpy.calib import Sun
    >>>
    >>> sun = Sun.from_default()
    >>> surface = LambertianSurface()
    >>> scattered = ScatteredLight(surface, sun)
    >>>
    >>> wave = 550 * u.nm
    >>> albedo = 0.1
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>> scattered.radiance(wave, albedo, i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 25.8915765 W / (sr um m2)>

    """

    def __init__(self, surface: Surface, source: SpectralSource):
        self.surface = surface
        self.source = source

    @u.quantity_input
    def radiance(
        self,
        wfb,
        albedo: u.physical.dimensionless,
        i: u.physical.angle,
        e: u.physical.angle,
        phi: u.physical.angle,
        unit: Optional[UnitLike] = None,
    ) -> SpectralRadianceQuantity:
        """Spectral radiance scattered from the surface.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, list
            Wavelengths, frequencies, bandpass, or list of
            bandpasses of the observation.  Bandpasses require
            `~synphot`.

        albedo : `~astropy.units.Quantity`
            Surface albedo for given `wfb`

        i : `~astropy.units.Quantity`
            Angle from normal of incident light.

        e : `~astropy.units.Quantity`
            Angle from normal of emitted light.

        phi : `~astropy.units.Quantity`
            Source-target-observer (phase) angle.

        unit : str or `~astropy.units.Unit`, optional
            Spectral radiance units of the result.


        Returns
        -------
        I : `~astropy.units.Quantity`

        """

        try:
            F_i = self.source.observe(wfb, unit=unit)
        except SinglePointSpectrumError:
            F_i = self.source(wfb, unit=unit)

        return F_i * self.surface.reflectance(albedo, i, e, phi)

    def radiance_from_vectors(
        self,
        wfb,
        albedo: u.physical.dimensionless,
        n: np.ndarray,
        r: u.physical.length,
        ro: u.physical.length,
        unit: Optional[UnitLike] = None,
    ) -> SpectralRadianceQuantity:
        """Vector-based alternative to `radiance`.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, list
            Wavelengths, frequencies, bandpass, or list of
            bandpasses of the observation.  Bandpasses require
            `~synphot`.

        albedo : `~astropy.units.Quantity`
            Surface albedo for given `wfb`

        n : `numpy.ndarray`
            Surface normal vector.

        r : `~astropy.units.Quantity`
            Radial vector from the surface to the light source.

        ro : `~astropy.units.Quantity` or ``None``
            Radial vector from the surface to the observer.  Parameter is unused
            and may be ``None``.

        unit : str or `~astropy.units.Unit`, optional
            Spectral radiance units of the result.


        Returns
        -------
        I : `~astropy.units.Quantity`

        """

        i = self._angle(n, r)
        e = self._angle(n, ro)
        phi = self._angle(r, ro)
        return self.radiance(wfb, albedo, i, e, phi, unit=unit)


class ScatteredSunlight(ScatteredLight):
    """Scatter sunlight off a surface.


    Examples
    --------

    Scatter sunlight from a Lambertian surface.

    >>> import astropy.units as u
    >>> from sbpy.surfaces import ScatteredSunlight, LambertianSurface
    >>>
    >>> surface = LambertianSurface()
    >>> scattered = ScatteredSunlight(surface)
    >>>
    >>> wave = 550 * u.nm
    >>> albedo = 0.1
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>> scattered.radiance(wave, albedo, i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 25.8915765 W / (sr um m2)>

    The solar spectrum used is controlled with the ``sbpy.calib`` module::

    >>> from sbpy.calib import solar_spectrum
    >>> with solar_spectrum.set("E490_2014LR"):
    ...     scattered = ScatteredSunlight(surface)
    >>> scattered.radiance(wave, albedo, i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 25.8915765 W / (sr um m2)>

    """

    def __init__(self, surface: Surface):
        self.surface = surface
        self.source = Sun.from_default()
