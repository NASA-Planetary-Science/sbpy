# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy Coma Spectroscopy

Spectrophotometric classes for cometary comae.

Uses the astropy modeling framework.

Requires synphot.

"""

import numpy as np
import astropy.units as u
from astropy.modeling import Fittable1DModel, Parameter
from .core import SpectralGradient
from .sources import BlackbodySource
from ..activity.dust import Afrho, Efrho
from .. import data as sbd
from .. import units as sbu
from ..calib import Sun

__all__ = ["Scattered", "Thermal"]


class Scattered(Fittable1DModel):
    """Light scattered from coma dust.

    Does not directly account for geometric albedo or phase function:

        F_lambda = cross_section * F_sun * R_lambda(S) / pi / rh**2 / delta**2

    where ``R`` is the spectral reddening function, ``F_sun`` is the flux
    density of the Sun at 1 au, ``rh`` is the heliocentric distance scale factor
    for the solar flux density, ``delta`` is the observer-coma distance, and
    ``cross_section`` is the effective cross sectional area of the dust.


    Parameters
    ----------
    eph: dictionary-like, `~sbpy.data.Ephem`
        Ephemerides of the comet.  Required fields: 'rh', 'delta'.

    unit : string or `astropy.units.Unit`
        Unit of the resultant spectrum (spectral flux density).

    cross_section : float or `~astropy.units.Quantity`, optional
        Cross-sectional area of the dust.  For float input, km**2 is assumed.

    S : float, `~astropy.units.Quantity`, or
    `~sbpy.spectroscopy.SpectralGradient`, optional
        Spectral gradient.  For float input, %/100 nm is assumed.  For float or
        quantity input, the wavelength of normalization is 550 nm.

    """

    n_inputs = 1
    n_outputs = 1

    input_units = {"x": u.um}
    _input_units_allow_dimensionless = True
    input_units_equivalencies = {"x": u.spectral()}

    cross_section = Parameter(
        description="cross-sectional area of dust", default=1, unit=u.km**2
    )
    S = Parameter(
        description="spectral gradient", default=5, unit=u.percent / sbu.hundred_nm
    )

    @sbd.dataclass_input(eph=sbd.Ephem)
    def __init__(self, eph, unit, *args, **kwargs):
        self.sun = Sun.from_default()
        self.eph = eph
        self.unit = u.Unit(unit)
        if "S" in kwargs and isinstance(kwargs["S"], SpectralGradient):
            kwargs["S"] = u.Quantity(
                kwargs["S"].renormalize(550 * u.nm), u.percent / sbu.hundred_nm
            )
        super().__init__(*args, **kwargs)

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "cross_section": u.km**2,
            "S": u.percent / sbu.hundred_nm,
        }

    @property
    def rh(self):
        return self.eph["rh"][0]

    @property
    def delta(self):
        return self.eph["delta"][0]

    def evaluate(self, wave, cross_section, S):
        _wave = u.Quantity(wave, u.um)
        _cross_section = u.Quantity(cross_section, u.km**2)
        # Keep S as a scalar value
        _S = np.atleast_1d(u.Quantity(S, u.percent / sbu.hundred_nm))[0]

        sun = self.sun.redden(SpectralGradient(_S, wave0=550 * u.nm))
        if np.size(_wave) == 1:
            F_sun = sun(_wave, unit=self.unit)
        else:
            F_sun = sun.observe(_wave, unit=self.unit)

        spec = (
            _cross_section
            * F_sun
            / np.pi
            / self.delta**2
            / self.rh.to_value(u.au) ** 2
        ).to(self.unit)
        return spec if hasattr(cross_section, "unit") else spec.value

    # def fit_deriv(self, wave, cross_section, S):
    #     spec = self.evaluate(wave, cross_section, S)
    #     dspec_dcross_section = spec / cross_section
    #     dspec_dS = spec * (wave - 0.55 * u.um)
    #     return [dspec_dcross_section, dspec_dS]

    @u.quantity_input(wave=u.m)
    def afrho(self, wave, aper):
        """Convert spectrum to Afρ quantity.


        Parameters
        ----------
        wave : `~astropy.units.Quantity`
            Spectral wavelengths.

        aper : `~astropy.units.Quantity` or `~sbpy.activity.Aperture`
            Aperture of the observation as a circular radius (length
            or angular units), or as an `~sbpy.activity.Aperture`.

        """

        return Afrho.from_fluxd(wave, self(wave), aper, self.eph)


class Thermal(Fittable1DModel):
    """Thermal emission from coma dust.

    Simple blackbody model, does not directly account for thermal emissivity.

        F_lambda = cross_section * B_lambda(T) / delta**2

    where ``B_lambda`` is the Planck function, ``T = Tscale * 278 / sqrt(rh)``
    is the dust temperature, and ``cross_section`` is the effective cross
    sectional area.


    Parameters
    ----------
    eph: dictionary-like, `~sbpy.data.Ephem`
        Ephemerides of the comet.  Required fields: 'rh', 'delta'

    unit : string or `astropy.units.Unit`
        Unit of the resultant spectrum (spectral flux density).

    cross_section : float or `~astropy.units.Quantity`, optional
        Effective cross-sectional area of the dust.  For float input, km**2 is
        assumed.

    Tscale : float or `~astropy.units.Quantity`, optional
        LTE blackbody sphere temperature scale factor: ``T = Tscale * 278 /
        sqrt(rh)``.

    """

    n_inputs = 1
    n_outputs = 1

    input_units = {"x": u.um}
    _input_units_allow_dimensionless = True
    input_units_equivalencies = {"x": u.spectral()}

    cross_section = Parameter(
        description="effective emission cross-sectional area of dust",
        default=1,
        unit=u.km**2,
    )
    Tscale = Parameter(
        description="LTE blackbody sphere temperature scale factor",
        default=1.0,
        min=0,
    )

    @sbd.dataclass_input(eph=sbd.Ephem)
    def __init__(self, eph, unit, *args, **kwargs):
        self.eph = eph
        self.unit = unit
        super().__init__(*args, **kwargs)

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "cross_section": u.km**2,
            "Tscale": "",
        }

    @property
    def rh(self):
        return self.eph["rh"][0]

    @property
    def delta(self):
        return self.eph["delta"][0]

    def evaluate(self, wave, cross_section, Tscale):
        """Evaluate the model.


        Parameters
        ----------
        x : float, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Frequency or wavelength at which to compute the blackbody. If no
            units are given, then micrometer is assumed.

        cross_section : float or `~astropy.units.Quantity`, optional
            Effective cross-sectional area of the dust.  If no units are given,
            then km**2 is assumed.

        Tscale : float or `~astropy.units.Quantity`, optional
            LTE blackbody sphere temperature scale factor: ``T = Tscale * 278 /
            sqrt(rh)``.


        Returns
        -------
        y : number, `~numpy.ndarray`, or `~astropy.units.Quantity`
            Spectrum.

        """

        _wave = u.Quantity(wave, u.um)
        _cross_section = u.Quantity(cross_section, u.km**2)
        T = Tscale * u.Quantity(278, u.K) / np.sqrt(self.rh.to_value("au"))

        B = BlackbodySource(T)
        if np.size(_wave) == 1:
            F_bb = B(_wave, unit=self.unit)
        else:
            F_bb = B.observe(_wave, unit=self.unit)

        spec = (_cross_section * (F_bb / np.pi) / self.delta**2).to(self.unit)
        return spec if hasattr(cross_section, "unit") else spec.value

    def efrho(self, epsilon, aper):
        """Convert to εfρ quantity.


        Parameters
        ----------
        epsilon : float
            Grain emissivity.

        aper : `~astropy.units.Quantity` or `~sbpy.activity.Aperture`
            Aperture of the observation as a circular radius (length or angular
            units), or as an `~sbpy.activity.Aperture` object.

        """

        return Efrho.from_cross_section(self.cross_section, epsilon, aper, self.delta)
