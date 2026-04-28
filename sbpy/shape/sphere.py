# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Asteroid and comet nucleus shape models
"""

__all__ = ["Sphere"]

from typing import Callable

import numpy as np
from numpy import pi
import astropy.units as u
from scipy.integrate import dblquad

from ..data.ephem import Ephem
from ..data.decorators import dataclass_input
from ..surfaces.surface import Surface
from .core import Shape
from .transformations import twovec


class Sphere(Shape):
    """A spherical object.


    Parameters
    ----------
    radius : |Quantity|
        The radius of the sphere.

    pole : tuple of |Quantity|, optional
        The direction of the pole as longitude and latitude.

    """

    @u.quantity_input
    @dataclass_input
    def __init__(
        self,
        radius: u.Quantity["length"],
        pole: tuple[u.Quantity["angle"]] | None = None,
    ):
        self.radius: u.Quantity = radius
        self.pole = None if pole is None else (pole[0], pole[1])

    @dataclass_input
    def integrate(
        self,
        func: Callable,
        **kwargs,
    ) -> tuple[u.Quantity, u.Quantity]:
        """Integrate arbitrary function over the surface of the sphere.


        Parameters
        ----------
        func : callable
            The function to integrate: ``func(phi, theta, *args)``.

        kwargs : dict, optional
            Keyword arguments passed to the integrator,
            `~scipy.integrate.dblquad`.


        Returns
        -------
        total : |Quantity|
            The result.

        err : |Quantity|
            Estimated integration error on ``total``.

        """

        result, err = dblquad(func, 0, pi, 0, 2 * pi, **kwargs) * self.radius**2

        return result, err

    def integrate_i_e_phi(
        self,
        func: Callable,
        phase: u.physical.angle,
        **kwargs,
    ) -> tuple[u.Quantity, u.Quantity]:
        """Integrate a function of $(i, e, \\phi)$ over the sphere.


        Parameters
        ----------
        func : callable
            The function to integrate: ``func(i, e, phase, *args)``,
            where :math:`i` is the angle of incidence, :math:`e` is the the
            angle of emittance, and :math:`phase` is the phase angle.

        phase : |Quantity|
            Sun-target-observer (phase) angle.

        kwargs : dict, optional
            Additional keyword arguments are passed to the integrator,
            `~scipy.integrate.dblquad`.


        Returns
        -------
        total : |Quantity|
            The result.

        err : |Quantity|
            Estimated integration error on ``total``.

        """

        # Let target-observer vector be +z
        r_obs = np.array([0, 0, 1])

        # Define target-Sun vector using the phase angle
        r_sun = np.array([np.sin(phase), 0, np.cos(phase)])

        # Define a matrix that transforms vectors in a reference frame with
        # target-Sun along +z
        if np.isclose(phase, 0):
            # if r_obs == r_sun, then use the identity matrix
            M = np.eye(3)
        else:
            M = twovec(r_sun, 2, r_obs, 0)

        def f(phi, theta, *args):
            """
            z-axis is to the observer
            theta is the angle from the z-axis
            phi is the angle from the x-axis
            """

            r = np.array(
                [
                    np.cos(phi) * np.sin(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(theta),
                ]
            )

            # angle of emittance
            e = theta * u.rad

            # transform to Sun-centered frame to find angle of incidence
            rprime = np.dot(M, r)
            i = np.arctan2(np.hypot(rprime[0], rprime[1]), rprime[2]) * u.rad

            a = func(i, e, phase, *args) * np.sin(theta)
            return a

        return self.integrate(f, **kwargs)

    @u.quantity_input
    @dataclass_input
    def absorption(
        self,
        I: u.Quantity,
        epsilon: u.Quantity["dimensionless"],
        surface: Surface,
        **kwargs,
    ) -> tuple[u.Quantity, u.Quantity]:
        """Surface total absorbed light.


        Parameters
        ----------
        I : |Quantity|
            Incident light.

        epsilon : |Quantity| or float
            Surface emissivity.

        surface : `~sbpy.surface.surface.Surface`
            The surface being illuminated.

        kwargs : dict, optional
            Additional keyword arguments are passed to the integrator,
            `~scipy.integrate.dblquad`.


        Returns
        -------
        A : |Quantity|
            Surface total absorbed energy.

        err : |Quantity|
            Estimated surface integration error on ``A``.

        """

        def f(i, e, phi):
            return surface.absorptance(1, i)

        # set phase to zero: absorption does not depend on an observer
        A, err = self.integrate_i_e_phi(f, 0 * u.deg, **kwargs)
        return (I * epsilon * A), (I * epsilon * err)

    @dataclass_input
    def absorption_of_sunlight(
        self,
        wfb,
        epsilon: u.Quantity["dimensionless"] | float,
        surface: Surface,
        eph: Ephem,
        interpolate: bool = False,
        unit: u.Unit = "W / um",
        **kwargs,
    ) -> tuple[u.Quantity, u.Quantity]:
        """Surface total absorbed sunlight.

        Uses the current default ``Sun``.  See :ref:`_sbpy-calib` and
        :ref:`_default-spectra` for more information.


        Parameters
        ----------
        wfb : |Quantity|, |SpectralElement|, str, or list
            Wavelengths, frequencies, or bandpasses.

        epsilon : |Quantity| or float
            Surface emissivity.

        surface : `~sbpy.surface.surface.Surface`
            The surface being illuminated.

        eph : |Ephem|
            Distance to the sun.

        interpolate : bool, optional
            Set to ``True`` to interpolate the solar spectrum for ``wfb``.  See
            :ref:`_calib-binning-vs-interpolation` for guidance.

        unit : |Unit| or str, optional
            Spectral density units for the return value (e.g., W/μm).

        kwargs : dict, optional
            Additional keyword arguments are passed to the integrator,
            `~scipy.integrate.dblquad`.


        Returns
        -------
        A : |Quantity|
            Surface total absorbed sunlight at each ``wfb``.

        err : |Quantity|
            Estimated surface integration error on ``A``.

        """

        unit = u.Unit(unit)
        S = Shape._incident_sunlight(wfb, eph["rh"][0], unit / u.m**2, interpolate)
        A, err = self.absorption(S, epsilon, surface, **kwargs)
        return A.to(unit), err.to(unit)

    @dataclass_input
    def reflected_light(
        self,
        I: u.Quantity,
        albedo: u.Quantity["dimensionless"],
        surface: Surface,
        eph: Ephem,
        **kwargs,
    ) -> tuple[u.Quantity, u.Quantity]:
        """Surface total reflected light.


        Parameters
        ----------
        I : |Quantity|
            Incident light.

        albedo : |Quantity| or float
            Surface albedo.

        surface : `~sbpy.surface.surface.Surface`
            The surface being illuminated.

        eph : |Ephem|
            Target-observer distance and Sun-target-observer (phase) angle.

        kwargs : dict, optional
            Additional keyword arguments are passed to the integrator,
            `~scipy.integrate.dblquad`.


        Returns
        -------
        R : |Quantity|
            Reflected light at observer.

        err : |Quantity|
            Estimated surface integration error on ``R``.

        """

        delta = eph["delta"][0]
        phase = eph["phase"][0]

        def f(i, e, phi):
            return surface.reflectance(1, i, e, phi)

        R, err = self.integrate_i_e_phi(f, phase, **kwargs)
        return (
            (I * albedo * R / delta**2).to(I.unit),
            (I * albedo * err / delta**2).to(I.unit),
        )

    @dataclass_input
    def reflected_sunlight(
        self,
        wfb,
        albedo: u.Quantity["dimensionless"],
        surface: Surface,
        eph: Ephem,
        unit: u.Unit | str = "W / (m2 um)",
        interpolate: bool = False,
        **kwargs,
    ) -> tuple[u.Quantity, u.Quantity]:
        """Surface total reflected sunlight.

        Uses the current default ``Sun``.  See :ref:`_sbpy-calib` and
        :ref:`_default-spectra` for more information.


        Parameters
        ----------
        wfb : |Quantity|, |SpectralElement|, str, or list
            Wavelengths, frequencies, or bandpasses.

        albedo : |Quantity| or float
            Surface albedo.

        surface : `~sbpy.surface.surface.Surface`
            The surface being illuminated.

        eph : |Ephem|
            Heliocentric distance, target-observer distance and
            Sun-target-observer (phase) angle.

        interpolate : bool, optional
            Set to ``True`` to interpolate the solar spectrum for ``wfb``.  See
            :ref:`_calib-binning-vs-interpolation` for guidance.

        kwargs : dict, optional
            Additional keyword arguments are passed to the integrator,
            `~scipy.integrate.dblquad`.

        """

        unit = u.Unit(unit)
        S = Shape._incident_sunlight(wfb, eph["rh"][0], unit, interpolate)
        R, err = self.reflected_light(S, albedo, surface, eph, **kwargs)
        return R.to(unit), err.to(unit)
