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

    """

    @dataclass_input
    def __init__(
        self,
        radius: u.physical.length,
        # phys: Phys | None = None,
    ):
        self.radius: u.Quantity = radius
        # self.phys: Phys = Phys() if phys is None else phys

    def to_faceted_model(self):
        raise NotImplemented

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
        """Integrate a function of $(i, e, \phi)$ over the sphere.


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

    @dataclass_input
    def absorption(
        self,
        wfb: u.Quantity | str,
        epsilon: u.physical.dimensionless,
        surface: Surface,
        eph: Ephem,
        unit: u.Unit | str = "W / um",
        interpolate: bool = False,
        **kwargs,
    ):
        """Total absorbed sunlight.

        Uses the current default ``Sun``.  See :ref:`_sbpy-calib` and
        :ref:`_default-spectra` for more information.


        Parameters
        ----------

        surface : `~sbpy.surface.surface.Surface`

        interpolate : bool, optional
            Set to ``True`` to interpolate the solar spectrum for ``wfb``.  See
            :ref:`_calib-binning-vs-interpolation` for guidance.

        kwargs : dict, optional
            Additional keyword arguments are passed to the integrator,
            `~scipy.integrate.dblquad`.

        """

        unit = u.Unit(unit)
        S = Shape._incident_sunlight(wfb, eph, unit / u.m**2, interpolate)

        def f(i, e, phi):
            return surface.absorptance(1, i)

        # set phase to zero: absorption does not depend on an observer
        A, err = self.integrate_i_e_phi(f, 0 * u.deg, **kwargs)
        return (S * epsilon * A).to(unit), (S * epsilon * err).to(unit)

    # def reflectance(self, albedo: float):
    #     def f(i, e, phi):
    #         return self.surface.absorptance(epsilon, i)

    #     return self.integrate_i_e_phi(f, 0 * u.deg)
