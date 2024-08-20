# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy activity.dust.dynamics
===========================

Dust dynamical models.

"""

__all__ = [
    "DynamicalModel",
    "FreeExpansion",
    "SolarGravity",
    "SolarGravityAndRadiationPressure",
    "SolverFailed",
]

import abc
from typing import TypeVar

import numpy as np

try:
    from scipy.integrate import solve_ivp
except ImportError:  # pragma: no cover
    pass

from astropy.time import Time
import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame
import astropy.constants as const

from ..exceptions import SbpyException
from ..utils.decorators import requires
from .. import time  # noqa: F401
from .state import State


class SolverFailed(SbpyException):
    """DynamicalModel solver failed."""


FrameInputTypes = TypeVar("FrameInputTypes", str, BaseCoordinateFrame)


class DynamicalModel(abc.ABC):
    """Super-class for dynamical models.


    Parameters
    ----------
    **kwargs
        Arguments passed on to `~scipy.integrate.solve_ivp`.  Units are seconds,
        km, and km/s, e.g., ``max_step`` is a float value in units of seconds.
        For relative and absolute tolerance keywords, ``rtol`` and ``atol``,
        6-element arrays may be used, where the first three elements are for
        position, and the last three are for velocity.

    """

    @requires("scipy")
    def __init__(self, **kwargs):
        self.solver_kwargs: dict = dict(
            rtol=2.3e-14,
            jac=self.df_drv,
            method="LSODA",
        )
        self.solver_kwargs.update(kwargs)

    @classmethod
    @abc.abstractmethod
    def dx_dt(cls, t: float, rv: np.ndarray, *args) -> np.ndarray:
        """Derivative of position and velocity.


        Parameters
        ----------
        t : float
            Time, s.  Not used.

        rv : ndarray
            First three elements are the position vector at time ``t``, km. The
            next three elements are the velocity vector at time ``t``, km/s.

        *args :
            Additional parameters.


        Returns
        -------
        dx_dt : `numpy.ndarray`
            First three elements for :math:`dr/dt`, next three for :math:`dv/dt`.

        """

    @classmethod
    @abc.abstractmethod
    def df_drv(cls, t: float, rv: np.ndarray, *args) -> np.ndarray:
        """Jacobian matrix, :math:`df/drv`.


        Parameters
        ----------
        t : float
            Time, s.  Not used.

        rv : ndarray
            First three elements are the position vector at time ``t``, km. The
            next three elements are the velocity vector at time ``t``, km/s.

        *args :
            Additional parameters.


        Returns
        -------
        df_drv : `numpy.ndarray`
            First three elements for :math:`df/dr`, next three for
            :math:`df/dv`.

        """

    def solve(
        self,
        initial: State,
        t_final: Time,
        *args,
    ) -> State:
        """Solve the equations of motion for a single particle.

        The solution is calculated with `scipy.integrate.solve_ivp`.


        Parameters
        ----------
        initial : `State`
            Initial state (position and velocity at time) of the particle.

        t_final : `~astropy.time.Time` or `~astropy.units.Quantity`
            Time at which the solution is desired.  Use of ``Time`` versus
            ``Quantity`` must match how time is defined in the ``initial``
            state.

        *args :
            Additional arguments passed to `dx_dt` and `df_drv`.


        Returns
        -------
        final : State

        """

        t0: float
        t1: float
        if initial.arbitrary_time:
            t0 = initial.t.to_value("s")
            t1 = t_final.to_value("s")
        else:
            t0 = initial.t.et
            t1 = t_final.et

        result = solve_ivp(
            self.dx_dt,
            (t0, t1),
            initial.rv,
            args=args,
            **self.solver_kwargs,
        )

        if not result.success:
            raise SolverFailed(result.message)

        final: State = State(
            result.y[:3, -1] * u.km,
            result.y[3:, -1] * u.km / u.s,
            t_final,
            frame=initial.frame,
        )
        return final


class FreeExpansion(DynamicalModel):
    """Equation of motion solver for particle motion in free space.


    Parameters
    ----------
    **kwargs
        Arguments passed on to `~scipy.integrate.solve_ivp`.  Units are seconds,
        km, and km/s, e.g., ``max_step`` is a float value in units of seconds.
        For relative and absolute tolerance keywords, ``rtol`` and ``atol``,
        6-element arrays may be used, where the first three elements are for
        position, and the last three are for velocity.

    """

    @classmethod
    def dx_dt(cls, t: float, rv: np.ndarray, *args) -> np.ndarray:
        dx_dt = np.empty(6)
        dx_dt[:3] = rv[3:]
        dx_dt[3:] = 0

        return dx_dt

    @classmethod
    def df_drv(cls, t: float, rv: np.ndarray, *args) -> np.ndarray:
        # df_drv[i, j] = df_i/drv_j
        df_drv = np.zeros((6, 6))

        df_drv[0, 3] = 1
        df_drv[1, 4] = 1
        df_drv[2, 5] = 1

        return df_drv


class SolarGravity(DynamicalModel):
    """Equation of motion solver for a particle orbiting the Sun.


    Parameters
    ----------
    **kwargs
        Arguments passed on to `~scipy.integrate.solve_ivp`.  Units are seconds,
        km, and km/s, e.g., ``max_step`` is a float value in units of seconds.
        For relative and absolute tolerance keywords, ``rtol`` and ``atol``,
        6-element arrays may be used, where the first three elements are for
        position, and the last three are for velocity.

    """

    _GM: float = (const.G * const.M_sun).to_value("km3/s2")

    @property
    def GM(self):
        """Gravitational constant times mass."""
        return u.Quantity(self._GM, "km3/s2")

    @classmethod
    def dx_dt(cls, t: float, rv: np.ndarray, *args) -> np.ndarray:
        r = rv[:3]
        v = rv[3:]

        r2 = (r**2).sum()
        r1 = np.sqrt(r2)
        r3 = r2 * r1
        GM_r3 = cls._GM / r3

        dx_dt = np.empty(6)
        dx_dt[:3] = v
        dx_dt[3:] = -r * GM_r3

        return dx_dt

    @classmethod
    def df_drv(cls, t: float, rv: np.ndarray, *args) -> np.ndarray:
        r = rv[:3]
        r2 = (r**2).sum()
        r1 = np.sqrt(r2)
        GM_r5 = cls._GM / (r2 * r2 * r1)

        # df_drv[i, j] = df_i/drv_j
        df_drv = np.zeros((6, 6))

        df_drv[0, 3] = 1
        df_drv[1, 4] = 1
        df_drv[2, 5] = 1

        df_drv[3, 0] = GM_r5 * (r2 - 3 * r[0] * r[0])
        df_drv[3, 1] = -GM_r5 * 3 * r[0] * r[1]
        df_drv[3, 2] = -GM_r5 * 3 * r[0] * r[2]

        df_drv[4, 0] = -GM_r5 * 3 * r[1] * r[0]
        df_drv[4, 1] = GM_r5 * (r2 - 3 * r[1] * r[1])
        df_drv[4, 2] = -GM_r5 * 3 * r[1] * r[2]

        df_drv[5, 0] = -GM_r5 * 3 * r[2] * r[0]
        df_drv[5, 1] = -GM_r5 * 3 * r[2] * r[1]
        df_drv[5, 2] = GM_r5 * (r2 - 3 * r[2] * r[2])

        return df_drv


class SolarGravityAndRadiationPressure(DynamicalModel):
    """Equation of motion solver for a particle orbiting the Sun, including radiation force.


    Dust is parameterized with ``beta``, the ratio of the force from solar
    radiation pressure (:math:`F_r`) to that from solar gravity (:math:`F_g`):

    .. math::
        \\beta = \\frac{{F_r}}{{F_g}}

    For spherical dust grains, ``beta`` reduces to:

    .. math::
        \\beta = \\frac{{0.57 Q_{{pr}}}}{{\\rho a}}

    where :math:`Q_{{pr}}` is the radiation pressure efficiency averaged over
    the solar spectrum, :math:`\\rho` is the mass density of the grain (g/cm3),
    and :math:`a` is the grain radius (μm) (Burns et al. 1979).

    Only Newtonian gravity and radiation pressure are considered.
    Poynting-Roberston drag and general relativity are not included.


    Parameters
    ----------
    **kwargs
        Arguments passed on to `~scipy.integrate.solve_ivp`.  Units are seconds,
        km, and km/s, e.g., ``max_step`` is a float value in units of seconds.
        For relative and absolute tolerance keywords, ``rtol`` and ``atol``,
        6-element arrays may be used, where the first three elements are for
        position, and the last three are for velocity.

    """

    # For quick reference
    _GM: float = (const.G * const.M_sun).to_value("km3/s2")

    @property
    def GM(self):
        """Gravitational constant times mass."""
        return u.Quantity(self._GM, "km3/s2")

    @classmethod
    def dx_dt(cls, t: float, rv: np.ndarray, beta: float, *args) -> np.ndarray:
        r = rv[:3]
        v = rv[3:]

        r2 = (r**2).sum()
        r1 = np.sqrt(r2)
        r3 = r2 * r1
        GM_r3 = cls._GM / r3 * (1 - beta)

        dx_dt = np.empty(6)
        dx_dt[:3] = v
        dx_dt[3:] = -r * GM_r3

        return dx_dt

    @classmethod
    def df_drv(cls, t: float, rv: np.ndarray, beta: float, *args) -> np.ndarray:
        r = rv[:3]
        r2 = (r**2).sum()
        r1 = np.sqrt(r2)
        r3 = r1 * r2
        GM_r5 = cls._GM * (1 - beta) / (r2 * r3)

        # df_drv[i, j] = df_i/drv_j
        df_drv = np.zeros((6, 6))

        df_drv[0, 3] = 1
        df_drv[1, 4] = 1
        df_drv[2, 5] = 1

        df_drv[3, 0] = GM_r5 * (r2 - 3 * r[0] * r[0])
        df_drv[3, 1] = -GM_r5 * 3 * r[0] * r[1]
        df_drv[3, 2] = -GM_r5 * 3 * r[0] * r[2]

        df_drv[4, 0] = -GM_r5 * 3 * r[1] * r[0]
        df_drv[4, 1] = GM_r5 * (r2 - 3 * r[1] * r[1])
        df_drv[4, 2] = -GM_r5 * 3 * r[1] * r[2]

        df_drv[5, 0] = -GM_r5 * 3 * r[2] * r[0]
        df_drv[5, 1] = -GM_r5 * 3 * r[2] * r[1]
        df_drv[5, 2] = GM_r5 * (r2 - 3 * r[2] * r[2])

        return df_drv
