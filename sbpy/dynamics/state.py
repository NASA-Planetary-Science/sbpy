# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy.dynamics.state
===================

Dynamical state.


"""

__all__ = [
    "ArbitraryFrame",
    "StateBase",
    "State",
    "SolverFailed",
]

import abc
from typing import Iterable, List, Optional, Tuple, TypeVar, Union
from packaging.version import Version

try:
    # python 3.11 feature
    from typing import Self
except ImportError:
    Self = TypeVar("Self", bound="StateBase")

import numpy as np
import astropy
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import frame_transform_graph, SkyCoord, BaseCoordinateFrame
import astropy.coordinates.representation as cr

from .. import data as sbd
from ..data.ephem import Ephem
from ..exceptions import SbpyException


class SolverFailed(SbpyException):
    """DynamicalModel solver failed."""


FrameInputTypes = TypeVar("FrameInputTypes", str, BaseCoordinateFrame)


class ArbitraryFrame(BaseCoordinateFrame):
    """Coordinate frame with arbitrary space and time references."""

    default_representation = cr.CartesianRepresentation
    default_differential = cr.CartesianDifferential


if Version("6.1.0") <= Version(astropy.__version__) <= Version("6.1.1"):
    ArbitraryFrame.separation.__doc__ = ArbitraryFrame.separation.__doc__.replace(
        ":doc:`astropy:/", ":doc:`astropy:"
    )
    ArbitraryFrame.separation_3d.__doc__ = ArbitraryFrame.separation_3d.__doc__.replace(
        ":doc:`astropy:/", ":doc:`astropy:"
    )


class StateBase(abc.ABC):
    """Abstract base class for dynamical state.


    Parameters
    ----------
    r : `~astropy.units.Quantity`
        Position (x, y, z), shape = (3,) or (N, 3).

    v : `~astropy.units.Quantity`
        Velocity (x, y, z), shape = (3,) or (N, 3).

    t : `~astropy.time.Time` or `~astropy.units.Quantity`
        Time, a scalar or shape = (N,).  Time as a `Quantity` may only be used
        with the `ArbitraryFrame` coordinate frame.

    frame : `~astropy.coordinates.BaseCoordinateFrame` class or string, optional
        Reference frame for ``r`` and ``v``.  Default:
        `~sbpy.dynamics.ArbitraryFrame`.


    Examples
    --------

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from sbpy.dynamics import State
    >>> r = [1e9, 1e9, 0] * u.km
    >>> v = [0, 0, 10] * u.km / u.s
    >>> t = Time("2022-07-24", scale="tdb")
    >>> state = State(r, v, t)

    Or, specify time relative to an arbitrary epoch:

    >>> t = 327 * u.day
    >>> state = State(r, v, t)


    Notes
    -----

    State data is immutable.

    """

    def __init__(
        self,
        r: u.Quantity,
        v: u.Quantity,
        t: Union[u.Quantity, Time],
        frame: Optional[FrameInputTypes] = None,
    ) -> None:
        frame_class: BaseCoordinateFrame = self._get_frame_class(frame)

        # some frames require observation time for coordinate transformations
        frame_kwargs: dict = {}
        if "obstime" in frame_class.frame_attributes:
            frame_kwargs["obstime"] = t

        xyz_axis: int = 1 if np.ndim(r) != 1 else 0
        self._data: BaseCoordinateFrame = frame_class(
            cr.CartesianRepresentation(
                u.Quantity(r),
                xyz_axis=xyz_axis,
                differentials={
                    "s": cr.CartesianDifferential(u.Quantity(v), xyz_axis=xyz_axis)
                },
            ),
            representation_type="cartesian",
            **frame_kwargs,
        )

        if self.r.ndim > 2:
            raise ValueError(
                "`State` only supports 1 and 2 dimensional `r` and `v` vectors."
            )

        if isinstance(t, u.Quantity) and not isinstance(self._data, ArbitraryFrame):
            raise TypeError(
                "`State` only supports time as a quantity with `ArbitraryFrame`."
            )

        t_: Union[Time, u.Quantity, list]
        if np.size(t) != len(self):
            t_ = [t] * len(self)
        else:
            t_ = t

        try:
            self._t = u.Quantity(t_)
        except TypeError:
            self._t = Time(t_)

    def __repr__(self) -> str:
        return (
            f"<{type(self).__name__} ({self.frame}):\n"
            + f" r\n  {self.r}\n"
            + f" v\n  {self.v}\n"
            + f" t\n  {self.t}>"
        )

    def __len__(self):
        """Number of state vectors in this object."""
        if self.r.ndim == 1:
            return 1
        else:
            return self.r.shape[0]

    def __getitem__(self, k: Union[int, tuple, slice]) -> Self:
        """Get the state(s) at ``k``."""
        return State(self.r[k], self.v[k], self.t[k], frame=self.frame)

    def __add__(self, other: Self) -> Self:
        """Vector addition of two states.

        Time is taken from the left operand.

        """

        return State(
            self.r + other.r,
            self.v + other.v,
            self.t,
            frame=self.frame,
        )

    def __sub__(self, other: Self) -> Self:
        """Vector subtraction of two states.

        Time is taken from the left operand.

        """

        return self + -other

    def __neg__(self) -> Self:
        """Invert the direction of the state vector position and velocity."""
        return State(
            -self.r,
            -self.v,
            self.t,
            frame=self.frame,
        )

    def __abs__(self) -> Tuple[u.Quantity, u.Quantity]:
        """Return the magnitude of the position and velocity."""
        r = np.sqrt(np.sum(self.r**2, axis=-1))
        v = np.sqrt(np.sum(self.v**2, axis=-1))
        return r, v

    @staticmethod
    def _get_frame_class(
        frame_input: Union[None, FrameInputTypes]
    ) -> Union[None, BaseCoordinateFrame]:
        """Get a frame class based on allowed ``State`` frame input."""

        frame_class: BaseCoordinateFrame
        if frame_input is None:
            frame_class = ArbitraryFrame
        elif isinstance(frame_input, str):
            frame_class = frame_transform_graph.lookup_name(frame_input)
            if frame_class is None:
                raise ValueError(f"Invalid frame name: {frame_input}")
        elif isinstance(frame_input, BaseCoordinateFrame):
            frame_class = type(frame_input)
        else:
            frame_class = frame_input

        return frame_class

    @property
    def r(self) -> u.Quantity:
        """Position vector."""
        return self._data.cartesian.get_xyz().T

    @property
    def x(self) -> u.Quantity:
        """x component of the position vector."""
        return self.r[..., 0]

    @property
    def y(self) -> u.Quantity:
        """y component of the position vector."""
        return self.r[..., 1]

    @property
    def z(self) -> u.Quantity:
        """z component of the position vector."""
        return self.r[..., 2]

    @property
    def v(self) -> u.Quantity:
        """Velocity vector."""
        return self._data.cartesian.differentials["s"].get_d_xyz().T

    @property
    def v_x(self) -> u.Quantity:
        """x component of the velocity vector."""
        return self.v[..., 0]

    @property
    def v_y(self) -> u.Quantity:
        """y component of the velocity vector."""
        return self.v[..., 1]

    @property
    def v_z(self) -> u.Quantity:
        """z component of the velocity vector."""
        return self.v[..., 2]

    @property
    def rv(self) -> np.ndarray:
        """Position in km, and velocity in km/s."""
        return np.hstack([self.r.to_value("km"), self.v.to_value("km/s")])

    @property
    def t(self) -> Union[Time, u.Quantity]:
        """Time."""
        return self._t

    @property
    def arbitrary_time(self) -> bool:
        """True if the time attribute is arbitrary."""
        return isinstance(self.t, u.Quantity)

    @property
    def frame(self) -> Union[BaseCoordinateFrame, None]:
        return self._data.replicate_without_data()

    def to_skycoord(self) -> SkyCoord:
        """Convert to a `~astropy.coordinates.SkyCoord` object."""

        if hasattr(self._data, "obstime"):
            return SkyCoord(self._data, representation_type="cartesian")
        else:
            return SkyCoord(self._data, obstime=self.t, representation_type="cartesian")

    def transform_to(self, frame: FrameInputTypes) -> Self:
        """Transform state into another reference frame.


        Parameters
        ----------
        frame : string or `~astropy.coordinates.BaseCoordinateFrame`
            Transform into this reference frame.


        Returns
        -------
        state : State
            The transformed state.

        """

        frame_class: BaseCoordinateFrame = self._get_frame_class(frame)
        frame_kwargs: dict = {}

        if "obstime" in frame_class.frame_attributes:
            frame_kwargs["obstime"] = self.t

        transformed: BaseCoordinateFrame = self._data.transform_to(
            frame_class(**frame_kwargs)
        )

        return State(
            transformed.cartesian.get_xyz().T,
            transformed.cartesian.differentials["s"].get_d_xyz().T,
            self.t,
            frame=transformed.replicate_without_data(),
        )


class State(StateBase):
    @classmethod
    def from_states(cls, states: Iterable[Self]) -> Self:
        """Initialize from a list of states.

        The resulting reference frame will be that of ``states[0]``.


        Parameters
        ----------
        states : array

        """

        frame: BaseCoordinateFrame = states[0].frame
        states_: List[State] = [state.transform_to(frame) for state in states]

        r: List[u.Quantity] = [state.r for state in states_]
        v: List[u.Quantity] = [state.v for state in states_]
        t: List[Union[u.Quantity, Time]] = [state.t for state in states_]

        return State(r, v, t, frame=frame)

    @classmethod
    def from_skycoord(cls, coords: SkyCoord) -> Self:
        """Initialize from astropy `~astropy.coordinates.SkyCoord`.


        Parameters
        ----------
        coords: ~astropy.coordinates.SkyCoord
            The object state.  Must have position and velocity, ``obstime``,
            and be convertible to cartesian (3D) coordinates.

        """

        _coords: SkyCoord = coords.copy()
        _coords.representation_type = "cartesian"

        r: u.Quantity = u.Quantity([_coords.x, _coords.y, _coords.z]).T
        v: u.Quantity = u.Quantity([_coords.v_x, _coords.v_y, _coords.v_z]).T
        t: Time = coords.obstime
        return cls(r, v, t, frame=coords.frame.replicate_without_data())

    @classmethod
    @sbd.dataclass_input
    def from_ephem(
        cls,
        eph: Ephem,
        frame: Optional[FrameInputTypes] = None,
    ) -> Self:
        """Initialize from an `~sbpy.data.Ephem` object.


        Parameters
        ----------
        eph : ~sbpy.data.ephem.Ephem
            Ephemeris object, must have time, position, and velocity.  Position
            and velocity may be specified using ("x", "y", "z", "vx", "vy", and
            "vz"), or ("ra", "dec", "Delta", "RA*cos(Dec)_rate", "Dec_rate",
            and "deltadot").

        frame : string or `~astropy.coordinates.BaseCoordinateFrame`, optional
            The reference frame for the ephemeris.

        """

        rectangular: Tuple[str] = ("x", "y", "z", "vx", "vy", "vz", "date")
        spherical: Tuple[str] = (
            "ra",
            "dec",
            "Delta",
            "RA*cos(Dec)_rate",
            "Dec_rate",
            "deltadot",
            "date",
        )

        if all([x in eph for x in rectangular]):
            r: u.Quantity = (
                u.Quantity([eph["x"], eph["y"], eph["z"]]).reshape((3, len(eph))).T
            )
            v: u.Quantity = (
                u.Quantity([eph["vx"], eph["vy"], eph["vz"]]).reshape((3, len(eph))).T
            )
            return cls(r, v, eph["date"], frame=frame)
        elif all([x in eph for x in spherical]):
            c: cr.SphericalRepresentation = cr.SphericalRepresentation(
                eph["ra"], eph["dec"], eph["Delta"]
            )
            d: cr.SphericalDifferential = cr.SphericalDifferential(
                eph["RA*cos(Dec)_rate"],
                eph["Dec_rate"],
                eph["deltadot"],
            ).to_cartesian(base=c)
            c = c.to_cartesian()

            r: u.Quantity = u.Quantity([c.x, c.y, c.z]).reshape((3, len(c))).T
            v: u.Quantity = u.Quantity([d.x, d.y, d.z]).reshape((3, len(c))).T
            return cls(r, v, eph["date"], frame=frame)

        raise ValueError(
            "`Ephem` does not have the required time, position, and/or"
            " velocity fields."
        )

    def observe(self, target: Self) -> SkyCoord:
        """Project a target's position onto the sky.


        Parameters
        ----------
        target : State
            The target to observe.


        Returns
        -------
        coords : SkyCoord


        """

        # transform into the observer's reference frame
        target_in_frame: State = target.transform_to(self.frame)

        coords: SkyCoord = (target_in_frame - self).to_skycoord()
        coords.representation_type = "spherical"
        return coords


State.__init__.__doc__ = """Dynamical state.""" + "\n".join(
    StateBase.__doc__.splitlines()[1:]
)
