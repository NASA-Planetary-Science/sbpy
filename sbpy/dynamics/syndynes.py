# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy activity.dust.syndynes
===========================

Generate cometary dust syndynes and synchrones.

"""

__all__ = [
    "SourceOrbit",
    "SynGenerator",
    "SynStates",
    "SynCollection",
    "Syndyne",
    "Syndynes",
    "Synchrone",
    "Synchrones",
]

import abc
import time
import logging
from typing import Iterable, Optional, TypeVar

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    pass

import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.time import Time
from astropy.table import vstack
from astropy.coordinates import SkyCoord

from ..data import Ephem
from ..utils.decorators import requires
from ..utils import unmasked
from .models import DynamicalModel, SolarGravityAndRadiationPressure
from .state import StateBase, State

SynGeneratorType = TypeVar("SynGeneratorType", bound="SynGenerator")
SynCollectionType = TypeVar("SynCollectionType", bound="SynCollection")


class SynStates(StateBase, abc.ABC):
    """Abstract base class for particle states that make up a syndyne or synchrone."""

    def __init__(
        self,
        source: State,
        betas: np.ndarray,
        ages: u.Quantity,
        r: u.Quantity,
        v: u.Quantity,
        t: Time,
        initial: State,
        observer: Optional[State] = None,
    ) -> None:
        if u.Quantity(r).ndim != 2 or u.Quantity(v).ndim != 2:
            raise ValueError("Syndyne only supports 2 dimensional r and v vectors.")

        self.source: State = source
        self.initial: State = initial
        self.observer: State | None = observer

        # syndynes will be single beta and array of ages, synchrones will be
        # single age and array of betas
        betas_: np.ndarray
        ages_: np.ndarray
        betas_, ages_ = np.broadcast_arrays(betas, ages)
        self.betas: np.ndarray = betas_
        self.ages: u.Quantity = u.Quantity(ages_, ages.unit)

        super().__init__(r, v, t, frame=initial.frame)

        # generate sky coordinates as needed
        self.coords: SkyCoord | None = (
            None if observer is None else observer.observe(self)
        )

    def to_ephem(self) -> Ephem:
        r"""Convert to an sbpy ephemeris object.


        Returns
        -------
        eph : Ephem


        Notes
        -----

        Source and observer states are stored in the `Ephem.meta` attribute.

        =========================   ====================
        Attribute or quantity       ``Ephem`` field name
        =========================   ====================
        beta(s)                     beta_rad
        age(s)                      age
        t, as `Time`                date
        t, as `Quantity`            t_relative
        :math:`|r|`                 r
        :math:`|v \cdot \hat{r}|`   rdot
        coords                      coords
        coords.ra                   ra
        coords.dec                  dec
        coords.lon                  lon
        coords.lat                  lat
        coords.distance             delta
        coords.radial_velocity      deltadot
        x                           x
        y                           y
        z                           z
        v_x                         vx
        v_y                         vy
        v_z                         vz
        initial.x                   x initial
        initial.y                   y initial
        initial.z                   z initial
        initial.v_x                 vx initial
        initial.v_y                 vy initial
        initial.v_z                 vz initial
        initial.t                   t initial
        =========================   ====================

        """

        data: dict = {}
        data["beta_rad"] = self.betas
        data["age"] = self.ages

        if isinstance(self.t, Time):
            data["date"] = self.t
        else:
            data["t_relative"] = self.t

        data["r"] = abs(self)[0]
        data["rdot"] = np.sum(self.r * self.v, 1) / np.sqrt(np.sum(self.r * self.r, 1))

        if self.observer is not None:
            for k, v in self.coords.representation_component_names.items():
                if v in ("lon", "lat"):
                    data[k] = getattr(self.coords, k)
                elif v == "distance":
                    data["delta"] = getattr(self.coords, k)
            data["deltadot"] = self.coords.radial_velocity
            data["coords"] = self.coords

        data["x"] = self.x
        data["y"] = self.y
        data["z"] = self.z
        data["vx"] = self.v_x
        data["vy"] = self.v_y
        data["vz"] = self.v_z
        data["x initial"] = self.initial.x
        data["y initial"] = self.initial.y
        data["z initial"] = self.initial.z
        data["vx initial"] = self.initial.v_x
        data["vy initial"] = self.initial.v_y
        data["vz initial"] = self.initial.v_z
        data["t initial"] = self.initial.t

        meta: dict = {}
        meta["source"] = {
            "r": self.source.r,
            "v": self.source.v,
            "t": self.source.t,
            "frame": self.source.frame,
        }
        if self.observer is None:
            meta["observer"] = None
        else:
            meta["observer"] = {
                "r": self.observer.r,
                "v": self.observer.v,
                "t": self.observer.t,
                "frame": self.observer.frame,
            }

        return Ephem.from_dict(data, meta=meta)

    @requires("matplotlib")
    def plot(
        self,
        ax: mpl.axes.Axes | None = None,
        *,
        wcs: WCS | None = None,
        unit: u.Unit | str = "arcsec",
        **kwargs,
    ) -> None:
        """Plot the coordinates.

        Requires `self.observer`.


        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`, optional
            Plot to this axis object.

        wcs : `~astropy.wcs.WCS`, optional
            Transform coordinates to the image plane using this world coordinate
            system object.

        unit : `~astropy.unit.Unit` or str, optional
            Plot in these angular units.

        **kwargs : dict
            Keyword arguments are passed along to ``ax.plot``.

        """

        if wcs is None:
            # plot offsets from the source in the tangent plane

            if self.observer is None:
                raise ValueError("observer is not defined")

            coords0 = self.observer.observe(self.source)

            dRA = self.coords.ra - coords0.ra
            dDec = self.coords.dec - coords0.dec
            x = dRA.to(unit) / np.cos(self.coords.dec)
            y = dDec.to(unit)
        else:
            # convert coordinates to plot units with the WCS object (avoid
            # passing a masked object, or wcs will complain)
            x, y = wcs.world_to_pixel(unmasked(self.coords))

        ax.plot(x, y, **kwargs)


class Syndyne(SynStates):
    """Collection of particle states that make up a syndyne.


    Parameters
    ----------
    source : `State`
        The source of the syndyne dust.

    beta : float
        The beta-value of this syndyne.

    ages : ~astropy.units.Quantity
        Array of particle ages, shape = (N,).

    r : `~astropy.units.Quantity`
        Position (x, y, z), shape = (N, 3).  Same coordinate frame as ``source``.

    v : `~astropy.units.Quantity`
        Velocity (x, y, z), shape = (N, 3).  Same coordinate frame as ``source``.

    t : `~astropy.time.Time` or `~astropy.units.Quantity`
        Time of observation.

    initial : `State`
        The initial states for each particle.

    observer : `~sbpy.dynamics.State`, optional
        The observer, used to generate sky coordinates.

    """

    def __init__(
        self,
        source: State,
        beta: float,
        ages: u.Quantity,
        r: u.Quantity,
        v: u.Quantity,
        t: Time,
        initial: State,
        observer: Optional[State] = None,
    ) -> None:
        super().__init__(source, beta, ages, r, v, t, initial, observer=observer)

    @property
    def beta(self) -> float:
        """Syndyne beta value."""
        return self.betas[0]


class Synchrone(SynStates):
    """Collection of particle states that make up a synchrone.


    Parameters
    ----------
    source : State
        The source of the synchrone dust.

    betas : array
        The beta-values of this synchrone, shape = (N,).

    age : ~astropy.units.Quantity
        The particle age of this synchrone.

    r : `~astropy.units.Quantity`
        Position (x, y, z), shape = (N, 3).

    v : `~astropy.units.Quantity`
        Velocity (x, y, z), shape = (N, 3).

    t : `~astropy.time.Time` or `~astropy.units.Quantity`
        Time of observation.

    initial : `State`
        The initial states for each particle.

    observer : `~sbpy.dynamics.State`, optional
        The observer, used to generate sky coordinates.

    """

    def __init__(
        self,
        source: State,
        betas: Iterable[float],
        age: u.Quantity,
        r: u.Quantity,
        v: u.Quantity,
        t: Time,
        initial: State,
        observer: Optional[State] = None,
    ) -> None:
        super().__init__(source, betas, age, r, v, t, initial, observer=observer)

    @property
    def age(self) -> u.Quantity:
        """Synchrone age."""
        return self.ages[0]

    @property
    def epoch(self) -> Time:
        """Epoch of synchrone ejection."""
        return self.source.t - self.ages[0]


class SourceOrbit(SynStates):
    """Collection of states that make up an orbit.

    Parameters
    ----------
    source : State
        The orbit is for this source state.

    dt : ~astropy.units.Quantity
        Each point's time relative to the observation time, shape = (N,).

    r : `~astropy.units.Quantity`
        Position (x, y, z), shape = (N, 3).  Same coordinate frame as ``source``.

    v : `~astropy.units.Quantity`
        Velocity (x, y, z), shape = (N, 3).  Same coordinate frame as ``source``.

    t : `~astropy.time.Time` or `~astropy.units.Quantity`
        Time of observation.

    observer : `~sbpy.dynamics.State`, optional
        The observer, used to generate sky coordinates.

    """

    def __init__(
        self,
        source: State,
        dt: u.Quantity,
        r: u.Quantity,
        v: u.Quantity,
        t: Time,
        observer: Optional[State] = None,
    ) -> None:
        super().__init__(source, 0, -dt, r, v, t, source, observer=observer)

    @property
    def dt(self) -> u.Quantity:
        return -self.ages

    @property
    def epoch(self) -> Time:
        """Epoch of each orbit point."""
        return self.source.t - self.ages[0]


class SynCollection:
    """Immutable collection of syndynes or synchrones.


    Parameters
    ----------

    items : array of `Syndyne` or `Synchone`
        The items.

    """

    data_type = SynStates

    def __init__(self, items: Iterable[SynStates]) -> None:
        self._data = []
        if type(items) == type(self._data):
            self._data[:] = items
        elif isinstance(items, SynCollection):
            self._data[:] = items._data
        else:
            self._data[:] = list(items)

        if any([not isinstance(s, self.data_type) for s in self._data]):
            raise TypeError(
                "All items must be an instance of {}".format(self.data_type)
            )

    def __init_subclass__(cls, /, data_type, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.data_type = data_type

    def __len__(self) -> int:
        """Number of items in the container."""
        return len(self._data)

    def __getitem__(self, k: int | tuple | slice) -> SynStates | SynCollectionType:
        # return a SynStates object
        if isinstance(k, int):
            return self._data[k]

        # return a SynCollection
        return type(self)(self._data[k])

    def to_ephem(self) -> Ephem:
        """Convert to an sbpy ephemeris object.

        Only the metadata from the first item is retained.

        """

        if len(self) == 0:
            return Ephem()

        tables: list[Ephem] = [s.to_ephem().table for s in self]

        return Ephem.from_table(
            vstack(tables, metadata_conflicts="error"),
            meta=tables[0].meta,
        )

    @requires("matplotlib")
    def plot(
        self,
        ax: mpl.axes.Axes | None = None,
        *,
        wcs: WCS | None = None,
        unit: u.Unit | str = "arcsec",
        label_format: str | None = None,
        **kwargs,
    ) -> None:
        """Plot the collection.


        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`, optional
            Plot to this axis object.

        wcs : `~astropy.wcs.WCS`, optional
            Transform coordinates to the image plane using this world coordinate
            system object.

        unit : `~astropy.unit.Unit` or str, optional
            Plot in these angular units.

        label_format : str, optional
            A string defining the label format.  ``{beta}`` or ``{age}`` will be
            replaced as appropriate for the syndyne or synchrone.

        **kwargs : dict
            Keyword arguments are passed along to ``ax.plot``.

        """

        ax = plt.gca() if ax is None else ax

        match (label_format, self):
            case (None, Syndynes()):
                label_format = "$\\beta={beta}$"
            case (None, Synchrones()):
                label_format = "$\\Delta t={age}$"
            case (None, _):
                label_format = ""

        for syn in self:
            label = label_format.format(
                beta=getattr(syn, "beta", None), age=getattr(syn, "age", None)
            )
            syn.plot(ax=ax, wcs=wcs, unit=unit, label=label, **kwargs)


class Syndynes(SynCollection, data_type=Syndyne):
    """Collection of syndynes.


    Parameters
    ----------

    items : array of `Syndyne`
        The items.

    """

    def __repr__(self) -> str:
        beta: np.ndarray = np.array([syndyne.beta for syndyne in self])
        return f"<{type(self).__name__}: betas={str(beta)}>"


class Synchrones(SynCollection, data_type=Synchrone):
    """Collection of synchrones.


    Parameters
    ----------

    items : array of `Synchrone`
        The items.

    """

    def __repr__(self) -> str:
        age: u.Quantity = u.Quantity([synchrone.age for synchrone in self])
        return f"<{type(self).__name__}: ages={str(age)}>"


class SynGenerator:
    """Syndyne / synchrone generator for cometary dust.


    Dust is parameterized with ``beta``, the ratio of the force from solar
    radiation pressure (:math:`F_r`) to that from solar gravity (:math:`F_g`):

    .. math::
        \\beta = \\frac{F_r}{F_g}

    For spherical dust grains, ``beta`` reduces to:

    .. math::
        \\beta = \\frac{0.57 Q_{pr}}{\\rho a}

    where :math:`Q_{pr}` is the radiation pressure efficiency averaged over the
    solar spectrum, :math:`\\rho` is the mass density of the grain (g/cm3), and
    :math:`a` is the grain radius (μm) (Burns et al. 1979).


    Parameters
    ----------
    source : State
        State vector (i.e., position and velocity at time) of the object
        producing dust at the time of the observation.  Must be with respect to
        the central mass (e.g., the Sun).

    betas : ~numpy.ndarray
        Array of beta-parameters to be simulated (dimensionless).

    ages : ~astropy.units.Quantity
        Array of particle ages (time).

    observer : State, optional
        State vector of the observer in the same reference frame as ``source``.

    solver : `~sbpy.activity.dust.dynamics.DynamicalModel`, optional
        Solve the equations of motion with this object.  The default solver is
        `SolarGravityAndRadiationPressure`.

    """

    def __init__(
        self,
        source: State,
        betas: Iterable | u.Quantity,
        ages: u.Quantity,
        observer: Optional[State] = None,
        solver: Optional[DynamicalModel] = None,
    ) -> None:
        if len(source) != 1:
            raise ValueError("Only one source state vector allowed.")

        self.source: State = source
        self.betas: u.Quantity = u.Quantity(betas, "").reshape((-1,))
        self.ages: u.Quantity = u.Quantity(ages).reshape((-1,))
        self.observer: State = observer
        self.solver: DynamicalModel = (
            SolarGravityAndRadiationPressure() if solver is None else solver
        )

        self.initialize_states()
        self.solve()

    def __repr__(self) -> str:
        return f"<SynGenerator:\n betas\n    {self.betas}\n ages\n    {self.ages}>"

    @classmethod
    def at_epochs(
        cls,
        source: State,
        betas: Iterable | u.Quantity,
        epochs: Time,
        **kwargs: dict,
    ) -> SynGeneratorType:
        """An alternative constructor that ejects dust at specific times.


        Parameters
        ----------

        source : State
            State vector (i.e., position and velocity at time) of the object
            producing dust at the time of the observation.  Must be with respect
            to the central mass (e.g., the Sun).

        betas : ~numpy.ndarray
            Array of beta-parameters to be simulated (dimensionless).

        epochs : ~astropy.units.Time
            Specific times to produce dust test particles.  The times will be
            converted to particle ages.

        **kwargs : dict
            Any other `SynGenerator` keyword argument.

        """

        ages: u.Quantity = (source.t - epochs).to("s")
        return cls(source, betas, ages, **kwargs)

    def initialize_states(self) -> None:
        """Generate the initial particle states.

        This method is automatically run on initialization.

        """

        states: list[State] = []
        for age in self.ages:
            t_i: Time = self.source.t - age
            state = self.solver.solve(self.source, t_i, 0)
            states.append(state)

        self.initial_states = State.from_states(states)

        logger: logging.Logger = logging.getLogger()
        logger.info("Initialized %d time steps.", self.ages.size)

    def solve(self) -> None:
        """Generate test particle positions by solving the equations of motion.

        This method is automatically run on initialization.

        """

        logger: logging.Logger = logging.getLogger()

        particles: list[State] = []
        t0: float = time.monotonic()
        for i in range(self.betas.size):
            for j in range(self.ages.size):
                particles.append(
                    self.solver.solve(
                        self.initial_states[j], self.source.t, self.betas[i]
                    )
                )
        t1: float = time.monotonic()
        self.particles = State.from_states(particles)

        logger.info(
            "Solved for %d syndynes, %d time steps each.",
            self.betas.size,
            self.ages.size,
        )
        logger.info(f"{(t1 - t0) / self.betas.size / self.ages.size} s/particle.")

    def syndyne(self, i: int) -> Syndyne:
        """Get a single syndyne.


        Parameters
        ----------

        i : int
            Index of the syndyne (same index as the `betas` array).


        Returns
        -------
        syndyne : Syndyne

        """

        n: int = self.ages.size
        indices: slice = slice(i * n, (i + 1) * n)
        state: State = self.particles[indices]

        return Syndyne(
            self.source,
            self.betas[i],
            self.ages,
            state.r,
            state.v,
            state.t,
            self.initial_states[indices],
            observer=self.observer,
        )

    def synchrone(self, i: int) -> Synchrone:
        """Get a single synchrone.


        Parameters
        ----------
        i : int
            Index of the synchrone (same index as the `ages` array).


        Returns
        -------
        synchrone : Synchrone

        """

        n: int = self.ages.size
        indices: slice = slice(i, None, n)
        state: State = self.particles[indices]

        return Synchrone(
            self.source,
            self.betas,
            self.ages[i],
            state.r,
            state.v,
            state.t,
            self.initial_states[indices],
            observer=self.observer,
        )

    def syndynes(self) -> Syndynes:
        """Get a collection of all syndynes.


        Returns
        -------
        syndynes : Syndynes

        """

        return Syndynes([self.syndyne(i) for i in range(len(self.betas))])

    def synchrones(self) -> Synchrones:
        """Get a collection of all synchrones.


        Returns
        -------
        synchrones : Synchrones

        """

        return Synchrones([self.synchrone(i) for i in range(len(self.ages))])

    def source_orbit(self, dt: u.Quantity) -> State | tuple[State, SkyCoord]:
        """Calculate and observe the orbit of the dust source.


        Parameters
        ----------

        dt : `astropy.units.Quantity`
            The times at which to calculate the orbit, relative to the
            observation time.


        Returns
        -------

        orbit : SourceOrbit
            The orbital states.

        """

        states: list[State] = []
        for i in range(len(dt)):
            t: Time = self.source.t + dt[i]
            states.append(self.solver.solve(self.source, t, 0))
        states: State = State.from_states(states)

        return SourceOrbit(self.source, dt, states.r, states.v, states.t, self.observer)
