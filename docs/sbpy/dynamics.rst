**************************
Dynamics (`sbpy.dynamics`)
**************************

`sbpy` has the capability to describe dynamical states (position and velocity vectors), and to integrate test particle orbits.  These are primarily intended to support the generation of :ref:`syndynes`.


Dynamical state
===============

`~sbpy.dynamics.state.State` objects encapsulate the position and velocity of an object at a given time.  Create a `State` for a comet at :math:`x=2` au, moving along the y-axis at a speed of 30 km/s:


.. doctest::

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from sbpy.dynamics import State
    >>> 
    >>> r = [2, 0, 0] * u.au
    >>> v = [0, 30, 0] * u.km / u.s
    >>> t = Time("2023-12-08")
    >>> comet = State(r, v, t)

`State` objects may also represent an array of objects:

.. doctest::

    >>> r = ([2, 0, 0], [0, 2, 0]) * u.au
    >>> v = ([0, 30, 0], [30, 0, 0]) * u.km / u.s
    >>> t = Time(["2023-12-08", "2023-12-09"])
    >>> comets = State(r, v, t)
    >>> len(comets)
    2

The `r`, `v`, and `t` attributes hold the position, velocity, and time for the object(s).  The first index iterates over the object, the second iterates over the x-, y-, and z-axes.

.. doctest::

    >>> comets.r.shape
    (2, 3)
    >>> comets.r[0]
    <Quantity [2., 0., 0.] AU>

Or, index the ``State`` object itself:

.. doctest::

    >>> comets[0].r  # equivalent to comets.r[0]
    <Quantity [2., 0., 0.] AU>


Convert to/from `Ephem` and `SkyCoord`
------------------------------------------

State objects may be initialized from ephemeris objects (`~sbpy.data.Ephem`), provided they contain time, and 3D position and velocity:

.. doctest-requires:: astroquery

    >>> from sbpy.data import Ephem
    >>> eph = Ephem.from_horizons("9P",
    ...                           epochs=Time("2005-07-04"),
    ...                           id_type="designation",
    ...                           closest_apparition=True)  # doctest: +REMOTE_DATA
    >>> tempel1 = State.from_ephem(eph)                     # doctest: +REMOTE_DATA

And `State` may be converted to an `Ephem` object:

.. doctest-requires:: astroquery

    >>> eph = tempel1.to_ephem()  # doctest: +REMOTE_DATA

`astropy`'s `~astropy.coordinates.SkyCoord` objects may also be used, assuming the time and 3D vectors are fully defined:

    >>> from astropy.coordinates import SkyCoord
    >>> coords = SkyCoord(ra="01:02:03h",
    ...                   dec="+04:05:06d",
    ...                   pm_ra_cosdec=100 * u.arcsec / u.hr,
    ...                   pm_dec=-10 * u.arcsec / u.hr,
    ...                   distance=1 * u.au,
    ...                   radial_velocity=15 * u.km / u.s,
    ...                   obstime=Time("2024-01-19"))
    >>> state = State.from_skycoord(coords)

And back to a `SkyCoord` object:

.. doctest::

    >>> state.to_skycoord()
    <SkyCoord (ICRS): (x, y, z) in AU
        (0.96112414, 0.26676914, 0.07123631)
     (v_x, v_y, v_z) in km / s
        (9.16701924, 23.45244422, -0.94097856)>


Reference frames
----------------

Coordinate reference frames can be specified with the `frame` keyword argument during initialization.  Most of `astropy` reference frames are supported (see `astropy's Built-in Frame Classes <https://docs.astropy.org/en/stable/coordinates/index.html#module-astropy.coordinates.builtin_frames>`_):

.. note::
    When working with heliocentric ecliptic coordinates from JPL Horizons or NAIF SPICE, you may want to use the `~astropy.coordinates.HeliocentricEclipticIAU76` reference frame.

.. doctest::

    >>> r = [2, 0, 0] * u.au
    >>> v = [0, 30, 0] * u.km / u.s
    >>> t = Time("2023-12-08")
    >>> state = State(r, v, t, frame="heliocentriceclipticiau76")
    >>> state
    <State (<HeliocentricEclipticIAU76 Frame (obstime=2023-12-08 00:00:00.000)>):
     r
      [2. 0. 0.] AU
     v
      [ 0. 30.  0.] km / s
     t
      2023-12-08 00:00:00.000>

Use :func:`~sbpy.dynamics.state.State.transform_to` to transform between reference frames:

.. doctest::

    >>> new_state = state.transform_to("icrs")
    >>> new_state
     <State (<ICRS Frame>):
      r
       [ 1.99191790e+00 -2.59236051e-03 -8.93633307e-04] AU
      v
       [8.11760173e-03 2.75129298e+01 1.19282350e+01] km / s
      t
       2023-12-08 00:00:00.000>


State lengths, subtraction, and observations
--------------------------------------------

A few mathematical operations are possible.  Get the magnitude of the heliocentric distance and velocity with `abs`:

.. doctest::

    >>> print(abs(comet))
    (<Quantity 2. AU>, <Quantity 30. km / s>)

Get the Earth-comet state vector by subtraction:

.. doctest::

    >>> earth = State([0, 1, 0] * u.km, [30, 0, 0] * u.km / u.s, comet.t)
    >>> print(comet - earth)
    <State (<ArbitraryFrame Frame>):
     r
      [ 2.00000000e+00 -6.68458712e-09  0.00000000e+00] AU
     v
      [-30.  30.   0.] km / s
     t
      2023-12-08 00:00:00.000>

Or, get the sky coordinates of the comet as seen from the Earth:

.. doctest::

    >>> earth.observe(comet)
    <SkyCoord (ArbitraryFrame): (lon, lat, distance) in (deg, deg, AU)
        (359.99999981, 0., 2.)
     (pm_lon, pm_lat, radial_velocity) in (mas / yr, mas / yr, km / s)
        (6.52671946e+08, 0., -30.0000001)>

The result, a `~astropy.coordinates.SkyCoord` object, is expressed in the reference frame of the observer.


Dynamical integrators
=====================

A state object may be propagated to a new time using a dynamical integrator.  Three integrators are defined.  Use `~sbpy.dynamics.FreeExpansion` for motion in free space, `~sbpy.dynamics.SolarGravity` for orbits around the Sun, and `~sbpy.dynamics.SolarGravityAndRadiationPressure`  for orbits around the Sun considering radiation pressure.

.. doctest-requires:: scipy

    >>> from sbpy.dynamics import SolarGravity
    >>> state = State([1, 0, 0] * u.au, [0, 30, 0] * u.km / u.s, 0 * u.s)
    >>> integrator = SolarGravity()
    >>> t_final = 1 * u.year
    >>> integrator.solve(state, t_final)
    <State (<ArbitraryFrame Frame>):
     r
      [ 1.48146925e+08 -2.09358771e+07  0.00000000e+00] km
     v
      [ 4.137801   29.70907176  0.        ] km / s
     t
      1.0 yr>

Other integrators may be defined.  Use the above classes as templates or as base classes.  For example, to simulate orbits around the Earth, update the `_GM` private attribute from `SolarGravity`:

.. doctest-requires:: scipy

  >>> class EarthGravity(SolarGravity):
  ...     _GM = 3.9860043543609598e5  # km3/s2


.. _syndynes:

Dust syndynes and synchrones
============================

Syndynes are lines in space connecting particles that are experiencing the same forces.  A syndyne is parameterized by :math:`\beta`, the ratio of the force from solar radiation to the force from solar gravity, :math:`F_r / F_g`, and age (or time of release).  Thus, all particles in a syndyne have a constant :math:`\beta` but variable age.  Similarly, synchrones are lines of constant particle age, but variable :math:`\beta`.


Syndynes
--------

Zero-ejection velocity syndynes are generated with the `~sbpy.dynamics.syndynes.SynGenerator` class.  The class requires a dust source described by a `~sbpy.dynamics.state.State` object, :math:`\beta` values, and particle ages from which to generate the syndynes.

First, define the source of the syndynes, a comet at 2 au from the Sun:

.. doctest::

   >>> from astropy.time import Time
   >>> import astropy.units as u
   >>> from sbpy.dynamics import State
   >>> 
   >>> r = [2, 0, 0] * u.au
   >>> v = [0, 30, 0] * u.km / u.s
   >>> t = Time("2023-12-08")
   >>> comet = State(r, v, t)

Next, initialize the syndyne object:

.. doctest-requires:: scipy

   >>> import numpy as np
   >>> from sbpy.dynamics import SynGenerator
   >>> 
   >>> betas = [1, 0.1, 0.01, 0]
   >>> ages = np.linspace(0, 100, 26) * u.day
   >>> dust = SynGenerator(comet, betas, ages)
   >>> dust
   <SynGenerator:
   betas
      [1.   0.1  0.01 0.  ]
   ages
      [  0.   4.   8.  12.  16.  20.  24.  28.  32.  36.  40.  44.  48.  52.
    56.  60.  64.  68.  72.  76.  80.  84.  88.  92.  96. 100.] d>

The computed particle positions are saved in the :attr:`~sbpy.dynamics.syndynes.SynGenerator.particles` attribute.  For our example, the 4 :math:`\beta`-values and the 26 ages produce 104 particles:

.. doctest-requires:: scipy

   >>> print(len(dust.particles))
   104

Get the results with the :func:`~sbpy.dynamics.syndynes.SynGenerator.syndynes` method, which returns a list-like collection of `~sbpy.dynamics.syndynes.Syndyne` objects.  `Syndyne` objects are specialized `State` objects.  For example, we can compute the linear distance from the comet to farthest particle in each syndyne:

.. doctest-requires:: scipy

   >>> for syndyne in dust.syndynes():
   ...     r, v = abs(syndyne - comet)
   ...     print(syndyne.beta, "{:.3f}".format(r.max().to("au")), sep=", ")
   1.0, 0.309 AU
   0.1, 0.032 AU
   0.01, 0.003 AU
   0.0, 0.000 AU

Individual syndynes may be retrieved with the :func:`~sbpy.dynamics.syndynes.SynGenerator.syndyne` method and a syndyne index.  The index for the syndyne matches the index of the `betas` array, i.e., to get the :math:`\beta=0.1` syndyne from our example:

.. doctest-requires:: scipy

   >>> print(dust.betas)
   [1.   0.1  0.01 0.  ]
   >>> syndyne = dust.syndyne(1)
   >>> print(syndyne.beta)
   0.1


Synchrones
----------

Synchrones are also simulated with the `SynGenerator` class, but instead retrieved with the :func:`~sbpy.dynamics.syndynes.SynGenerator.synchrone` and :func:`~sbpy.dynamics.syndynes.SynGenerator.synchrones` methods, which return `~sbpy.dynamics.syndynes.Synchrone` objects, or collections thereof:

.. doctest-requires:: scipy

   >>> synchrone = dust.synchrone(24)
   >>> r, v = abs(synchrone)
   >>> print(synchrone.age, "{:.3g}".format(r.max().to("au")), sep=", ")
   96.0 d, 2.26 AU

To support generating synchrones at specific times, the :func:`~sbpy.dynamics.syndynes.SynGenerator.at_epochs` method may be used to initialize a `SynGenerator`:

.. doctest-requires:: scipy

   >>> dates = Time(["2023-10-01", "2023-11-01", "2023-12-01"])
   >>> dust = SynGenerator.at_epochs(comet, betas, dates)
   >>> synchrone = dust.synchrone(0)
   >>> synchrone.epoch
   <Time object: scale='utc' format='iso' value=2023-10-01 00:00:00.000>


Adding ejection velocity
------------------------

An ejection velocity can be added to the dust particles.  This is not implemented in `sbpy`, but is possible by sub-classing `SynGenerator` and overriding the `~sbpy.dynamics.syndynes.SynGenerator.initialize_states` method.  In the following example we add 10 m/s to the z-component of the velocity:

.. doctest-requires:: scipy

   >>> class AlternativeSynGenerator(SynGenerator):
   ...     def initialize_states(self):
   ...         super().initialize_states()
   ...         delta_v = State([0, 0, 0] * u.km, [0, 0, 0.01] * u.km / u.s, 0)
   ...         self.initial_states = self.initial_states + delta_v
   >>>
   >>> dust_with_delta_v = AlternativeSynGenerator.at_epochs(comet, betas, dates)
   >>> dust_with_delta_v.initial_states - dust.initial_states
    <State (<ArbitraryFrame Frame>):
    r
     [[0. 0. 0.]
    [0. 0. 0.]
    [0. 0. 0.]] km
    v
     [[0.   0.   0.01]
    [0.   0.   0.01]
    [0.   0.   0.01]] km / s
    t
     ['2023-10-01 00:00:00.000' '2023-11-01 00:00:00.000'
    '2023-12-01 00:00:00.000']>


Projecting onto the sky
-----------------------

Syndynes and synchrones may be projected onto the sky as seen by a telescope.  This requires an observer.  For precision work, the `SynGenerator` source object's `State` should be defined in a heliocentric reference frame.  Typically, the observer will observe in an equatorial reference frame, but your needs may vary.  Here, we generate syndynes in the `~astropy.coordinates.HeliocentricEclipticIAU76` coordinate frame, which is the same J2000 heliocentric ecliptic coordinate frame that `JPL Horizons <https://ssd.jpl.nasa.gov/horizons/manual.html#frames>`_ and the `NAIF SPICE toolkit <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#Frames%20Supported%20in%20SPICE>`_ use.  We then observe the results in the `~astropy.coordinates.ICRS` reference frame:

.. doctest-requires:: scipy

   >>> r = [2, 0, 0] * u.au
   >>> v = [0, 30, 0] * u.km / u.s
   >>> t = Time("2023-12-08")
   >>> comet = State(r, v, t, frame="heliocentriceclipticiau76")
   >>> observer = State(
   ...     r=[0, 2, 2] * u.au,
   ...     v=[0, 0, 0] * u.km / u.s,
   ...     t=comet.t,
   ...     frame="icrs",
   ... )
   >>> dust = SynGenerator(comet, betas, ages, observer=observer)

With the observer and coordinate frames defined, the generated syndynes and synchrones will have a `coords` attribute, which is an `astropy.coordinates.SkyCoord` object representing the sky positions of the test particles.  Here, we print a sample of the coordinates:

.. doctest-requires:: scipy

   >>> syndyne = dust.syndyne(0)
   >>> print("\n".join(syndyne.coords[::5].to_string("hmsdms", precision=0)))
   20h59m23s -35d18m48s
   21h00m08s -35d12m50s
   21h01m53s -34d55m44s
   21h03m49s -34d30m22s
   21h05m19s -34d00m41s
   21h05m59s -33d30m26s


Source object's projected orbit
-------------------------------

Calculating the positions of the projected orbit of the source object may be helpful for interpreting an observation or a set of syndynes.  They are calculated with the :func:`~sbpy.dynamics.syndynes.SynGenerator.source_orbit` method:

.. doctest-requires:: scipy

   >>> dt = np.linspace(-2, 2) * u.d
   >>> orbit, coords = dust.source_orbit(dt)


Using other dynamical models
----------------------------

`sbpy`'s built-in models solve the equations of motion for dust grains given two-body dynamics.  Users may provide their own models in order to, e.g., improve code performance, add planetary perturbations, model grain fragmentation, etc..  Use the `~sbpy.dynamics.models.SolarGravityAndRadiationPressure` class as a template.

In this example, we compute the syndynes of a comet orbiting β Pic (1.8 solar masses) by sub-classing `SolarGravityAndRadiationPressure` and updating :math:`GM`, the mass of the star times the gravitational constant:

.. doctest-requires:: scipy

   >>> import astropy.constants as const
   >>> from sbpy.dynamics import SolarGravityAndRadiationPressure
   >>>
   >>> class BetaPicGravityAndRadiationPressure(SolarGravityAndRadiationPressure):
   ...     _GM = 1.8 * SolarGravityAndRadiationPressure._GM
   >>>
   >>> solver = BetaPicGravityAndRadiationPressure()
   >>> betapic_dust = SynGenerator(comet, [1, 0], [0, 100] * u.d, solver=solver)


Plotting syndynes and synchrones
--------------------------------

Generally, we are interested in plotting syndynes and synchrones on an image of a comet.  The accuracy of the coordinates object depends on the the comet and observer states, but also on whether or not light travel time is accounted for, and the accuracy of the orbit integrator.  The `sbpy` testing suite shows that arcsecond-level accuracy is possible, but this is generally not accurate enough for direct comparison to typical images of comets.  Instead, it helps to compute the positions of the syndynes and synchrone coordinate objects relative to the comet, and plot the results.

.. doctest-requires:: scipy,matplotlib

   >>> import matplotlib.pyplot as plt
   >>> 
   >>> coords0 = observer.observe(comet)
   >>> def plot(ax, coords, **kwargs):
   ...     dRA = coords.ra - coords0.ra
   ...     dDec = coords.dec - coords0.dec
   ...     ax.plot(dRA.arcsec, dDec.arcsec, **kwargs)
   >>> 
   >>> fig, ax = plt.subplots()
   >>> 
   >>> # plot all but the last (beta=0) syndyne
   >>> for syndyne in dust.syndynes()[:-1]:
   ...     plot(ax, syndyne.coords, label=f"$\\beta={syndyne.beta:.2g}$")
   >>> 
   >>> # plot every 5th synchrone
   >>> for synchrone in dust.synchrones()[4::5]:
   ...     plot(
   ...         ax,
   ...         synchrone.coords,
   ...         ls="--",
   ...         label=f"$\\Delta t={synchrone.age.to(u.d):.2g}$"
   ...     )
   >>> 
   >>> # and plot the orbit
   >>> dt = np.linspace(-2, 2) * u.d
   >>> states, coords = dust.source_orbit(dt)
   >>> plot(ax, coords, color="k", ls=":", label="Orbit")
   >>>
   >>> ax.invert_xaxis()
   >>> plt.setp(ax,
   ...          xlim=[100, -10],
   ...          ylim=[-10, 100],
   ...          xlabel="$\\Delta$RA (arcsec)",
   ...          ylabel="$\\Delta$Dec (arcsec)",
   ... ) # doctest: +SKIP
   >>> plt.legend() # doctest: +SKIP
   >>> plt.tight_layout()

.. plot::

   import numpy as np
   import matplotlib.pyplot as plt

   import astropy.units as u
   from astropy.time import Time
   from sbpy.dynamics import State, SynGenerator

   r = [2, 0, 0] * u.au
   v = [0, 30, 0] * u.km / u.s
   t = Time("2023-12-08")
   frame = "heliocentriceclipticiau76"
   comet = State(r, v, t, frame=frame)

   betas = [1, 0.1, 0.01, 0]
   ages = np.linspace(0, 100, 25) * u.day
   observer = State(
       r=[0, 2, 2] * u.au,
       v=[0, 0, 0] * u.km / u.s,
       t=comet.t,
       frame="icrs",
   )
   dust = SynGenerator(comet, betas, ages, observer=observer)

   coords0 = observer.observe(comet)
   def plot(ax, coords, **kwargs):
       dRA = coords.ra - coords0.ra
       dDec = coords.dec - coords0.dec
       ax.plot(dRA.arcsec, dDec.arcsec, **kwargs)
   
   fig, ax = plt.subplots()
   
   for syndyne in dust.syndynes():
       # don't draw the beta = 0 syndyne
       if syndyne.beta == 0:
           continue
       plot(ax, syndyne.coords, label=f"$\\beta={syndyne.beta:.2g}$")
   
   # plot every 5th synchrone
   for synchrone in dust.synchrones()[4::5]:
       plot(
           ax,
           synchrone.coords,
           ls="--",
           label=f"$\Delta t={synchrone.age.to(u.d):.2g}$",
       )
   
   # and plot the orbit
   dt = np.linspace(-2, 2) * u.d
   states, coords = dust.source_orbit(dt)
   plot(ax, coords, color="k", ls=":", label="Orbit")

   ax.invert_xaxis()
   plt.setp(ax,
            xlim=[100, -10],
            ylim=[-10, 100],
            xlabel="$\Delta$RA (arcsec)",
            ylabel="$\Delta$Dec (arcsec)",
   )
   plt.legend()
   plt.tight_layout()


Reference/API
=============

.. automodapi:: sbpy.dynamics.state
  :no-main-docstr:
  :inherited-members:

.. automodapi:: sbpy.dynamics.models
  :no-main-docstr:
  :inherited-members:

.. automodapi:: sbpy.dynamics.syndynes
  :no-main-docstr:
  :inherited-members:
