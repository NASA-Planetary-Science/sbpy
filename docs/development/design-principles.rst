.. _design principles:

Design Principles
=================

All `sbpy` code must be written according to the following design principles.


Physical parameters are quantities
----------------------------------

* If a variable has units of length, time, flux density, magnitude, etc., then it must be an `astropy` `~astropy.units.Quantity`.

* Inputs may be type and unit checked with the `~astropy.units.quantity_input` decorator.

* Note magnitudes may also carry physical units.  Compare `astropy.units.mag` (unitless) to `astropy.units.ABmag` (flux density per unit frequency), and `sbpy.units.VEGAmag` (units of Vega spectral flux density).

Epochs must be Time objects
---------------------------

* Any kind of epoch or point in time must be of type `~astropy.time.Time`; time scales must be properly set and propagated through all functions.


Use sbpy ``DataClass`` objects
------------------------------

* `~sbpy.data.Orbit`, `~sbpy.data.Phys`, `~sbpy.data.Ephem`, and `~sbpy.data.Obs` are the `sbpy` `~sbpy.data.DataClass` objects.  See :doc:`../sbpy/data/index` for details.

* All inputs based on ephemeris, orbit, and physical parameters must use these classes.

* The classes enable easy parameter passing from online sources.  Compare the following:

    .. code-block:: python

        eph = Ephem.from_horizons('2P')
        # rh, delta required, phase angle is optional:
        Afrho(wave, fluxd, aper, eph['rh'], eph['delta'], phase=eph['phase'])
        # more to the point:
        Afrho(wave, fluxd, aper, eph)

    Carefully document which fields are used by your function or method.
     
* Dictionary-like objects may be allowed for user input, but should be internally converted to a ``DataClass`` object with the `~sbpy.data.dataclass_input` decorator:

    .. code-block:: python

        @dataclass_input(eph=Ephem)
        def H11(eph):
            ...

    The same, but using function annotations:

    .. code-block:: python

        @dataclass_input
        def H11(eph: Ephem):
            ...

* Exceptions are allowed when only one parameter is needed, e.g., ``phase_func(phase)``.  But instead consider using the relevant ``DataClass`` object, and decorating the function with `~sbpy.data.quantity_to_dataclass`:

    .. code-block:: python

        @quantity_to_dataclass(eph=(Ephem, 'phase'))
        def phase_func(eph):
            ...

    The decorator may be stacked with ``dataclass_input`` for maximum
    flexibility:

    .. code-block:: python

        @dataclass_input
        @quantity_to_dataclass(eph=(Ephem, 'phase'))
        def phase_func(eph):
            ...


Append fields to ``DataClass`` at the user's request
----------------------------------------------------

* If a function takes a `~sbpy.data.DataClass` object as input, the
  results of this function may be optionally appended to this object,
  if appropriate.

* Use the keyword argument ``append_results``.

* If ``True``, add the data to the ``DataClass`` object as new fields, and return the result.


Cite relevant works
-------------------

* All important references for methods, data sources, parameters, software, etc., must be cited.

* Citations may be executed internally with :func:`sbpy.bib.register`, or via the `~sbpy.bib.cite` decorator:

    .. code-block:: python

        @cite({'method': '1687pnpm.book.....N'})
        def force(mass, acceleration):
            return mass * acceleration

* Labels describing references (``'method'`` in the above example) are
  required to start with the following strings: ``'method'`` (for
  general method descriptions), ``'data source'`` (for data provenance
  referencing), ``'parameter'`` (for use of specific parameters),
  ``'software'`` (for use of specific software packages). Each of
  these labels can be suffixed at the users discretion (e.g.,
  ```parameter: beaming parameter'`` defines a reference for the value
  used for the beaming parameter).

Exceptions for private functions or speed
------------------------------------------

* ``Quantity`` and ``DataClass`` objects are not required for private methods or functions requiring high performance.

* If a high-performance method is needed, consider writing two methods: one that uses the ``Quantity`` and/or ``DataClass`` objects, and a second that is unitless.

* To simplify code maintenance and testing, the ``Quantity``-loaded method should call the unitless method.
