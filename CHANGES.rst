0.5.0 (2024-08-28)
==================

- Revised required and optional packages:

  - Only numpy and astropy are required; scipy, synphot, ads, and astroquery are
    now optional dependencies.  [#383]

  - Created an option to install a recommended list of packages, e.g., ``pip
    install sbpy[recommended]``.  [#383]

  - Require `astropy` >= 5.3.3 to avoid a security vulnerability.  [#400]


New Features
------------

sbpy.activity
^^^^^^^^^^^^^
- New `sbpy.activity.CircularAperture.from_coma_equivalent()` to immediately
  create a `CircularAperture` from any other `Aperture` given a nominal coma
  surface brightness distribution. [#393]

sbpy.dynamics
^^^^^^^^^^^^^

- New class to encapsulate dynamical state vectors: `sbpy.dynamics.state.State`. [#394]

- New support for dynamical integrations: `sbpy.dynamics.models`.  Three models
  are implemented: `FreeExpansion`, `SolarGravity`, and `SolarGravityAndRadiationPressure`.  [#394]

- New ability to produce dust syndynes and synchrones: `sbpy.dynamics.syndynes`.  [#394]

sbpy.time
^^^^^^^^^

- New sub-module to support conversions to NAIF SPICE ephemeris time. [#394]

sbpy.utils
^^^^^^^^^^
- New `required_packages` and `optional_packages` functions to test for the
  presence of required and optional packages.  [#383]

sbpy.utils.decorators
^^^^^^^^^^^^^^^^^^^^^
- New `requires` and  `optionally_uses` function decorators to simplify testing
  for required and optional packages.  [#383]


API Changes
-----------

sbpy.sources
^^^^^^^^^^^^
* Deprecated ``SynphotRequired``.  Use
  ``sbpy.execptions.RequiredPackageUnavailable``.  [#383]


Bug Fixes
---------
* ``sbpy.sources.SpectralSource`` now correctly raises
  ``RequiredPackageUnavailable`` when ``synphot`` is not available, replacing a
  locally defined ``SynphotRequired`` or the generic ``ImportError``. [#383]

* Fix `sbpy.photometry.bandpass` crash when reading in SDSS filter transmission
  curves using `synphot` 1.4. [#400]

sbpy.calib
^^^^^^^^^^
* Fixed: the CALSPEC solar spectrum was missing from `solar_sources`. [#387]


0.4.0 (2023-06-30)
==================

* Updated minimum supported versions:
  - Python 3.8
  - `numpy` 1.18
  - `astropy` 4.3
  - `synphot` 1.1
  - `astroquery` 0.4.5


New Features
------------

sbpy.activity
^^^^^^^^^^^^^

- Added ``VectorialModel.binned_production`` constructor for compatibility with
  time-dependent production implemented in the original FORTRAN vectorial model
  code by Festou. [#336]

- Added ``VMResult``, ``VMFragmentSputterPolar``, ``VMParams``,
  ``VMGridParams``, ``VMFragment``, and ``VMParent`` dataclasses to expose
  details of ``VectorialModel`` results that may be of interest. [#336]

sbpy.calib
^^^^^^^^^^

- Added a model spectrum of the Sun from STScI's CALSPEC database (Bohlin et al.
  2014, PASP 126, 711, DOI:10.1086/677655). [#371]

sbpy.data
^^^^^^^^^

- Added ``Orbit.tisserand`` to calculate the Tisserand parameter of small body's
  orbits with respect to planets. [#325]

- Added ``Orbit.D_criterion`` to evaluate the D-criterion between two sets of
  orbital elements. [#325]

- Added ``DataClass.__contains__`` to enable `in` operator for ``DataClass``
  objects. [#357]
  
- Added ``DataClass.add_row``, ``DataClass.vstack`` methods. [#367]

sbpy.photometry
^^^^^^^^^^^^^^^

- Added parameter constraints to the IAU disk-integrated phase function models,
  such as ``HG``, ``HG1G2``, ``HG12``, and ``HG12_Pen16``. [#366]

Documentation
^^^^^^^^^^^^^

- Index page has been reorganized. [#337]


API Changes
-----------

sbpy.activity
^^^^^^^^^^^^^

- ``VectorialModel`` now no longer takes an ``angular_substeps`` parameter. [#336]

sbpy.data
^^^^^^^^^

- IAU HG series functions moved from `sbpy.photometry.core` to `sbpy.photometry.iau`. [#354]

sbpy.photometry
^^^^^^^^^^^^^^^

- Replaced ``NonmonotonicPhaseFunctionWarning`` with
  ``InvalidPhaseFunctionWarning``. [#366]


Bug Fixes
---------

sbpy.calib
^^^^^^^^^^

- Updated STScI URLs for solar spectra (Castelli and Kurucz models). [#345]

sbpy.data
^^^^^^^^^

- Cometary magnitudes obtained via ``Phys.from_sbdb`` (i.e., M1 and M2) now have
  appropriate units. [#349]

- Asteroids with A/ designations (e.g., A/2019 G2) are correctly identified by
  ``Names`` as asteroids.  Improved handling of interstellar object (I/)
  designations: they do not parse as cometary or asteroidal. [#334, #340]


0.3.0 (2022-04-28)
==================

New Features
------------

sbpy.activity
^^^^^^^^^^^^^

- New ``VectorialModel`` to implement the Festou (1981) model of the same name.
  The code reproduces tests based on the literature within 20%, but the causes
  of the differences are unknown.  Help testing this new feature is appreciated.
  [#278, #305]

sbpy.data
^^^^^^^^^

- ``DataClass`` fields are now checked for physically consistent units (e.g.,
  heliocentric distance in units of length), or that they are ``Time`` objects,
  as appropriate. [#275]

sbpy.photometry
^^^^^^^^^^^^^^^

- Add ATLAS c and o bandpasses to ``bandpass``. [#258]

sbpy.spectroscopy
^^^^^^^^^^^^^^^^^

- Add the ability to redden ``SpectralSource`` (such as the ``Sun`` model in
  ``sbpy.calib``) with a new ``.redden()`` method. [#289]


Bug Fixes
---------

sbpy.activity
^^^^^^^^^^^^^

- Allow apertures to be astropy ``Quantity`` objects in ``GasComa`` models,
  e.g., ``Haser``.  [#306]

sbpy.data
^^^^^^^^^
- Corrected ``Orbit.oo_propagate`` units on angles from degrees to radians.
  [#262]
- Corrected ``Orbit`` fields from openorb to use ``'Tp'`` for perihelion date
  epochs as astropy ``Time`` objects, instead of ``'Tp_jd'``. [#262]
- Corrected ``Name.from_packed`` which could not unpack strings including "j".
  [#271]
- Remove hard-coded URL for JPL Horizons and astroquery's ``Horizons`` objects.
  [#295]
- NaNs no longer crash ``Phys.from_sbdb``. [#297]
- When units are not included in the ``Phys.from_sbdb`` results returned from
  NASA JPL, return unit-less values (and any description of the units, such as
  ``'density_sig'``) to the user. [#297]
- ``Names.parse_comet`` now correctly parses Pan-STARRS if included in a comet
  name string, and corrected the label for fragment names in C/ objects:
  ``'fragm'`` --> ``'fragment'`` . [#279]
- Preserve the order of the user's requested epochs in ``Ephem.from_horizons``.
  [#318]

sbpy.photometry
^^^^^^^^^^^^^^^

- Corrected PS1 filter wavelength units in ``bandpass`` from Å to nm. [#258]
- Fix ``HG1G2`` to respect the units on phase angle ``ph`` or else assume
  radians. [#288]

API Changes
-----------

sbpy.data
^^^^^^^^^

- ``DataClass.field_names`` now returns a list of field names rather than a list
  of internal column names. [#275]

Other Changes and Additions
---------------------------

- Improved compatibility with Python 3.8 [#259]
- Added support for astropy 4.0, drop support for astropy 3. [#260, #322]
- Infrastructure updated to use contemporary astropy project standards. [#284]
- Tests may be run in parallel with pytest, e.g., using ``-n auto``. [#297]


0.2.2 (2020-04-27)
==================

New Features
------------
None


Bug Fixes
---------


sbpy.activity
^^^^^^^^^^^^^

- Fix exception from ``Haser`` when ``CircularAperture`` in linear units is
  used. [#240]


sbpy.data
^^^^^^^^^

- ``DataClass.__getitem__`` now always returns a new object of the same
  class, unless a single field name is provided in which case an
  astropy.Table.Column (no units provided) or astropy.units.Quantity
  (units provided) is returned. [#238]

- Fixed ``Ephem.from_horizons`` to skip adding units to the ``'siderealtime'``
  field if it is missing.  Now, the only required field is ``'epoch'``. [#242]

- ``Ephem.from_horizons`` no longer modifies the ``epochs`` parameter in-place.
  [#247]


sbpy.photometry
^^^^^^^^^^^^^^^

- Fixed ``HG12_Pen16`` calculations, which were using the 2010 G1 and G2
  definitions. [#233]

- Use "Partner" under NASA logo. [#249]


API Changes
-----------
None


Other Changes and Additions
---------------------------

sbpy.activity
^^^^^^^^^^^^^

- Test ``Haser.column_density`` output for angular apertures << lengthscale.
  [#243]


website
-------

- Use HTTPS everywhere. [#244]


0.2.1
=====
This version was not released.


Notes
=====

This changelog tracks changes to sbpy starting with version v0.2.  Recommended
subsection titles: New Features, Bug Fixes, API Changes, and Other Changes and
Additions.  Recommended sub-sub-section titles: sbpy submodules, in alphabetical
order.
