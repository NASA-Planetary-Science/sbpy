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