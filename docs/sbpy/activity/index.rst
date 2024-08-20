Activity Module (`sbpy.activity`)
=================================

`sbpy.activity` models cometary dust and gas activity.  It is separated into two main sub-modules: :doc:`dust` and :doc:`gas`.  The base module itself defines photometric apertures that may be useful for observations of comets.

.. toctree::
  :maxdepth: 2
  :caption: Table of Contents

  dust
  gas


Apertures
---------

Four photometric aperture classes are defined, primarily for use with cometary comae:

  * `~sbpy.activity.CircularAperture`: a circle,
  * `~sbpy.activity.AnnularAperture`: an annulus,
  * `~sbpy.activity.RectangularAperture`: a rectangle,
  * `~sbpy.activity.GaussianAperture`: a Gaussian-weighted beam.

All apertures assume they are centered on the source, i.e., offsets are not supported, but may be added in the future.

Each object is initialized with the dimensions of the aperture, which may be in linear or angular units:

  >>> import astropy.units as u
  >>> import sbpy.activity as sba
  >>> 
  >>> sba.CircularAperture(10 * u.arcsec)
  <CircularAperture: radius 10.0 arcsec>
  >>>
  >>> sba.RectangularAperture((2000, 5000) * u.km)
  <RectangularAperture: dimensions 2000.0×5000.0 km>

Apertures may be converted between linear and angular units using :func:`~sbpy.activity.Aperture.as_angle` and :func:`~sbpy.activity.Aperture.as_length`.  The conversion requires the observer-target distance (``'delta'``) as a `~astropy.units.Quantity` or `~sbpy.data.Ephem`.

  >>> ap = sba.CircularAperture(1 * u.arcsec)
  >>> ap.as_length(1 * u.au)  # doctest: +FLOAT_CMP
  <CircularAperture: radius [725.27094381] km>

Ideal comae (constant production rate, free-expansion, infinite lifetime) have *1/ρ* surface brightness distributions.  With :func:`~sbpy.activity.Aperture.coma_equivalent_radius`, we may convert the aperture into a circular aperture that would contain the same total flux at the telescope:

  >>> ap = sba.RectangularAperture((2000, 5000) * u.km)
  >>> sba.CircularAperture(ap.coma_equivalent_radius())  # doctest: +FLOAT_CMP
  <CircularAperture: radius 1669.4204086589311 km>

Through the ``coma_equivalent_radius()`` method, all apertures may be used to initialize a ``CircularAperture`` instance using the :func:`~sbpy.activity.CircularAperture.from_coma_equivalent` method:

  >>> sba.CircularAperture.from_coma_equivalent(ap)
  <CircularAperture: radius 1669.4204086589311 km>

Reference/API
-------------
.. automodapi:: sbpy.activity
    :no-heading:
    :include: Aperture
    :include: CircularAperture
    :include: AnnularAperture
    :include: RectangularAperture
    :include: GaussianAperture
