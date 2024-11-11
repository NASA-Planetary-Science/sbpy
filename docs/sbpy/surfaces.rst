Surfaces Module (`sbpy.surfaces`)
=================================

The ``surfaces`` module describes the interaction of electromagnetic radiation with surfaces.  Sbpy uses the :math:`(i, e, \phi)` model (angle of incidence, angle of emittance, and phase angle) to describe how light scatters and emits light.  It has a flexible system that can incorporate any surface scattering model that can be described with these three angles.  A few built-in surface models are provided.


.. figure:: ../static/scattering-vectors.svg
    :alt: Diagram of surface scattering and emission vectors

    Sbpy's geometric basis for surface scattering and emission: :math:`n` is the surface normal vector, :math:`r_s` is the radial vector to the light source, and :math:`r_o` is the radial vector to the observer.  The angle of incidence (:math:`i`), angle of emittance (:math:`e`), phase angle (:math:`\phi`) are labeled.

A instance of a ``Surface`` will have methods to calculate absorptance, emittance, and reflectance.  A radiance method is used to calculate the observed spectral radiance of a surface given incident light.

Surfaces are expected to require albedo and/or emissivity.  Conventions on which property is used and when should be defined by each class.  For example, a surface that only calculates reflectance may only require albedo, but one that calculates thermal emission may use the convention of albedo for absorbed sunlight and emissivity for emitted thermal radiation.


Built-in surface models
-----------------------

The model `sbpy.surfaces.scattered.LambertianSurfaceScatteredSunlight` is used to observe sunlight scattered from a Lambertian surface (light scattered uniformly in all directions).

Create an instance of the ``LambertianSurfaceScatteredSunlight`` model, and calculate the absorptance, emittance, and reflectance for :math:`(i, e, \phi) = (30^\circ, 60^\circ, 90^\circ)`.

.. code:: python

    >>> import astropy.units as u
    >>> import matplotlib.pyplot as plt
    >>> from sbpy.surfaces.scattered import LambertianSurfaceScatteredSunlight
    >>>
    >>> surface = LambertianSurfaceScatteredSunlight({"albedo": 0.1})
    >>>
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>> surface.absorptance(i)  # doctest: +FLOAT_CMP
    <Quantity [0.77942286]>
    >>> surface.emittance(e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.5>
    >>> surface.reflectance(i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity [0.04330127]>


Building your own surface models
--------------------------------
