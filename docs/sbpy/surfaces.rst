Surfaces Module (`sbpy.surfaces`)
=================================

The ``surfaces`` module describes the interaction of electromagnetic radiation with surfaces.  Sbpy uses the :math:`(i, e, \phi)` model (angle of incidence, angle of emittance, and phase angle) to describe how absorptance, emittance, and reflectance vary with incoming and outgoing radiation.  It has a flexible system that can incorporate any surface scattering model that can be described with these three angles.  However, most users will use the built-in surface models.


.. figure:: ../static/scattering-vectors.svg
    :alt: Diagram of surface scattering and emission vectors

    Sbpy's geometric basis for surface scattering and emission: :math:`n` is the surface normal vector, :math:`r_s` is the radial vector to the light source, and :math:`r_o` is the radial vector to the observer.  The angle of incidence (:math:`i`), angle of emittance (:math:`e`), phase angle (:math:`\phi`) are labeled.

A ``Surface`` as methods to 


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

