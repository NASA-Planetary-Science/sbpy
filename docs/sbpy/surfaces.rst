Surfaces Module (`sbpy.surfaces`)
=================================

Introduction
------------

.. admonition:: warning

    The surfaces module is being made available on a preview basis.  The API is subject to change.  Feedback on the approach is welcome.

The ``surfaces`` module describes the interaction of electromagnetic radiation with surfaces.  Sbpy uses the :math:`(i, e, \phi)` model (angle of incidence, angle of emittance, and phase angle) to describe how light scatters and emits light.  It has a flexible system that can incorporate most surface scattering Models that can be described with these three angles.

.. figure:: ../static/scattering-vectors.svg
    :alt: Diagram of surface scattering and emittance vectors

    Sbpy's geometric basis for surface scattering and emittance: :math:`n` is the surface normal vector, :math:`r_s` is the radial vector to the light source, and :math:`r_o` is the radial vector to the observer.  The angle of incidence (:math:`i`), angle of emittance (:math:`e`), phase angle (:math:`\phi`) are labeled.

An implementation of the ``Surface`` model has methods to calculate electromagnetic absorptance, emittance, and bidirectional reflectance.


Getting Started
---------------

Currently the Lambertian surface model is implemented.  A Lambertian surface absorbs and emits light uniformly in all directions.

Create an instance of the ``LambertianSurface`` model, and calculate the absorptance and emittance scale factors for :math:`(i, e, \phi) = (30^\circ, 60^\circ, 90^\circ)`.  Let the albedo be 0.1 (emissivity = 0.9)::

    >>> import astropy.units as u
    >>> from sbpy.surfaces import LambertianSurface
    >>>
    >>> albedo = 0.1
    >>> epsilon = 1 - albedo
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>> 
    >>> surface = LambertianSurface()
    >>> surface.absorptance(epsilon, i)  # doctest: +FLOAT_CMP
    <Quantity 0.77942286>
    >>> surface.emittance(epsilon, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.45>

Calculate the bidirectional reflectance for :math:`e=0`, and a range of incident angles::

.. plot::
    :context: reset
    :nofigs:

    >>> import numpy as np
    >>>
    >>> e = 0 * u.deg
    >>> i = np.linspace(-90, 90) * u.deg
    >>> phi = np.abs(e - i)  # calculate phase angle
    >>> r = surface.reflectance(albedo, i, e, phi)

.. plot::
    :include-source:
    :context:

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(i, r)
    ax.set(xlabel="$i$ (deg)", ylabel="$r$ (1 / sr)")


Using Vectors Instead of Angles
-------------------------------

As an alternative to using :math:`(i, e, \phi)`, results may be calculated using vectors that define the normal direction, radial vector to the light source, and radial vector to the observer::

    >>> # vectors equivalent to (i, e, phi) = (30, 60, 90) deg
    >>> n = [1, 0, 0]
    >>> r = [0.8660254, 0.5, 0] * u.au
    >>> ro = [0.5, -0.8660254, 0] * u.au
    >>>
    >>> surface.reflectance_from_vectors(albedo, n, r, ro)  # doctest: +FLOAT_CMP
    <Quantity 0.0137832 1 / sr>


Build Your Own Surface Models
-----------------------------

To define your own surface model create a new class based on `~sbpy.surfaces.surface.Surface`, and define the methods for ``absorptance``, ``emittance``, and ``reflectance``.  The `~sbpy.surfaces.lambertian.LambertianSurface` model may serve as a reference.

Here, we define a new surface model with surface properties proportional to :math:`\cos^2`::

    >>> # use min_zero_cos(a) to ensure cos(a >= 90 deg) = 0
    >>> from sbpy.surfaces.surface import Surface, min_zero_cos
    >>>
    >>> class Cos2Surface(Surface):
    ...     """Surface properties proportional to :math:`\\cos^2`."""
    ...
    ...     def absorptance(self, epsilon, i):
    ...         return epsilon * min_zero_cos(i)**2
    ...
    ...     def emittance(self, epsilon, e, phi):
    ...         return epsilon * min_zero_cos(e)**2
    ...
    ...     def reflectance(self, albedo, i, e, phi):
    ...         return albedo * min_zero_cos(i)**2 * min_zero_cos(e)**2 / np.pi / u.sr

Create and use an instance of our new model::

    >>> surface = Cos2Surface()
    >>> albedo = 0.1
    >>> epsilon = 1 - albedo
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>>
    >>> surface.absorptance(epsilon, i)  # doctest: +FLOAT_CMP
    <Quantity 0.675>
    >>> surface.emittance(epsilon, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.225>
    >>> surface.reflectance(albedo, i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.00596831 1 / sr>


Reference/API
-------------
.. automodapi:: sbpy.surfaces
   :no-heading:
   :inherited-members:
