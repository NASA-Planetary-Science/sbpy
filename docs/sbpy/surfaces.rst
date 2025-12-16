Surfaces Module (`sbpy.surfaces`)
=================================

Introduction
------------

.. admonition:: warning

    The surfaces module is being made available on a preview basis.  The API is subject to change.  Feedback on the approach is welcome.

The ``surfaces`` module describes the interaction of electromagnetic radiation with surfaces.  Sbpy uses the :math:`(i, e, \phi)` model (angle of incidence, angle of emittance, and phase angle) to describe how light scatters and emits light.  It has a flexible system that can incorporate any surface scattering model that can be described with these three angles.

.. figure:: ../static/scattering-vectors.svg
    :alt: Diagram of surface scattering and emission vectors

    Sbpy's geometric basis for surface scattering and emission: :math:`n` is the surface normal vector, :math:`r_s` is the radial vector to the light source, and :math:`r_o` is the radial vector to the observer.  The angle of incidence (:math:`i`), angle of emittance (:math:`e`), phase angle (:math:`\phi`) are labeled.

An implementation of the ``Surface`` model has methods to calculate electromagnetic absorption, emission, and bidirectional reflectance.


Getting Started
---------------

Currently the Lambertian surface model is implemented.  A Lambertian surface absorbs and emits light uniformly in all directions.

Create an instance of the ``LambertianSurface`` model, and calculate the absorption and emission scale factors for :math:`(i, e, \phi) = (30^\circ, 60^\circ, 90^\circ)`::

    >>> import astropy.units as u
    >>> from sbpy.surfaces import LambertianSurface
    >>>
    >>> surface = LambertianSurface()
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>> 
    >>> surface.absorption(i)  # doctest: +FLOAT_CMP
    <Quantity 0.8660254>
    >>> surface.emission(e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.5>

Calculate the bidirectional reflectance for :math:`e=0`, and a range of incident angles::

.. plot::
    :context: reset
    :nofigs:

    >>> import numpy as np
    >>>
    >>> e = 0 * u.deg
    >>> i = np.linspace(-90, 90) * u.deg
    >>> phi = np.abs(e - i)  # calculate phase angle
    >>> r = surface.reflectance(i, e, phi)

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
    >>> surface.reflectance_from_vectors(n, r, ro)  # doctest: +FLOAT_CMP
    <Quantity 0.13783222 1 / sr>


Build Your Own Surface Models
-----------------------------

To define your own surface model create a new class based on `~sbpy.surfaces.surface.Surface`, and define the methods for ``absorption``, ``emission``, and ``reflectance``.  The `~sbpy.surfaces.lambertian.LambertianSurface` model may serve as a reference.

Here, we define a new surface model with surface properties proportional to :math:`\cos^2`::

    >>> # use min_zero_cos(a) to ensure cos(a >= 90 deg) = 0
    >>> from sbpy.surfaces.surface import Surface, min_zero_cos
    >>>
    >>> class Cos2Surface(Surface):
    ...     """Absorption and emission proportional to :math:`\\cos^2`."""
    ...
    ...     def absorption(self, i):
    ...         return min_zero_cos(i)**2
    ...
    ...     def emission(self, e, phi):
    ...         return min_zero_cos(e)**2
    ...
    ...     def reflectance(self, i, e, phi):
    ...         return self.absorption(i) * self.emission(e, phi) / np.pi / u.sr

Create and use an instance of our new model::

    >>> surface = Cos2Surface()
    >>> i, e, phi = [30, 60, 90] * u.deg
    >>>
    >>> surface.absorption(i)  # doctest: +FLOAT_CMP
    <Quantity 0.75>
    >>> surface.emission(e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.25>
    >>> surface.reflectance(i, e, phi)  # doctest: +FLOAT_CMP
    <Quantity 0.0596831 1 / sr>


Reference/API
-------------
.. automodapi:: sbpy.surfaces
   :no-heading:
   :inherited-members:
