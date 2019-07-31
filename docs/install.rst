
Installation
------------

We recommend that you install the latest version of
`Anaconda Python 3.x <https://www.anaconda.com/download/>`__ on your
system before installing `sbpy`. Make sure that Anaconda Python is
your default Python (this will be asked during the installation process).

Requirements
^^^^^^^^^^^^

`sbpy` has the following requirements that will be automatically taken
of with installation using pip:

* Python 3.6 or later
* `numpy <https://www.numpy.org/>`__ 1.13.0 or later
* pytest 3.1 or later
* `astropy <https://www.astropy.org/>`__
* `astroquery <https://astroquery.readthedocs.io/en/latest/>`__ 0.3.9.dev5089 or later: For retrieval of online data, e.g., ephemerides and orbits.
* `scipy <https://scipy.org/>`__: For numerical integrations in `sbpy.activity.gas` and `sbpy.photometry`, among others.
* `synphot <https://github.com/spacetelescope/synphot_refactor>`__ 0.1.3 or later: For calibration to Sun and Vega.
* `ginga <https://ejeschke.github.io/ginga/>`__ and `photutils <https://photutils.readthedocs.io/en/stable/>`__: For use with `sbpy.imageanalysis`

The following packages will have to be installed manually, if the user
wants to use them:

* `oorb <https://github.com/oorb/oorb/tree/master/python>`__: For
  orbit transformations (`~sbpy.data.Orbit.oo_transform`) and
  propagations (`~sbpy.data.Orbit.oo_propagate`), as well as
  ephemerides calculations (`~sbpy.data.Ephem.from_oo`).


Using pip
^^^^^^^^^

The latest development version of `sbpy` can be easily installed using

.. code-block:: bash

    $ pip install git+https://github.com/NASA-Planetary-Science/sbpy.git


Using GitHub
^^^^^^^^^^^^

This way of installing `sbpy` is recommended if you plan to contribute
to the module. The current development version of `sbpy` can be
obtained from `GitHub <https://github.com/NASA-Planetary-Science/sbpy>`__ using

.. code-block:: bash

    $ git clone https://github.com/NASA-Planetary-Science/sbpy.git

This will create a new directory (``sbpy/``). In this directory, run

.. code-block:: bash

    $ python setup.py install --user

in order to use `sbpy` in your default Python environment. If you plan to work on the code and always want to use the latest version of your code, you can install it with


.. code-block:: bash

    $ python setup.py develop --user
