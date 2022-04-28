
Installation
------------

Requirements
^^^^^^^^^^^^

`sbpy` has the following requirements that will be automatically taken
care of with installation using pip:

* Python 3.7 or later
* `numpy <https://numpy.org/>`__ 1.17.0 or later
* `astropy <https://www.astropy.org/>`__ 4.0 or later
* `astroquery <https://astroquery.readthedocs.io/en/latest/>`__ 0.4.5 or later: For retrieval of online data, e.g., ephemerides and orbits.
* `scipy <https://www.scipy.org/>`__: For numerical integrations in `sbpy.activity.gas` and `sbpy.photometry`, among others.
* `synphot <https://github.com/spacetelescope/synphot_refactor>`__ 1.0.0 or later: For calibration with respect to the Sun and Vega, filtering spectra through bandpasses.

Optional dependencies
^^^^^^^^^^^^^^^^^^^^^

* Python extensions for `oorb <https://github.com/oorb/oorb/>`__: For orbit
  transformations (`~sbpy.data.Orbit.oo_transform`) and propagations
  (`~sbpy.data.Orbit.oo_propagate`), as well as ephemerides calculations
  (`~sbpy.data.Ephem.from_oo`).
* `pyradex <https://github.com/keflavich/pyradex>`__: For non-LTE production
  rate calculation related to cometary activity (`~sbpy.activity.gas.NonLTE`).
* `ginga <https://ejeschke.github.io/ginga/>`__ and `photutils
  <https://photutils.readthedocs.io/en/stable/>`__: To interactively enhance
  images of comets with the `~sbpy.imageanalysis.CometaryEnhancement` Ginga
  plugin.


Using pip
^^^^^^^^^

The latest stable version of `sbpy` can be installed with:

.. code-block:: bash

    $ pip install sbpy

Most optional dependencies may be installed via:

.. code-block:: bash

    $ pip install sbpy[all]

`oorb` and `pyradex` are left for the user to install manually.

The latest development version of `sbpy` can be easily installed using:

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
