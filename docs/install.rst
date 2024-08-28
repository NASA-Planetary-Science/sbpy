
Installation
------------

Requirements
^^^^^^^^^^^^

`sbpy` has the following requirements that will be automatically taken
care of with installation using pip:

* Python 3.8 or later
* `astropy <https://www.astropy.org/>`__ 5.3.3 or later.
* `numpy <https://numpy.org/>`__ 1.21 or later.

Optional dependencies
^^^^^^^^^^^^^^^^^^^^^

* `ads <https://github.com/andycasey/ads/>`__ 0.12 or later, to fetch citation details for bibliography tracking.  **Recommended**
* `astroquery <https://astroquery.readthedocs.io/en/latest/>`__ 0.4.5 or later, for retrieval of online data, e.g., ephemerides and orbits.  **Recommended**
* Python extensions for `oorb <https://github.com/oorb/oorb/>`__: For orbit
  transformations (`~sbpy.data.Orbit.oo_transform`) and propagations
  (`~sbpy.data.Orbit.oo_propagate`), as well as ephemerides calculations
  (`~sbpy.data.Ephem.from_oo`).
* `pyradex <https://github.com/keflavich/pyradex>`__: For non-LTE production
  rate calculations related to cometary activity (`~sbpy.activity.gas.NonLTE`).
* `scipy <https://scipy.org/>`__: 1.6 or later, for numerical integrations in `sbpy.activity.gas` and `sbpy.photometry`, among others.  **Recommended**
* `synphot <https://github.com/spacetelescope/synphot_refactor>`__ 1.1.1 or later, for calibration with respect to the Sun and Vega, filtering spectra through bandpasses.  **Recommended**
* `ginga <https://ejeschke.github.io/ginga/>`__ : To interactively enhance
  images of comets with the `~sbpy.imageanalysis.CometaryEnhancement` Ginga
  plugin.
* `photutils <https://photutils.readthedocs.io/en/stable/>`__: For centroiding within the Cometary Enhancements Ginga plugin.


Using pip
^^^^^^^^^

The latest stable version of `sbpy` can be installed with:

.. code-block:: bash

    $ pip install sbpy

Recommended dependencies may be installed via:

.. code-block:: bash

    $ pip install sbpy[recommended]

Most optional dependencies may be installed via:

.. code-block:: bash

    $ pip install sbpy[all]

`pyradex` and the ephemeris data for `oorb` are left for the user to install
manually.

The latest development version of `sbpy` can be installed using:

.. code-block:: bash

    $ pip install git+https://github.com/NASA-Planetary-Science/sbpy.git


Using conda
^^^^^^^^^^^

The latest stable version of `sbpy` can be installed with `Anaconda
<https://www.anaconda.com/>`__ via the `conda-forge <https://conda-forge.org/>`__
channel:

.. code-block:: bash

    $ conda install sbpy --channel=conda-forge

If you do not have the conda-forge channel available, add it and re-run the
installation command:

.. code-block:: bash

    $ conda config --add channels conda-forge
    $ conda install sbpy --channel=conda-forge


Using Git+Pip
^^^^^^^^^^^^^

This way of installing `sbpy` is recommended if you plan to contribute to the
module. The current development version of `sbpy` can be obtained from `GitHub
<https://github.com/NASA-Planetary-Science/sbpy>`__ using:

.. code-block:: bash

    $ git clone https://github.com/NASA-Planetary-Science/sbpy.git

This will create a new directory (``sbpy/``). In this directory, run:

.. code-block:: bash

    $ pip install .

As above, to install optional dependencies, instead use ``pip install .[all]``.

If you plan to work on the code and always want to use the latest version of
your code, we recommend installing in "editable" mode with the optional
dependences and the testing dependencies:

.. code-block:: bash

    $ pip install -e .[all,test]
