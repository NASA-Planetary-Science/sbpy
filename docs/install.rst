
Installation
------------

`sbpy` requires Python 3.5 or later - compatibility with Python 2.x is not
supported. We hence recommend that you install the latest version of
`Anaconda Python 3.x <https://www.anaconda.com/download/>`__ on your
system before installing `sbpy`. Make sure that Anaconda Python is
your default Python (this will be asked during the installation process).

Requirements
^^^^^^^^^^^^

`sbpy` has the following requirements (incomplete):

* Python 3.5 or later
* `numpy <https://www.numpy.org/>`__ 1.4.0 or later
* pytest 3.1 or later
* `astropy <https://www.astropy.org/>`__

`sbpy` also depends on the following packages for optional features (incomplete list):

* `astroquery <https://astroquery.readthedocs.io/en/latest/>`__ 0.3.9.dev5089 or later: For retrieval of online data, e.g., ephemerides and orbits.
* `scipy <https://scipy.org/>`__: For numerical integration of `activity.GasComa` distributions, e.g., in order to compute gas column density.
* `synphot <https://github.com/spacetelescope/synphot_refactor>`__: For calibration to Sun and Vega.
* `ginga <https://ejeschke.github.io/ginga/>`__: To use the ``CometaryEnhancements`` Ginga plug-in.
* `photutils <https://photutils.readthedocs.io/en/stable/>`__: For centroiding within ``CometaryEnhancements``.
* `oorb <https://github.com/oorb/oorb>`__: For orbit calculations that utilize ``pyoorb``.

Most requirements should be resolved during the installation process. However, we recommend to install the latest development version of `astroquery` using

.. code-block:: bash

    $ pip install --pre astroquery

Also, if you want to use `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`__, you will have to
install it using the instructions provided on that page.


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
