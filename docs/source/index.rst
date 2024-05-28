.. grasp2alm documentation master file, created by
   sphinx-quickstart on Mon May 27 15:48:41 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to grasp2alm's documentation!
=====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

This package supports the conversion from beam data calculated using `GRASP <https://www.ticra.com/software/grasp/>`_ for CMB experiments to spherical harmonic coefficients :math:`a_{\ell m}` based on the `HEALPix <https://healpix.sourceforge.io/>`_ framework.
The code is designed based on `Beam <https://github.com/zonca/planck-levelS/tree/master/Beam/>`_, which is part of `LevelS <https://github.com/zonca/planck-levelS>`_, the pipleline of the Planck experiment.

Installation
------------

You can use `pip <https://pypi.org/project/pip/>`_ by:

.. code-block:: python

   pip install grasp2alm


Or you can install it from the GitHub source by:

.. code-block:: python

   git clone https://github.com/yusuke-takase/grasp2alm
   cd grasp2alm
   pip install -e .

Tutorial
--------
.. toctree::
   :maxdepth: 1
   
   tutorial

Changelog
---------

Review the changes in each release in the `CHANGELOG on Github <https://github.com/yusuke-takase/grasp2alm/blob/main/CHANGELOG.rst>`_.

Reference
---------
.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
