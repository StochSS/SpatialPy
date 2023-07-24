Documentation for SpatialPy |release|
=====================================



SpatialPy is a Python 3 package for simulation of spatial deterministic/stochastic reaction-diffusion-advection problems embedded in Lagrangian reference frame particle based fluid dynamics domain

This package is intended to replace the PyURDME software https://github.com/pyurdme/pyurdme and will feature both a NSM solver for RDME simulation on static domains and a sSSA-SDPD particle based fluid dynamics solver as described in the publication "A hybrid smoothed dissipative particle dynamics (SDPD) spatial stochastic simulation algorithm (sSSA) for advection–diffusion–reaction problems" by Drawert, Jacob, Li, Yi, Petzold https://www.sciencedirect.com/science/article/pii/S0021999118307101



Getting a copy of SpatialPy
***************************

The latest version of SpatialPy can be found on `PyPI <https://pypi.org/project/spatialpy>`_.  The source code is available on `GitHub <https://github.com/StochSS/SpatialPy>`_.  SpatialPy is licensed under the GNU General Public License version 3.

.. raw:: html

      <p style="width: 80%; margin: auto; padding: 0.5em; border: 1px solid purple"><a href="https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform"><b>Please register as a user!</b></a>  SpatialPy's development is funded by NIH grant 2R01EB014877, and to continue support, we need to prove SpatialPy has users. <a href="https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform">Please fill out our short registration form!</a></p>

Examples
********
See our `Example Notebook - Start Here <https://github.com/StochSS/SpatialPy/blob/main/examples/Start_Here.ipynb>`_ for more information on how to build and simulate your models with SpatialPy. For an example of how to use SpatialPy to simulate a spatial stochastic reaction-diffusion system, see the `3D Cylinder Demo <https://github.com/StochSS/SpatialPy/blob/main/examples/3D_Cylinder_Demo.ipynb>`_. We also provide examples of how to use SpatialPy to simulate physics (`Gravity Demo <https://github.com/StochSS/SpatialPy/blob/main/examples/Gravity.ipynb>`_) and fluid flow (`Weir Model <https://github.com/StochSS/SpatialPy/blob/main/examples/Weir.ipynb>`_).


Reporting Issues
****************

If you find any problem with SpatialPy or this documentation, please report it using `the GitHub issue tracker <https://github.com/StochSS/SpatialPy/issues>`_ for the project.  You can also contact the main author, Dr. `Brian Drawert <https://github.com/briandrawert>`_, directly with questions and suggestions.


Documentation
*************

.. toctree::
   :maxdepth: 3
   :caption: API reference

   classes/spatialpy



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
