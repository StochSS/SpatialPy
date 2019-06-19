.. SpatialPy documentation master file, created by
   sphinx-quickstart on Mon Apr 21 11:09:45 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SpatialPy's documentation!
===================================

SpatialPy is a general software framework for modeling and simulation of stochastic reaction-diffusion processes on unstructured, tetrahedral (3D) and triangular (2D) meshes. Unstructured meshes allow for a more flexible handling of complex geometries compared to structured, Cartesian meshes. The current core simulation algorithm is based on the mesoscopic reaction-diffusion master equation (RDME) model.

SpatialPy was originally based on the PyURDME software package, and before that the URDME software. It has been rewritted to use the Smoothed Disapative Particles Dynamics formunlation for the discretization of the diffusion equation.  The core simulation routines are implemented in C, and requires GCC for compilation.  The default solver is an efficient implementation of the Next Subvolume Method (NSM).

Contents:

.. toctree::
   :maxdepth: 4

   install
   pyurdme
   examples


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

