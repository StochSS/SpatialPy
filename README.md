<p align="left">
<img src="https://raw.githubusercontent.com/StochSS/SpatialPy/main/logo/SpatialPy_logo.png">
</p>

SpatialPy is a Python 3 package for simulation of spatial deterministic/stochastic reaction-diffusion-advection problems embedded in Lagrangian reference frame particle based fluid dynamics domain

This package is intended to replace the PyURDME software https://github.com/pyurdme/pyurdme and will feature both a NSM solver for RDME simulation on static domains and a sSSA-SDPD particle based fluid dynamics solver as described in the publication "A hybrid smoothed dissipative particle dynamics (SDPD) spatial stochastic simulation algorithm (sSSA) for advection–diffusion–reaction problems" by Drawert, Jacob, Li, Yi, Petzold https://www.sciencedirect.com/science/article/pii/S0021999118307101

<table><tr><td><b>
<img width="20%" align="right" src="https://raw.githubusercontent.com/SpatialPy/SpatialPy/develop/.graphics/stochss-logo.png">
<a href="https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform">PLEASE REGISTER AS A USER</a>, so that we can prove SpatialPy has many users when we seek funding to support development. SpatialPy is part of the <a href="http://www.stochss.org">StochSS</a> project.
</td></tr></table>

![PyPI - License](https://img.shields.io/pypi/l/spatialpy.svg?color=informational)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spatialpy.svg)
[![PyPI](https://img.shields.io/pypi/v/spatialpy.svg)](https://pypi.org/project/spatialpy)
![PyPI - Downloads](https://img.shields.io/pypi/dm/SpatialPy?color=informational&label=pypi%20downloads)

Table of contents
-----------------

- [Table of contents](#table-of-contents)
- [Installation](#installation)
  - [_Using PyPI_](#using-pypi)
  - [_Using the source code repository_](#using-the-source-code-repository)
- [Usage](#usage)
  - [_Simple example to illustrate the use of SpatialPy_](#simple-example-to-illustrate-the-use-of-spatialpy)
  - [_Docker environment_](#docker-environment)
  - [_Debugging_](#debugging)
  - [_Profiling_](#profiling)
- [Getting help](#getting-help)
- [Contributing](#contributing)
- [License](#license)
- [Authors and history](#authors-and-history)
- [Acknowledgments](#acknowledgments)

Installation
------------

SpatialPy can be installed on your computer using different methods, as described below.

### _Using PyPI_

On **Linux**, **macOS**, and **Windows** operating systems, you should be able to install SpatialPy with `pip`. Please review the official pip [documentation](https://pip.pypa.io/en/stable/installing/) for installation instructions and additional information.

Then, to install SpatialPy from the Python package repository, run the following command:
```sh
python3 -m pip install spatialpy --user --upgrade
```

### _Using the source code repository_

As an alternative to getting it from PyPI, you can instruct `pip` to install SpatialPy directly from the GitHub repository:
```sh
python3 -m pip install https://github.com/StochSS/SpatialPy/archive/main.zip --user --upgrade
```

As a final alternative, you can first use `git` to clone a copy of the SpatialPy source tree from the GitHub repository to your local computer disk, and then install SpatialPy using that copy:
```sh
git clone https://github.com/StochSS/SpatialPy.git
cd SpatialPy
python3 -m pip install  .  --user --upgrade
```

Usage
-----

SpatialPy provides simple object-oriented abstractions for defining a model of a biochemical system and simulating that model using efficient stochastic simulation algorithms.  The basic steps to use SpatialPy are:

1. Create a `SpatialPy.Model` containing molecular species, parameters, and reactions.
2. Invoke the model's `.run()` method.

The `run()` method can be customized using keyword arguments to select different solvers, random seed, data return type and more.  For more detailed examples on how to use SpatialPy, please see the Jupyter notebooks contained in the [examples](https://github.com/StochSS/SpatialPy/tree/main/examples) subdirectory.


### Docker environment

You can use Docker to create a repeatable environment for developing and debugging SpatialPy. The supplied Dockerfile starts a jupyter server with SpatialPy dependencies installed.

If you have Docker Compose: `docker-compose build && docker-compose up`

Otherwise:

```bash
docker build -t spatialpy:latest .
docker run -v ./:/home/jovyan/spatialpy -v ./tmp:/tmp -p 8888:8888 -p 5000:5000
```

The SpatialPy repo is mounted into /home/jovyan so you can import it in the usual way for development (see examples).

Any changes you make to your local codebase are reflected in the docker container. Note that you DO NOT need to restart the docker container or even re-import spatialpy to see source changes take effect in jupyter notebooks.

The `/tmp` directory is also mounted for easy access to build and result directories.

### Debugging

In order to compile the solver binary for use by the debugger, run `solver.compile()` with `debug=True`. This will inject the `-g` flag into the `gcc` command that compiles the solver, enabling gdb debug information.

You can invoke `solver.run_debugger()` anytime after you instantiate a solver in Python to start up a new session of gdbgui. The debugger will be available at http://127.0.0.1:5000.


### Profiling

To enable profiling, both `solver.compile()` and `solver.run()` need to be invoked with `profile=True`. If you don't run `solver.compile()` explicitly, invoking `solver.run()` with `profile=True` will run `compile()` correctly for you.

Getting help
------------

SpatialPy's [online documentation](https://stochss.github.io/SpatialPy) provides more details about using the software.  If you find any problem with SpatialPy or the documentation, please report it using the [GitHub issue tracker](https://github.com/StochSS/SpatialPy/issues) for this repository.  You can also contact Dr. [Brian Drawert](http://www.cs.unca.edu/~drawert) directly with questions and suggestions.

Contributing
------------

We would be happy to receive your help and participation with enhancing SpatialPy! Please follow the guidelines described in [CONTRIBUTING.md](https://github.com/StochSS/SpatialPy/tree/main/CONTRIBUTING.md).

New developments happen primarily in the [`develop`](https://github.com/StochSS/SpatialPy/commits/develop) branch.  New releases are put in the `main` branch.

<p align="center">

| Main Branch   | Develop Branch |
|:---------------:|:--------------:|
| [![Build Status](https://github.com/StochSS/SpatialPy/actions/workflows/run-tests.yml/badge.svg)](https://github.com/StochSS/SpatialPy/actions/workflows/run-tests.yml) | [![Build Status](https://github.com/StochSS/SpatialPy/actions/workflows/run-tests.yml/badge.svg?branch=develop)](https://github.com/StochSS/SpatialPy/actions/workflows/run-tests.yml)

License
-------

SpatialPy is licensed under the GNU General Public License version 3.  Please see the file [LICENSE](https://github.com/StochSS/SpatialPy/blob/main/LICENSE) for more information.

Authors and history
---------------------------

* [**Dr. Brian Drawert** ](https://github.com/briandrawert)
* [**Dr. Kevin Sanft**](https://github.com/kevinsanft)
* [**Sean Matthew**](https://github.com/seanebum)
* [**Evie Hilton**](https://github.com/eviehilton)
* [**Bryan Rumsey**](https://github.com/BryanRumsey)
* [**Bruno Jacob**](https://github.com/brunopjacob)
* [**Mason Kidwell**](https://github.com/makdl)
* [**Matthew Geiger**](https://github.com/popensesame)

Acknowledgments
---------------

This work has been funded by National Institutes of Health (NIH) NIBIB Award No. 2R01EB014877-04A1.

SpatialPy uses numerous open-source packages, without which it would have been effectively impossible to develop this software with the resources we had.  We want to acknowledge this debt.  In alphabetical order, the packages are:

* [Jupyter](https://jupyter.org) &ndash; web application for creating documents containing code, visualizations and narrative text
* [MatplotLib](https://matplotlib.org/index.html) &ndash; Python plotting library
* [Plotly](https://plot.ly/) &ndash; Graphing library for making interactive, publication-quality graphs
* [Numpy](https://www.numpy.org/) &ndash; the fundamental package for scientific computing with Python
* [Scipy](https://www.scipy.org/) &ndash; Python-based ecosystem of open-source software for mathematics, science, and engineering

Finally, we are grateful for institutional resources made available by the [University of North Carolina at Asheville](https://www.unca.edu), the [University of California at Santa Barbara](https://ucsb.edu), [Uppsala University](https://www.it.uu.se), and the [California Institute of Technology](https://www.caltech.edu).

<div align="center">
  <a href="https://www.nigms.nih.gov">
    <img width="100" height="100" src="https://raw.githubusercontent.com/StochSS/SpatialPy/develop/.graphics/US-NIH-NIGMS-Logo.png">
  </a>
  &nbsp;&nbsp;
  <a href="https://www.unca.edu">
    <img height="102" src="https://raw.githubusercontent.com/StochSS/SpatialPy/develop/.graphics/UNCASEAL_blue.png">
  </a>
  &nbsp;&nbsp;
  <a href="https://www.ucsb.edu">
    <img height="108" src="https://raw.githubusercontent.com/StochSS/SpatialPy/develop/.graphics/ucsb-seal-navy.jpg">
  </a>
  &nbsp;&nbsp;
  <a href="https://www.it.uu.se">
    <img height="115" src="https://raw.githubusercontent.com/StochSS/SpatialPy/develop/.graphics/uppsala-universitet-logo-svg-vector.png">
  </a>
  &nbsp;&nbsp;
  <a href="https://www.caltech.edu">
    <img width="115" src="https://raw.githubusercontent.com/StochSS/SpatialPy/develop/.graphics/caltech-round.png">
  </a>
</div>
