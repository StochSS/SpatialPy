<p align="left">
<img src="https://raw.githubusercontent.com/spatialpy/SpatialPy/develop/logo/SpatialPy_logo.png">
</p>

SpatialPy is a Python 3 package for simulation of spatial deterministic/stochastic reaction-diffusion-advection problems embedded in Lagrangian reference frame particle based fluid dynamics domain

This package is intended to replace the PyURDME software https://github.com/pyurdme/pyurdme and will feature both a NSM solver for RDME simulation on static domains and a sSSA-SDPD particle based fluid dynamics solver as described in the publication "A hybrid smoothed dissipative particle dynamics (SDPD) spatial stochastic simulation algorithm (sSSA) for advection–diffusion–reaction problems" by Drawert, Jacob, Li, Yi, Petzold https://www.sciencedirect.com/science/article/pii/S0021999118307101

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

<table><tr><td><b>
<img width="20%" align="right" src="https://raw.githubusercontent.com/GillesPy2/GillesPy2/develop/.graphics/stochss-logo.png">
<a href="https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform">PLEASE REGISTER AS A USER</a>, so that we can prove SpatialPy has many users when we seek funding to support development. SpatialPy is part of the <a href="http://www.stochss.org">StochSS</a> project.
</td></tr></table>

![GitHub](https://img.shields.io/github/license/spatialpy/SpatialPy?color=informational)
