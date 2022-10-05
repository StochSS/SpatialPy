# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2022 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy
from collections import OrderedDict
import scipy
from spatialpy.core.spatialpyerror import ModelError, SimulationError, SimulationTimeout
from spatialpy.solvers.solver import Solver

class HybridSolver(Solver):
    """
    SpatialPy solver object, implemented in python.  Uses the tau-hybrid algorithm.

    :param model: Target model of solver simulation.
    :type model: spatialpy.core.model.Model

    :param debug_level: Target level of debugging.
    :type debug_level: int
    """
    def __init__(self, model, debug=False):
        from spatialpy.core.model import Model # pylint: disable=import-outside-toplevel
        if not (isinstance(model, Model) or type(model).__name__ == 'Model'):
            raise SimulationError("Model must be of type spatialpy.Model.")
#        if not issubclass(self.__class__, Solver):
#            print(issubclass(self.__class__, Solver))
#            print(Solver)
#            print(self.__class__)
#            raise SimulationError("Solver classes must be a subclass of spatialpy.Solver.")

        self.model = model
        self.is_compiled = False
        self.debug = debug
        if not self.is_compiled:
            self.compile(debug=debug)

    def __find_neighbors(self,vox_i):
        neighbors = OrderedDict()
        x_i = self.model.domain.coordinates()[vox_i,:]
        for j in range(self.model.domain.coordinates().shape[0]):
            if vox_i == j: continue
            r = scipy.linalg.norm( self.model.domain.coordinates()[j,:] - x_i )
            if r < self.h:
                neighbors[j]=r
        return neighbors
            

    def compile(self, debug=False):
        """
        Compile the model.

        :param debug: If True, will print additional build debugging
        :type debug: bool

        :param profile: If True, will print additional profiling information
        :type profile: bool

        :raises SimulationError: Failed to compile
        """
        stoich_matrix, dep_graph = self.model.compile_prep()

        equations = OrderedDict()
        def _add_equation(key, term):
            if key not in equations:
                equations[key] = []
            equations[key].append(term)

        if self.h is None:
            self.h = self.model.domain.find_h()
        h = self.h

        aD = 105.0/16*numpy.pi*h*h*h # 3D

        for i in range(self.model.domain.coordinates().shape[0]):
            neighbors = self.__find_neighbors(i)
            for k in range(3):
                _add_equation(f"x[{i},{k}]", f"v[{i},{k}]")
                for j,r in neighbors.items():
                    _add_equation(f"rho[{i}]",f"-{self.model.domain.mass[i]}*(v[{i},{k}]-v[{j},{k}])*(x[{i},{k}]-x[{j},{k}]/{r}*dW ")
                #_add_equation["v[i,0]", 
        
        for k,v in equations.items():
            print(f"d/dt[ {k} ] = {' + '.join(v)}")



    def run(self, number_of_trajectories=1, seed=None, timeout=None,
                number_of_threads=None, debug=False, verbose=True):
        """
        Run one simulation of the model.

        :param number_of_trajectories: How many trajectories should be simulated.
        :type number_of_trajectories: int

        :param seed: the random number seed (incremented by one for multiple runs).
        :type seed: int

        :param timeout: maximum number of seconds the solver can run.
        :type timeout: int

        :param number_of_threads: the number threads the solver will use.
        :type number_of_threads: int

        :param debug: run in debug mode
        :type debug: bool

        :param verbose: If true, prints addtional data to console
        :type verbose: bool

        :returns: A SpatialPy Result object containing spatial and time series data from simulation.
        :rtype: spatialpy.Result.Result | list(spatialpy.Result.Result)

        :raises SimulationTimeout: Simulation exceeded timeout.
        :raises SimulationError: Simulation execution failed.
        """
        from spatialpy.core.result import Result # pylint: disable=import-outside-toplevel
        if number_of_trajectories > 1:
            result_list = []
        # Check if compiled, call compile() if not.
        if not self.is_compiled:
            self.compile(debug=debug)

