'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2019 - 2022 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#This module defines a model that simulates a discrete, stoachastic, mixed biochemical reaction network in python.

import math
from collections import OrderedDict

import numpy
import scipy

from spatialpy.solvers.build_expression import BuildExpression

from spatialpy.core.spatialpyerror import ModelError

def export_StochSS(spatialpy_model, filename=None, return_stochss_model=False):
    """
    SpatialPy model to StochSS converter

    :param spatailpy_model: SpatialPy model to be converted to StochSS
    :type spatialpy_model: spatialpy.Model

    :param filename: Path to the exported stochss model
    :type filename: str

    :param return_stochss_model: Whether or not to return the model
    :type return_stochss_model: bool

    :returns: Filename for JSON-formatted .smdl file for use with StochSS platform.
    :rtype: string
    """
    try:
        from spatialpy.stochss.stochss_export import export # pylint: disable=import-outside-toplevel
    except ImportError as err:
        raise ImportError('StochSS export conversion not imported successfully') from err

    return export(spatialpy_model, path=filename, return_stochss_model=return_stochss_model)

class Model():
    """
    Representation of a spatial biochemical model.

    :param name: Name of the model
    :type name: str
    """
    reserved_names = ['vol', 't']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']

    def __init__(self, name="spatialpy"):
        # The name that the model is referenced by (should be a String)
        self.name = name

        ######################
        # Dictionaries with Species, Reactions and Parameter objects.
        # Species, Reaction, and Parameter names are used as keys.
        self.listOfParameters = OrderedDict()
        self.listOfSpecies    = OrderedDict()
        self.listOfReactions  = OrderedDict()

        ######################
        # Dict that holds flattended parameters and species for
        # evaluation of expressions in the scope of the model.
        self.namespace = OrderedDict([])
        self.species_map = {}

        ######################
        self.domain = None
        self.listOfDiffusionRestrictions = {}
        self.listOfDataFunctions = []
        self.listOfInitialConditions = []
        self.listOfBoundaryConditions = []

        ######################
        self.staticDomain = True
        self.enable_rdme = True
        self.enable_pde = True

        ######################
        self.tspan = None
        self.timestep_size = None
        self.num_timesteps = None
        self.output_freq = None
        self.output_steps = None

        ######################
        # Expression utility used by the solver
        # Should be None until the solver is compiled
        self.expr = None
        self.u0 = None

    def __str__(self):
        try:
            self.__update_diffusion_restrictions()
        except Exception:
            pass
        self.__resolve_parameters()
        divider = f"\n{'*'*10}\n"

        def decorate(header):
            return f"\n{divider}{header}{divider}"

        print_string = f"{self.name}"
        if len(self.listOfSpecies):
            print_string += decorate("Species")
            for _, species in self.listOfSpecies.items():
                print_string += f"\n{str(species)}"
        if self.listOfInitialConditions:
            print_string += decorate("Initial Conditions")
            for initial_condition in self.listOfInitialConditions:
                print_string += f"\n{str(initial_condition)}"
        if len(self.listOfDiffusionRestrictions):
            print_string += decorate("Diffusion Restrictions")
            for species, types in self.listOfDiffusionRestrictions.items():
                print_string += f"\n{species.name} is restricted to: {str(types)}"
        if len(self.listOfParameters):
            print_string += decorate("Parameters")
            for _, parameter in self.listOfParameters.items():
                print_string += f"\n{str(parameter)}"
        if len(self.listOfReactions):
            print_string += decorate("Reactions")
            for _, reaction in self.listOfReactions.items():
                print_string += f"\n{str(reaction)}"
        print(print_string)

        if self.domain is not None:
            print(decorate("Domain"))
            print(f"\n{str(self.domain)}")

        return ""

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return self.listOfParameters == other.listOfParameters and \
            self.listOfSpecies == other.listOfSpecies and \
            self.listOfReactions == other.listOfReactions and \
            self.name == other.name

    def __check_if_complete(self):
        if self.timestep_size is None or self.num_timesteps is None:
            raise ModelError("The model's timespan is not set.  Use 'timespan()' or 'set_timesteps()'.")
        if self.domain is None:
            raise ModelError("The model's domain is not set.  Use 'add_domain()'.")

    def __problem_with_name(self, name):
        if name in Model.reserved_names:
            errmsg = f'Name "{name}" is unavailable. It is reserved for internal SpatialPy use. '
            errmsg += f'Reserved Names: ({Model.reserved_names}).'
            raise ModelError(errmsg)
        if name in self.listOfSpecies:
            raise ModelError(f'Name "{name}" is unavailable. A species with that name exists.')
        if name in self.listOfParameters:
            raise ModelError(f'Name "{name}" is unavailable. A parameter with that name exists.')
        if name.isdigit():
            raise ModelError(f'Name "{name}" is unavailable. Names must not be numeric strings.')
        for special_character in Model.special_characters:
            if special_character in name:
                errmsg = f'Name "{name}" is unavailable. '
                errmsg += f'Names must not contain special characters: {Model.special_characters}.'
                raise ModelError(errmsg)

    def __apply_initial_conditions(self):
        # initalize
        num_spec = self.get_num_species()
        num_vox = self.domain.get_num_voxels()
        self.u0 = numpy.zeros((num_spec, num_vox))
        # apply initial condition functions
        for init_cond in self.listOfInitialConditions:
            init_cond.apply(self)

    def __create_dependency_graph(self):
        # We cannot safely generate a dependency graph (without attempting to analyze the
        # propensity string itself) if the model contains custom propensities.
        mass_action_model = True
        for reaction in self.listOfReactions.values():
            if not reaction.massaction:
                raw_graph = numpy.ones((self.get_num_reactions(),
                    self.get_num_reactions() + self.get_num_species()))
                mass_action_model = False

        if mass_action_model:
            raw_graph = numpy.zeros((self.get_num_reactions(),
                self.get_num_reactions() + self.get_num_species()))
            species_map = self.species_map

            involved_species = []
            reactants = []
            for reaction in self.listOfReactions.values():
                temp = []
                temp2 = []
                for species in reaction.reactants:
                    temp.append(species_map[species])
                    temp2.append(species_map[species])
                for species in reaction.products:
                    temp.append(species_map[species])
                involved_species.append(temp)
                reactants.append(temp2)

            species_to_reactions = []
            for species in self.listOfSpecies.values():
                temp = []
                for j, x in enumerate(reactants):
                    if species_map[species] in x:
                        temp.append(j)
                species_to_reactions.append(temp)

            reaction_to_reaction = []
            for reaction in self.listOfReactions.values():
                temp = []
                for species in reaction.reactants:
                    if species_to_reactions[species_map[species]] not in temp:
                        temp = temp + species_to_reactions[species_map[species]]

                for species in reaction.products:
                    if species_to_reactions[species_map[species]] not in temp:
                        temp = temp + species_to_reactions[species_map[species]]

                temp = list(set(temp))
                reaction_to_reaction.append(temp)

            # Populate raw_graph
            for j, spec in enumerate(species_to_reactions):
                for species in spec:
                    raw_graph[species, j] = 1

            for i, reac in enumerate(reaction_to_reaction):
                for reaction in reac:
                    raw_graph[reaction, self.get_num_species() + i] = 1
        try:
            dep_graph = scipy.sparse.csc_matrix(raw_graph)
        except Exception:
            dep_graph = raw_graph

        return dep_graph

    def __create_stoichiometric_matrix(self):
        if self.get_num_reactions() > 0:
            raw_matrix = numpy.zeros((self.get_num_species(), self.get_num_reactions()))
            for i, reaction in enumerate(self.listOfReactions.values()):
                reactants = reaction.reactants
                products  = reaction.products

                for species in reactants:
                    raw_matrix[self.species_map[species], i] -= reactants[species]
                for species in products:
                    raw_matrix[self.species_map[species], i] += products[species]

            matrix = scipy.sparse.csc_matrix(raw_matrix)
        else:
            matrix = numpy.zeros((self.get_num_species(), self.get_num_reactions()))

        return matrix

    def __get_expression_utility(self):
        base_namespace = {
            **{name: name for name in math.__dict__},
            **self.sanitized_species_names(),
            **self.sanitized_parameter_names(),
            **self.sanitized_data_function_names(),
            **{name: name for name in self.reserved_names}
        }
        self.expr = BuildExpression(namespace=base_namespace, blacklist=["="], sanitize=True)

    def _ipython_display_(self, use_matplotlib=False):
        if self.domain is None:
            print(self)
        else:
            self.domain.plot_types(width="auto", height="auto", use_matplotlib=use_matplotlib)

    def __resolve_parameters(self):
        self.update_namespace()
        for param in self.listOfParameters:
            self.listOfParameters[param]._evaluate(self.namespace) # pylint: disable=protected-access

    def __update_diffusion_restrictions(self):
        for species in self.listOfSpecies.values():
            if isinstance(species.restrict_to, list):
                self.listOfDiffusionRestrictions[species] = species.restrict_to
            elif isinstance(species.restrict_to, str):
                self.listOfDiffusionRestrictions[species] = [species.restrict_to]

    def compile_prep(self):
        """
        Make sure all paramters are evaluated to scalars, update the models diffusion restrictions,
        create the models expression utility, and generate the domain list of type ids in preperation
        of compiling the simulation files.

        :returns: The stoichiometric and dependency_graph
        :rtype: tuple

        :raises ModelError: Timestep size exceeds output frequency or Model is missing a domain
        """
        if self.timestep_size is None:
            self.timestep_size = 1e-5
        if self.output_freq < self.timestep_size:
            raise ModelError("Timestep size exceeds output frequency.")

        self.__check_if_complete()

        self.domain.compile_prep()
        
        self.__update_diffusion_restrictions()
        self.__apply_initial_conditions()
        self.__resolve_parameters()
        self.__get_expression_utility()
        stoich_matrix = self.__create_stoichiometric_matrix()
        dep_graph = self.__create_dependency_graph()

        return stoich_matrix, dep_graph

    def add_species(self, obj):
        """
        Adds a species, or list of species to the model. Will return the added object upon success.

        :param obj: The species or list of species to be added to the model object.
        :type obj: spatialpy.Model.Species | list(spatialpy.Model.Species

        :returns: Species object which was added to the model.
        :rtype: spatialpy.Species | list(spatialpy.Species)

        :raises ModelError: If obj is not a spatialpy.Species
        """
        from spatialpy.core.species import Species # pylint: disable=import-outside-toplevel
        if isinstance(obj, list):
            for species in obj:
                self.add_species(species)
        elif isinstance(obj, Species) or type(obj).__name__ == 'Species':
            self.__problem_with_name(obj.name)
            self.species_map[obj] = len(self.listOfSpecies)
            self.listOfSpecies[obj.name] = obj
        else:
            raise ModelError("Unexpected parameter for add_species. Parameter must be Species or list of Species.")
        return obj

    def delete_all_species(self):
        """
        Remove all species from model.listOfSpecies.
        """
        self.listOfSpecies.clear()

    def delete_species(self, obj):
        """
        Remove a Species from model.listOfSpecies.

        :param obj: Species object to be removed
        :type obj: spatialpy.Model.Species
        """
        self.listOfSpecies.pop(obj) # raises key error if param is missing

    def get_all_species(self):
        """
        Returns a dictionary of all species in the model using names as keys.

        :returns: A dictionary of all species in the form of {"species_name":Species_object}
        :rtype: dict
        """
        return self.listOfSpecies

    def get_num_species(self):
        """
        Returns total number of different species contained in the model.

        :returns: Number of different species in the model.
        :rtype: int
        """
        return len(self.listOfSpecies)

    def get_species(self, sname):
        """
        Returns target species from model as object.

        :param sname: name of species to be returned.
        :type sname: str

        :returns: The Species objected represented by given 'sname'
        :rtype: spatialpy.Model.Species

        :raises ModelError: if the model does not contain the requested species
        """
        try:
            return self.listOfSpecies[sname]
        except KeyError as err:
            raise ModelError(f"No species named {sname}") from err

    def sanitized_species_names(self):
        """
        Generate a dictionary mapping user chosen species names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user species names to their internal SpatialPy notation.
        :rtype: dict
        """
        species_name_mapping = OrderedDict([])
        for i, name in enumerate(self.listOfSpecies.keys()):
            species_name_mapping[name] = f'x[{i}]'
        return species_name_mapping

    def add_parameter(self,params):
        """
        Add Parameter(s) to model.listOfParameters. Input can be either a single
        Parameter object or a list of Parameter objects.

        :param params: Parameter object or list of Parameters to be added.
        :type params: spatialpy.Model.Parameter | list(spatialpy.Model.Parameter)

        :returns: Parameter object which has been added to the model.
        :rtype: spatialpy.Parameter | list(spatialpy.Parameter)

        :raises ModelError: if obj is not a spatialpy.Parameter
        """
        from spatialpy.core.parameter import Parameter # pylint: disable=import-outside-toplevel
        if isinstance(params,list):
            for param in params:
                self.add_parameter(param)
        elif isinstance(params, Parameter) or  type(params).__name__ == 'Parameter':
            self.__problem_with_name(params.name)
            self.update_namespace()
            params._evaluate(self.namespace) # pylint: disable=protected-access
            self.listOfParameters[params.name] = params
        else:
            errmsg = f"Parameter '{params.name}' needs to be of type '{Parameter}', it is of type '{params}'"
            raise ModelError(errmsg)
        return params

    def delete_all_parameters(self):
        """
        Remove all parameters from model.listOfParameters.
        """
        self.listOfParameters.clear()

    def delete_parameter(self, obj):
        """
        Remove a Parameter from model.listOfParameters.

        :param obj: Parameter object to be removed
        :type obj: spatialpy.Model.Parameter
        """
        self.listOfParameters.pop(obj)

    def get_all_parameters(self):
        """
        Return a dictionary of all model parameters, indexed by name.

        :returns: A dictionary of all model parameters in the form {'param_name':param_obj}
        :rtype: dict
        """
        return self.listOfParameters

    def get_parameter(self, pname):
        """
        Return the Parameter object from model associated with 'pname'

        :param pname: Name of parameter to be removed
        :type pname: spatialpy.Model.Parameter

        :returns: The Parameter object represented in the model by 'pname'
        :rtype: Spatialpy.Model.Parameter

        :raises ModelError: No parameter named {pname}
        """
        try:
            return self.listOfParameters[pname]
        except KeyError as err:
            raise ModelError(f"No parameter named {pname}") from err

    def sanitized_parameter_names(self):
        """
        Generate a dictionary mapping user chosen parameter names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user parameter names to their internal SpatialPy notation.
        :rtype: dict
        """
        parameter_name_mapping = OrderedDict()
        for i, name in enumerate(self.listOfParameters.keys()):
            if name not in parameter_name_mapping:
                parameter_name_mapping[name] = f'P{i}'
        return parameter_name_mapping

    def add_reaction(self, reacs):
        """
        Add Reaction(s) to the model. Input can be single instance, a list of instances
        or a dict with name, instance pairs.

        :param reacs: Reaction or list of Reactions to be added.
        :type reacs: spatialpy.Model.Reaction | list(spatialpy.Model.Reaction)

        :returns: The Reaction object(s) added to the model
        :rtype: spatialpy.Model.Reaction

        :raises ModelError: Invalid input/reaction to add_reaction()
        """
        from spatialpy.core.reaction import Reaction # pylint: disable=import-outside-toplevel
        if isinstance(reacs, list):
            for reaction in reacs:
                self.add_reaction(reaction)
        elif isinstance(reacs, Reaction) or type(reacs).__name__ == "Reaction":
            self.__problem_with_name(reacs.name)
            reacs.initialize(self)
            self.listOfReactions[reacs.name] = reacs
        else:
            raise ModelError("add_reaction() takes a spatialpy.Reaction object or list of objects")
        return reacs

    def delete_all_reactions(self):
        """
        Remove all reactions from model.listOfReactions.
        """
        self.listOfReactions.clear()

    def delete_reaction(self, obj):
        """
        Remove reaction from model.listOfReactions

        :param obj: Reaction to be removed.
        :type obj: spatialpy.Model.Reaction
        """
        self.listOfReactions.pop(obj)

    def get_all_reactions(self):
        """
        Returns a dictionary of all model reactions using names as keys.

        :returns: A dictionary of reactions in the form of {'react_name':react_obj}
        :rtype: dict
        """
        return self.listOfReactions

    def get_num_reactions(self):
        """
        Returns the number of reactions in this model.

        :returns: The total number of different reactions in the model.
        :rtype: int
        """
        return len(self.listOfReactions)

    def get_reaction(self, rname):
        """
        Retrieve a reaction object from the model by name

        :param rname: name of Reaction to retrieve
        :type rname: str

        :returns: The Reaction Object in the model represented by 'rname'
        :rtype: spatialpy.Model.Reaction

        :raises ModelError: Could not find reaction
        """
        try:
            return self.listOfReactions[rname]
        except KeyError as err:
            raise ModelError(f"No reaction named {rname}") from err

    def run(self, number_of_trajectories=1, seed=None, timeout=None,
            number_of_threads=None, debug_level=0, debug=False, profile=False):
        """
        Simulate the model. Returns a result object containing simulation results.

        :param number_of_trajectories: How many trajectories should be run.
        :type number_of_trajectories: int

        :param seed: The random seed given to the solver.
        :type seed: int

        :param timeout: Number of seconds for simulation to run.  Simulation will be
                        killed upon reaching timeout.
        :type timeout: int

        :param number_of_threads: The number threads the solver will use.
        :type number_of_threads: int

        :param debug_level: Level of output from the solver: 0, 1, or 2. Default: 0.
        :type debug_level: int

        :param debug: Optional flag to print out additional debug info during simulation.
        :type debug: bool

        :param profile: Optional flag to print out addtional performance profiling for simulation.
        :type profile: bool

        :returns: A SpatialPy Result object containing simulation data.
        :rtype: spatialpy.Result.Result
        """
        from spatialpy.solvers.solver import Solver # pylint: disable=import-outside-toplevel

        sol = Solver(self, debug_level=debug_level)

        return sol.run(number_of_trajectories=number_of_trajectories, seed=seed, timeout=timeout,
                       number_of_threads=number_of_threads, debug=debug, profile=profile)



    def set_timesteps(self, output_interval, num_steps, timestep_size=None):
        """"
        Set the simlation time span parameters.

        :param output_interval: size of each output timestep in seconds
        :type output_interval:  float

        :param num_steps: total number of steps to take. Note: the number of output times will be num_steps+1 \
        as the first output will be at time zero.
        :type num_steps: int

        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: float

        :raises ModelError: Incompatible combination of timestep_size and output_interval
        """
        if timestep_size is not None:
            self.timestep_size = timestep_size
        if self.timestep_size is None:
            self.timestep_size = output_interval

        self.output_freq = output_interval/self.timestep_size
        if self.output_freq < self.timestep_size:
            raise ModelError("Timestep size exceeds output frequency.")

        self.num_timesteps = math.ceil(num_steps * self.output_freq)

        # array of step numbers corresponding to the simulation times in the timespan
        output_steps = numpy.arange(0, self.num_timesteps + self.timestep_size, self.output_freq)
        self.output_steps = numpy.unique(numpy.round(output_steps).astype(int))
        self.tspan = numpy.zeros((self.output_steps.size), dtype=float)
        for i, step in enumerate(self.output_steps):
            self.tspan[i] = step*self.timestep_size

    def timespan(self, time_span, timestep_size=None):
        """
        Set the time span of simulation. The SSA-SDPD engine does not support
        non-uniform timespans.

        :param tspan: Evenly-spaced list of times at which to sample the species populations during the simulation.
        :type tspan: numpy.ndarray

        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: float

        :raises ModelError: non-uniform timespan not supported
        """
        items_diff = numpy.diff(time_span)
        items = map(lambda x: round(x, 10), items_diff)
        isuniform = (len(set(items)) == 1)

        if isuniform:
            self.set_timesteps(items_diff[0], len(items_diff), timestep_size=timestep_size)
        else:
            raise ModelError("Only uniform timespans are supported")

    def add_domain(self, domain):
        """
        Add a spatial domain to the model

        :param domain: The Domain object to be added to the model
        :type domain: spatialpy.Domain.Domain

        :raises ModelError: Invalid Domain object
        """
        from spatialpy.core.domain import Domain # pylint: disable=import-outside-toplevel
        if not isinstance(domain,Domain) and type(domain).__name__ != 'Domain':
            raise ModelError("Unexpected parameter for add_domain. Parameter must be a Domain.")

        self.domain = domain

    def add_data_function(self, data_function):
        """
        Add a scalar spatial function to the simulation. This is useful if you have a
        spatially varying in put to your model. Argument is a instances of subclass of the
        spatialpy.DataFunction class. It must implement a function 'map(x)' which takes a
        the spatial positon 'x' as an array, and it returns a float value.

        :param data_function: Data function to be added.
        :type data_function: spatialpy.DataFunction

        :returns: DataFunction object(s) added tothe model.
        :rtype: spatialpy.DataFunction | list(spatialpy.DataFunction)

        :raises ModelError: Invalid DataFunction
        """
        from spatialpy.core.datafunction import DataFunction # pylint: disable=import-outside-toplevel
        if isinstance(data_function, list):
            for data_fn in data_function:
                self.add_data_function(data_fn)
        elif isinstance(data_function, DataFunction) or type(data_function).__name__ == 'DataFunction':
            self.__problem_with_name(data_function.name)
            self.listOfDataFunctions.append(data_function)
        else:
            errmsg = "Unexpected parameter for add_data_function. "
            errmsg += "Parameter must be DataFunction or list of DataFunctions."
            raise ModelError(errmsg)
        return data_function

    def add_initial_condition(self, init_cond):
        """
        Add an initial condition object to the initialization of the model.

        :param init_cond: Initial condition to be added
        :type init_cond: spatialpy.InitialCondition.InitialCondition
        """
        self.listOfInitialConditions.append(init_cond)

    def add_boundary_condition(self, bound_cond):
        """
        Add an BoundaryCondition object to the model.

        :param bound_cond: Boundary condition to be added
        :type bound_cond: spatialpy.BoundaryCondition
        """
        bound_cond.model = self
        self.listOfBoundaryConditions.append(bound_cond)

    def update_namespace(self):
        """
        Create a dict with flattened parameter and species objects.
        """
        for param in self.listOfParameters:
            self.namespace[param]=self.listOfParameters[param].value

    def sanitized_data_function_names(self):
        """
        Generate a dictionary mapping user chosen data function names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user data function names to their internal SpatialPy notation.
        :rtype: dict
        """
        data_fn_name_mapping = OrderedDict([])
        for i, data_fn in enumerate(self.listOfDataFunctions):
            data_fn_name_mapping[data_fn.name] = f'data_fn[{i}]'
        return data_fn_name_mapping
