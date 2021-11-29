'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

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

import numpy
import scipy
import warnings
import math
import uuid
from collections import OrderedDict

from spatialpy.core.Species import Species
from spatialpy.core.Parameter import Parameter
from spatialpy.core.Reaction import Reaction
from spatialpy.core.DataFunction import DataFunction
from spatialpy.core.Domain import Domain
from spatialpy.solvers.Solver import Solver
from spatialpy.solvers.build.expression import Expression

from spatialpy.core.spatialpyError import *


def export_StochSS(spatialpy_model, filename=None, return_stochss_model=False):
    """
    SpatialPy model to StochSS converter

        :param spatailpy_model: SpatialPy model to be converted to StochSS
        :type spatialpy_model: spatialpy.Model
        :param filename: Path to the exported stochss model
        :type filename: str
        :param return_stochss_model: Whether or not to return the model
        :type return_stochss_model: bool
    """
    try:
        from spatialpy.stochss.StochSSexport import export
    except ImportError:
        raise ImportError('StochSS export conversion not imported successfully')

    return export(spatialpy_model, path=filename, return_stochss_model=return_stochss_model)


class Model():
    """ Representation of a spatial biochemical model.

            :param name: Name of the model
            :type name: str

    """
    reserved_names = ['vol', 't']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']



    def __init__(self, name="spatialpy"):
        """ Create an empty SpatialPy model. """

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
        self.listOfTypeIDs = [1] # starts with type '1'
        self.listOfDiffusionRestrictions = {}
        self.listOfDataFunctions = []
        self.listOfInitialConditions = []
        self.listOfBoundaryConditions = []

        ######################
        self.staticDomain = True
        self.enable_rdme = True
        self.enable_pde = True #TODO

        ######################
        self.tspan = None
        self.timestep_size = None
        self.num_timesteps = None
        self.output_freq = None

        ######################
        # Expression utility used by the solver
        # Should be None until the solver is compiled
        self.expr = None


    def __str__(self):
        self._resolve_parameters()
        divider = f"\n{'*'*10}\n"

        def decorate(header):
            return f"\n{divider}{header}{divider}"

        print_string = f"{self.name}"
        if len(self.listOfSpecies):
            print_string += decorate("Species")
            for _, species in self.listOfSpecies.items():
                print_string += f"\n{str(species)}"
        if len(self.listOfInitialConditions):
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


    def _ipython_display_(self, use_matplotlib=False):
        if self.domain is None:
            print(self)
        else:
            self.domain.plot_types(width="auto", height="auto", use_matplotlib=use_matplotlib)


    def get_expression_utility(self):
        """
        Create a new expression.
        """
        base_namespace = {
            **{name: name for name in math.__dict__.keys()},
            **self.sanitized_species_names(),
            **self.sanitized_parameter_names(),
            **self.sanitized_data_function_names(),
            **{name: name for name in self.reserved_names}
        }
        self.expr = Expression(namespace=base_namespace, blacklist=["="], sanitize=True)


    def run(self, number_of_trajectories=1, seed=None, timeout=None, number_of_threads=None, debug_level=0, debug=False, profile=False):
        """ Simulate the model. Returns a result object containing simulation results.

            :param number_of_trajectories: How many trajectories should be run.
            :type number_of_trajectories: int
            :param seed: The random seed given to the solver.
            :type seed: int
            :param number_of_threads: The number threads the solver will use.
            :type number_of_threads: int
            :param debug_level: Level of output from the solver: 0, 1, or 2. Default: 0.
            :type debug_level: int

            :rtype: spatialpy.Result.Result
        """
        self.__check_if_complete()

        sol = Solver(self, debug_level=debug_level)

        return sol.run(number_of_trajectories=number_of_trajectories, seed=seed, timeout=timeout,
                       number_of_threads=number_of_threads, debug=debug, profile=profile)



    def __check_if_complete(self):
        """ Check if the model is complete, otherwise raise an approprate exception.
        Raises:
            ModelError
        """
        if self.timestep_size is None or self.num_timesteps is None:
            raise ModelError("The model's timespan is not set.  Use 'timespan()' or 'set_timesteps()'.")
        if self.domain is None:
            raise ModelError("The model's domain is not set.  Use 'add_domain()'.")
            
    def set_timesteps(self, output_interval, num_steps, timestep_size=None):
        """" Set the simlation time span parameters
        :param output_interval: size of each output timestep in seconds
        :type output_interval:  float
        :param num_steps: total number of steps to take. Note: the number of output times will be num_steps+1 as the first
          output will be at time zero.
        :type num_steps: int
        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: float
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
        sim_steps = numpy.arange(0, self.num_timesteps + self.timestep_size, self.timestep_size)
        self.tspan = numpy.zeros((self.output_steps.size), dtype=float)
        for i, step in enumerate(self.output_steps):
            self.tspan[i] = sim_steps[step]
    
    def timespan(self, time_span, timestep_size=None):
        """
        Set the time span of simulation. The SSA-SDPD engine does not support
        non-uniform timespans.

        :param tspan: Evenly-spaced list of times at which to sample the species populations during the simulation.
        :type tspan: numpy.ndarray
        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: float
        """

        items_diff = numpy.diff(time_span)
        items = map(lambda x: round(x, 10), items_diff)
        isuniform = (len(set(items)) == 1)

        if isuniform:
            self.set_timesteps(items_diff[0], len(items_diff), timestep_size=timestep_size)
        else:
            raise ModelError("Only uniform timespans are supported")



    def set_type(self, geometry_ivar, type_id, vol=None, mass=None, nu=None, rho=None, c=None, fixed=False):
        """ Add a type definition to the model.  By default, all regions are set to
        type 0. Returns the number of domain points that were tagged with this type_id

            :param geometry_ivar: an instance of a 'spatialpy.Geometry' subclass.  The 'inside()' method
                       of this object will be used to assign type_id to points.
            :type geometry_ivar: spatialpy.Geometry.Geometry
            :param type_id: (usually an int) the identifier for this type
            :type type_id: int
            :param vol: The volume of each particle in the type
            :type vol: float
            :param mass: The mass of each particle in the type
            :type mass: float
            :param rho: The density of each particle in the type
            :type rho: float
            :param nu: The viscosity of each particle in the type
            :type nu: float
            :param c: The artificial speed of sound of each particle in the type
            :type c: float
            :param fixed: Are the particles in this type immobile
            :type fixed: bool

            :rtype: int
        """

        if self.domain is None:
            raise Exception("SpatialPy models must have a domain before types can be attached");
        if type_id not in self.listOfTypeIDs:
            # index is the "particle type", value is the "type ID"
            self.listOfTypeIDs.append(type_id)
        # apply the type to all points, set type for any points that match
        count = 0
        on_boundary = self.domain.find_boundary_points()
        for v_ndx in range(self.domain.get_num_voxels()):
            if geometry_ivar.inside( self.domain.coordinates()[v_ndx,:], on_boundary[v_ndx]):
                self.domain.type[v_ndx] = type_id
                if vol is not None:
                    self.domain.vol[v_ndx] = vol
                if mass is not None:
                    self.domain.mass[v_ndx] = mass
                if rho is not None:
                    self.domain.rho[v_ndx] = rho
                if nu is not None:
                    self.domain.nu[v_ndx] = nu
                if c is not None:
                    self.domain.c[v_ndx] = c
                self.domain.fixed[v_ndx] = fixed
                count +=1
        if count == 0:
            warnings.warn("Type with type_id={0} has zero particles in it".format(type_id))
        return count

    def restrict(self, species, listOfTypes):
        """ Set the diffusion coefficient to zero for 'species' in all types not in
            'listOfTypes'. This effectively restricts the movement of 'species' to
            the types specified in 'listOfTypes'.

            :param species: Target species to restrict
            :type species: spatialpy.Model.Species
            :param listOfTypes: a list, each object in the list should be a 'type_id'
            :type listOfTypes: list(int)
        """
        #x = Species()
        #if not isinstance(species, Species):
        #if str(type(species)) != 'Species':
        #    raise ModelError("First argument to restrict() must be a Species object, not {0}".format(str(type(species))))
        if not isinstance(listOfTypes,list):
            self.listOfDiffusionRestrictions[species] = [listOfTypes]
        else:
            self.listOfDiffusionRestrictions[species] = listOfTypes

    def add_domain(self, domain):
        '''
        Add a spatial domain to the model

        :param domain: The Domain object to be added to the model
        :type domain: spatialpy.Domain.Domain
        '''
        if not isinstance(domain,Domain) and type(domain).__name__ != 'Domain':
            raise ModelError("Unexpected parameter for add_domain. Parameter must be a Domain.")

        self.domain = domain
        self.listOfTypeIDs = list(set(domain.type))

    def add_data_function(self, data_function):
        """ Add a scalar spatial function to the simulation. This is useful if you have a
            spatially varying in put to your model. Argument is a instances of subclass of the
            spatialpy.DataFunction class. It must implement a function 'map(x)' which takes a
            the spatial positon 'x' as an array, and it returns a float value.

            :param data_function: Data function to be added.
            :type data_function: spatialpy.DataFunction
        """

        if isinstance(data_function, list):
            for S in data_function:
                self.add_data_function(S)
        elif isinstance(data_function, DataFunction) or type(data_function).__name__ == 'DataFunction':
            problem = self.__problem_with_name(data_function.name)
            if problem is not None:
                raise problem
            self.listOfDataFunctions.append(data_function)
        else:
            raise ModelError("Unexpected parameter for add_data_function. Parameter must be DataFunction or list of DataFunctions.")
        return data_function

    def add_initial_condition(self, ic):
        """ Add an initial condition object to the initialization of the model.

                :param ic: Initial condition to be added
                :type ic: spatialpy.InitialCondition.InitialCondition

        """
        self.listOfInitialConditions.append(ic)

    def add_boundary_condition(self, bc):
        """ Add an BoundaryCondition object to the model.

            :param bc: Boundary condition to be added
            :type bc: spatialpy.BoundaryCondition.BoundaryCondition

        """
        bc.model = self
        self.listOfBoundaryConditions.append(bc)

    def update_namespace(self):
        """ Create a dict with flattened parameter and species objects. """

        for param in self.listOfParameters:
            self.namespace[param]=self.listOfParameters[param].value

    def sanitized_species_names(self):
        """
        Generate a dictionary mapping user chosen species names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user species names to their internal SpatialPy notation.
        """
        species_name_mapping = OrderedDict([])
        for i, name in enumerate(self.listOfSpecies.keys()):
            species_name_mapping[name] = 'x[{}]'.format(i)
        return species_name_mapping

    def sanitized_data_function_names(self):
        """
        Generate a dictionary mapping user chosen data function names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user data function names to their internal SpatialPy notation.
        """
        data_fn_name_mapping = OrderedDict([])
        for i, data_fn in enumerate(self.listOfDataFunctions):
            data_fn_name_mapping[data_fn.name] = 'data_fn[{}]'.format(i)
        return data_fn_name_mapping

    def sanitized_parameter_names(self):
        """
        Generate a dictionary mapping user chosen parameter names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user parameter names to their internal SpatialPy notation.
        """
        parameter_name_mapping = OrderedDict()
        for i, name in enumerate(self.listOfParameters.keys()):
            if name not in parameter_name_mapping:
                parameter_name_mapping[name] = 'P{}'.format(i)
        return parameter_name_mapping

    def get_species(self, sname):
        """ Returns target species from model as object.

                :param sname: name of species to be returned
                :type sname: str

                :rtype: spatialpy.Model.Species

        """
        return self.listOfSpecies[sname]

    def get_num_species(self):
        """ Returns total number of species contained in the model.

                :rtype: int

        """
        return len(self.listOfSpecies)

    def get_all_species(self):
        """ Returns a dictionary of all species in the model using names as keys.

                :rtype: dict

        """
        return self.listOfSpecies

    def add_species(self, obj):
        """
        Adds a species, or list of species to the model. Will return the added object upon success.

        :param obj: The species or list of species to be added to the model object.
        :type obj: spatialpy.Model.Species | list(spatialpy.Model.Species

        :rtype: spatialpy.Model.Species | list(spatialpy.Model.Species
        """


        if isinstance(obj, list):
            for S in obj:
                self.add_species(S)
        elif isinstance(obj, Species) or type(obj).__name__ == 'Species':
            problem = self.__problem_with_name(obj.name)
            if problem is not None:
                raise problem
            self.species_map[obj] = len(self.listOfSpecies)
            self.listOfSpecies[obj.name] = obj
        else:
            raise ModelError("Unexpected parameter for add_species. Parameter must be Species or list of Species.")
        return obj


    def delete_species(self, obj):
        """ Remove a Species from model.listOfSpecies. 

                :param obj: Species object to be removed
                :type obj: spatialpy.Model.Species

        """
        self.listOfSpecies.pop(obj)

    def delete_all_species(self):
        """ Remove all species from model.listOfSpecies.
        """

        self.listOfSpecies.clear()

    def get_parameter(self,pname):
        """ Remove a Parameter from model.listOfParameters. 

                :param pname: Name of parameter to be removed
                :type pname: spatialpy.Model.Parameter

        """
        try:
            return self.listOfParameters[pname]
        except:
            raise ModelError("No parameter named "+pname)

    def get_all_parameters(self):
        """ Return a dictionary of all model parameters, indexed by name.

            :rtype: dict

        """

        return self.listOfParameters

    def __problem_with_name(self, name):
        if name in Model.reserved_names:
            return ModelError('Name "{}" is unavailable. It is reserved for internal SpatialPy use. Reserved Names: ({}).'.format(name, Model.reserved_names))
        if name in self.listOfSpecies:
            return ModelError('Name "{}" is unavailable. A species with that name exists.'.format(name))
        if name in self.listOfParameters:
            return ModelError('Name "{}" is unavailable. A parameter with that name exists.'.format(name))
        if name.isdigit():
            return ModelError('Name "{}" is unavailable. Names must not be numeric strings.'.format(name))
        for special_character in Model.special_characters:
            if special_character in name:
                return ModelError('Name "{}" is unavailable. Names must not contain special characters: {}.'.format(name, Model.special_characters))



    def add_parameter(self,params):
        """ Add Parameter(s) to model.listOfParameters. Input can be either a single
            Parameter object or a list of Parameter objects.

                :param params: Parameter object or list of Parameters to be added.
                :type params: spatialpy.Model.Parameter | list(spatialpy.Model.Parameter)
        """
        if isinstance(params,list):
            for p in params:
                self.add_parameter(p)
        else:
            if isinstance(params, Parameter) or  type(params).__name__ == 'Parameter':
                problem = self.__problem_with_name(params.name)
                if problem is not None:
                    raise problem
                self.update_namespace()
                params._evaluate(self.namespace)
                self.listOfParameters[params.name] = params
            else:
                raise ParameterError("Parameter '{0}' needs to be of type '{2}', it is of type '{1}'".format(params.name,str(params),str(type(Parameter))))
        return params

    def delete_parameter(self, obj):
        self.listOfParameters.pop(obj)

    def _resolve_parameters(self):
        """ Attempt to resolve all parameter expressions to scalar floating point values.
            Must be called prior to exporting the model.  """
        self.update_namespace()
        for param in self.listOfParameters:
            self.listOfParameters[param]._evaluate(self.namespace)

    def delete_all_parameters(self):
        """ Remove all parameters from model.listOfParameters
        """

        self.listOfParameters.clear()

    def add_reaction(self,reacs):
        """ Add Reaction(s) to the model. Input can be single instance, a list of instances
            or a dict with name, instance pairs. 

            :param reacs: Reaction or list of Reactions to be added.
            :type reacs: spatialpy.Model.Reaction | list(spatialpy.Model.Reaction)

        """
        if isinstance(reacs, list):
            for r in reacs:
                r.initialize(self)
                if r.name is None or r.name == "":
                    r.name = 'rxn' + str(uuid.uuid4()).replace('-', '_')
                self.listOfReactions[r.name] = r
        elif isinstance(reacs, Reaction) or type(reacs).__name__ == "Reaction":
                reacs.initialize(self)
                if reacs.name is None or reacs.name == "":
                    reacs.name = 'rxn' + str(uuid.uuid4()).replace('-', '_')
                self.listOfReactions[reacs.name] = reacs
        else:
            raise ModelError("add_reaction() takes a spatialpy.Reaction object or list of objects")

    def get_reaction(self, rname):
        """ Retrieve a reaction object from the model by name

                :param rname: name of Reaction to retrieve
                :type rname: str

                :rtype: spatialpy.Model.Reaction

        """
        return self.listOfReactions[rname]

    def get_num_reactions(self):
        """ Returns the number of reactions in this model.

                :rtype: int
        """
        return len(self.listOfReactions)

    def get_all_reactions(self):
        """ Returns a dictionary of all model reactions using names as keys.

            :rtype: dict
        """

        return self.listOfReactions

    def delete_reaction(self, obj):
        """ Remove reaction from model.listOfReactions

            :param obj: Reaction to be removed.
            :type obj: spatialpy.Model.Reaction

        """
        self.listOfReactions.pop(obj)

    def delete_all_reactions(self):
        """ Remove all reactions from model.listOfReactions
        """
        self.listOfReactions.clear()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return (self.listOfParameters == other.listOfParameters and \
            self.listOfSpecies == other.listOfSpecies and \
            self.listOfReactions == other.listOfReactions and \
            self.name == other.name)

    def _create_stoichiometric_matrix(self):
        """ Generate a stoichiometric matrix in sparse CSC format. """

        if self.get_num_reactions() > 0:
            ND = numpy.zeros((self.get_num_species(), self.get_num_reactions()))
            for i, r in enumerate(self.listOfReactions):
                R = self.listOfReactions[r]
                reactants = R.reactants
                products  = R.products

                for s in reactants:
                    ND[self.species_map[s], i] -= reactants[s]
                for s in products:
                    ND[self.species_map[s], i] += products[s]

            N = scipy.sparse.csc_matrix(ND)
        else:
            N = numpy.zeros((self.get_num_species(), self.get_num_reactions()))

        return N


    def _create_dependency_graph(self):
        """ Construct the sparse dependency graph. """
        # We cannot safely generate a dependency graph (without attempting to analyze the
        # propensity string itself) if the model contains custom propensities.
        mass_action_model = True
        for name, reaction in self.listOfReactions.items():
            if not reaction.massaction:
                GF = numpy.ones((self.get_num_reactions(),
                    self.get_num_reactions() + self.get_num_species()))
                mass_action_model = False

        if mass_action_model:
            GF = numpy.zeros((self.get_num_reactions(),
                self.get_num_reactions() + self.get_num_species()))
            species_map = self.species_map

            involved_species = []
            reactants = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                temp2 = []
                for s in reaction.reactants:
                    temp.append(species_map[s])
                    temp2.append(species_map[s])
                for s in reaction.products:
                    temp.append(species_map[s])
                involved_species.append(temp)
                reactants.append(temp2)

            species_to_reactions = []
            for sname,species in self.listOfSpecies.items():
                temp = []
                for j, x in enumerate(reactants):
                    if species_map[species] in x:
                        temp.append(j)
                species_to_reactions.append(temp)

            reaction_to_reaction = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                for s in reaction.reactants:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+species_to_reactions[species_map[s]]

                for s in reaction.products:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+ species_to_reactions[species_map[s]]

                temp = list(set(temp))
                reaction_to_reaction.append(temp)

            # Populate G
            for j, spec in enumerate(species_to_reactions):
                for s in spec:
                    GF[s, j] = 1

            for i,reac in enumerate(reaction_to_reaction):
                for r in reac:
                    GF[r, self.get_num_species()+i] = 1


        try:
            G = scipy.sparse.csc_matrix(GF)
        except Exception as e:
            G = GF

        return G

    def _apply_initial_conditions(self):
        """ Initalize the u0 matrix (zeros) and then apply each initial condition"""
        # initalize
        ns = self.get_num_species()
        nv = self.domain.get_num_voxels()
        self.u0 = numpy.zeros((ns, nv))
        # apply initial condition functions
        for ic in self.listOfInitialConditions:
            ic.apply(self)
