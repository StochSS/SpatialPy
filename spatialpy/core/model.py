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

#This module defines a model that simulates a discrete, stoachastic, mixed biochemical reaction network in python.

import math
from typing import Set, Type
from collections import OrderedDict

import numpy
import scipy

from spatialpy.core.domain import Domain
from spatialpy.core.species import Species
from spatialpy.core.initialcondition import (
    InitialCondition,
    PlaceInitialCondition,
    ScatterInitialCondition,
    UniformInitialCondition
)
from spatialpy.core.parameter import Parameter
from spatialpy.core.reaction import Reaction
from spatialpy.core.boundarycondition import BoundaryCondition
from spatialpy.core.datafunction import DataFunction
from spatialpy.core.timespan import TimeSpan
from spatialpy.solvers.build_expression import BuildExpression
from spatialpy.core.spatialpyerror import ModelError

def export_StochSS(spatialpy_model, filename=None, return_stochss_model=False):
    """
    SpatialPy model to StochSS converter

    :param spatialpy_model: SpatialPy model to be converted to StochSS
    :type spatialpy_model: spatialpy.core.model.Model

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

        # Dictionaries with model element objects.
        # Model element names are used as keys, and values are
        # sanitized versions of the names/formulas.
        # These dictionaries contain sanitized values and are for
        # Internal use only
        self._listOfParameters = OrderedDict()
        self._listOfSpecies = OrderedDict()
        self._listOfReactions = OrderedDict()
        
        ######################
        # Dict that holds flattended parameters and species for
        # evaluation of expressions in the scope of the model.
        self.namespace = OrderedDict([])
        self.species_map = {}

        ######################
        self.domain = None
        self.listOfDiffusionRestrictions = OrderedDict([])
        self.listOfDataFunctions = OrderedDict([])
        self.listOfInitialConditions = []
        self.listOfBoundaryConditions = []

        ######################
        self.staticDomain = True
        self.enable_rdme = True
        self.enable_pde = True

        ######################
        self.tspan = None

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
        self._resolve_all_parameters()
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

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.get_element(key)
        if hasattr(self.__class__, "__missing__"):
            return self.__class__.__missing__(self, key)
        raise KeyError(f"{key} is an invalid key.")

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

    def __update_diffusion_restrictions(self):
        for species in self.listOfSpecies.values():
            if isinstance(species.restrict_to, list):
                self.listOfDiffusionRestrictions[species] = species.restrict_to
            elif isinstance(species.restrict_to, str):
                self.listOfDiffusionRestrictions[species] = [species.restrict_to]

    def _ipython_display_(self, use_matplotlib=False):
        if self.domain is None:
            print(self)
        else:
            self.domain.plot_types(width="auto", height="auto", use_matplotlib=use_matplotlib)

    def _problem_with_name(self, name):
        names = Model.reserved_names
        if name in Model.reserved_names:
            raise ModelError(
                f'Name "{name}" is unavailable. It is reserved for internal GillesPy use. Reserved Names: ({names}).'
            )
        if name in self.listOfSpecies:
            raise ModelError(f'Name "{name}" is unavailable. A species with that name exists.')
        if name in self.listOfParameters:
            raise ModelError(f'Name "{name}" is unavailable. A parameter with that name exists.')
        if name in self.listOfReactions:
            raise ModelError(f'Name "{name}" is unavailable. A reaction with that name exists.')
        if name.isdigit():
            raise ModelError(f'Name "{name}" is unavailable. Names must not be numeric strings.')
        for special_character in Model.special_characters:
            if special_character in name:
                chars = Model.special_characters
                raise ModelError(
                    f'Name "{name}" is unavailable. Names must not contain special characters: {chars}.'
                )

    def _resolve_parameter(self, parameter):
        try:
            parameter.validate()
            self.update_namespace()
            parameter._evaluate(self.namespace) # pylint: disable=protected-access
        except ModelError as err:
            raise ModelError(
                f"Could not add/resolve parameter: {parameter.name}, Reason given: {err}"
            ) from err

    def _resolve_all_parameters(self):
        for _, parameter in self.listOfParameters.items():
            self._resolve_parameter(parameter)

    def _resolve_reaction(self, reaction):
        try:
            reaction.validate()

            # If the rate parameter exists in the reaction, confirm that it is a part of the model
            if reaction.marate is not None:
                name = reaction.marate if isinstance(reaction.marate, str) else reaction.marate.name
                if isinstance(reaction.marate, str) and not reaction.marate.replace(".", "", 1).isdigit():
                    reaction.marate = self.get_parameter(name)

            # Confirm that all species in reactants are part of the model
            for species in list(reaction.reactants.keys()):
                stoichiometry = reaction.reactants[species]
                name = species if isinstance(species, str) else species.name
                stoich_spec = self.get_species(name)
                if stoich_spec not in reaction.reactants:
                    reaction.reactants[stoich_spec] = stoichiometry
                    del reaction.reactants[species]

            # Confirm that all species in products are part of the model
            for species in list(reaction.products.keys()):
                stoichiometry = reaction.products[species]
                name = species if isinstance(species, str) else species.name
                stoich_spec = self.get_species(name)
                if stoich_spec not in reaction.products:
                    reaction.products[stoich_spec] = stoichiometry
                    del reaction.products[species]
        except ModelError as err:
            raise ModelError(f"Could not add/resolve reaction: {reaction.name}, Reason given: {err}") from err

    def _resolve_all_reactions(self):
        for _, reaction in self.listOfReactions.items():
            self._resolve_reaction(reaction)

    def update_namespace(self):
        """
        Create a dict with flattened parameter and species objects.
        """
        self.namespace = OrderedDict([])
        for param in self.listOfParameters:
            self.namespace[param] = self.listOfParameters[param].value

    def add(self, components):
        """
        Adds a component, or list of components to the model. If a list is provided, Species
        and Parameters are added before other components.  Lists may contain any combination
        of accepted types other than lists and do not need to be in any particular order.

        :param components: The component or list of components to be added the the model.
        :type components: Species, Parameters, Reactions, Domain, Data Function, \
                          Initial Conditions, Boundary Conditions, and TimeSpan or list

        :returns: The components that were added to the model.
        :rtype: Species, Parameters, Reactions, Domain, Data Function, \
                Initial Conditions, Boundary Conditions, TimeSpan, or list

        :raises ModelError: Component is invalid.
        """
        initialcondition_names = [
            PlaceInitialCondition.__name__,
            ScatterInitialCondition.__name__,
            UniformInitialCondition.__name__
        ]
        if isinstance(components, list):
            params = []
            others = []
            for component in components:
                if isinstance(component, Species) or type(component).__name__ in Species.__name__:
                    self.add_species(component)
                elif isinstance(component, Parameter) or type(component).__name__ in Parameter.__name__:
                    params.append(component)
                else:
                    others.append(component)

            for param in params:
                self.add_parameter(param)
            for component in others:
                self.add(component)
        elif isinstance(components, BoundaryCondition) or type(components).__name__ == BoundaryCondition.__name__:
            self.add_boundary_condition(components)
        elif isinstance(components, DataFunction) or type(components).__name__ == DataFunction.__name__:
            self.add_data_function(components)
        elif isinstance(components, Domain) or type(components).__name__ == Domain.__name__:
            self.add_domain(components)
        elif isinstance(components, InitialCondition) or type(components).__name__ in initialcondition_names:
            self.add_initial_condition(components)
        elif isinstance(components, Parameter) or type(components).__name__ == Parameter.__name__:
            self.add_parameter(components)
        elif isinstance(components, Reaction) or type(components).__name__ == Reaction.__name__:
            self.add_reaction(components)
        elif isinstance(components, Species) or type(components).__name__ == Species.__name__:
            self.add_species(components)
        elif isinstance(components, TimeSpan) or type(components).__name__ == TimeSpan.__name__:
            self.timespan(components)
        else:
            raise ModelError(f"Unsupported component: {type(components)} is not a valid component.")
        return components

    def get_element(self, name):
        """
        Get a model element specified by name.

        :param name: Name of the element to be returned.
        :type name: str

        :returns: The specified spatialpy.Model element.
        :rtype: Species, Parameters, Reactions, Domain, Data Function, or TimeSpan
        """
        if name in ("tspan", "timespan"):
            return self.tspan
        if name == "domain":
            return self.domain
        if name in self.listOfSpecies:
            return self.get_species(name)
        if name in self.listOfParameters:
            return self.get_parameter(name)
        if name in self.listOfReactions:
            return self.get_reaction(name)
        if name in self.listOfDataFunctions:
            return self.get_data_function(name)
        raise ModelError(f"{self.name} does not contain an element named {name}.")

    def add_domain(self, domain):
        """
        Add a spatial domain to the model

        :param domain: The Domain object to be added to the model
        :type domain: spatialpy.core.domain.Domain

        :raises ModelError: Invalid Domain object
        """
        if not isinstance(domain, Domain) and type(domain).__name__ != Domain.__name__:
            raise ModelError(
                "Unexpected parameter for add_domain. Parameter must be of type SpatialPy.Domain."
            )

        domain.compile_prep()
        self.domain = domain
        return domain

    def add_species(self, species):
        """
        Adds a species, or list of species to the model.

        :param species: The species or list of species to be added to the model object.
        :type species: spatialpy.core.species.Species | list(spatialpy.core.species.Species

        :returns: The species or list of species that were added to the model.
        :rtype: spatialpy.core.species.Species | list(spatialpy.core.species.Species)

        :raises ModelError: If an invalid species is provided or if Species.validate fails.
        """
        if isinstance(species, list):
            for spec in species:
                self.add_species(spec)
        elif isinstance(species, Species) or type(species).__name__ == "Species":
            try:
                species.validate()
                self._problem_with_name(species.name)
                self.species_map[species] = self.get_num_species()
                self.listOfSpecies[species.name] = species
                self._listOfSpecies[species.name] = f'S{len(self._listOfSpecies)}'
            except ModelError as err:
                errmsg = f"Could not add species: {species.name}, Reason given: {err}"
                raise ModelError(errmsg) from err
        else:
            errmsg = f"species must be of type Species or list of Species not {type(species)}"
            raise ModelError(errmsg)
        return species

    def delete_species(self, name):
        """
        Removes a species object by name.

        :param name: Name of the species object to be removed
        :type name: str

        :raises ModelError: If species is not part of the model.
        """
        try:
            self.listOfSpecies.pop(name)
            if name in self._listOfSpecies:
                self._listOfSpecies.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a species named {name}."
            ) from err

    def delete_all_species(self):
        """
        Removes all species from the model object.
        """
        self.listOfSpecies.clear()
        self._listOfSpecies.clear()

    def get_species(self, name):
        """
        Returns a species object by name.

        :param name: Name of the species object to be returned.
        :type name: str

        :returns: The specified species object.
        :rtype: spatialpy.core.species.Species

        :raises ModelError: If the species is not part of the model.
        """
        if name not in self.listOfSpecies:
            raise ModelError(f"Species {name} could not be found in the model.")
        return self.listOfSpecies[name]

    def get_all_species(self):
        """
        Returns a dictionary of all species in the model using names as keys.

        :returns: A dict of all species in the model, in the form: {name : species object}.
        :rtype: OrderedDict
        """
        return self.listOfSpecies

    def get_num_species(self):
        """
        Returns total number of different species contained in the model.

        :returns: Number of different species in the model.
        :rtype: int
        """
        return len(self.listOfSpecies)

    def sanitized_species_names(self):
        """
        Generate a dictionary mapping user chosen species names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user species names to their internal SpatialPy notation.
        :rtype: OrderedDict
        """
        species_name_mapping = OrderedDict([])
        for i, name in enumerate(self.listOfSpecies.keys()):
            species_name_mapping[name] = f'x[{i}]'
        return species_name_mapping

    def add_initial_condition(self, init_cond):
        """
        Add an initial condition object to the initialization of the model.

        :param init_cond: Initial condition to be added.
        :type init_cond: spatialpy.core.initialcondition.InitialCondition

        :returns: The initial condition or list of initial conditions that were added to the model.
        :rtype: spatialpy.core.initialcondition.InitialCondition | \
                list(spatialpy.core.initialcondition.InitialCondition)

        :raises ModelError: If an invalid initial condition is provided.
        """
        names = [
            PlaceInitialCondition.__name__,
            ScatterInitialCondition.__name__,
            UniformInitialCondition.__name__
        ]
        if isinstance(init_cond, list):
            for initial_condition in init_cond:
                self.add_initial_condition(initial_condition)
        elif isinstance(init_cond, InitialCondition) or type(init_cond).__name__ in names:
            self.listOfInitialConditions.append(init_cond)
        else:
            errmsg = f"init_cond must be of type InitialCondition or list of InitialCondition not {type(init_cond)}"
            raise ModelError(errmsg)
        return init_cond

    def delete_initial_condition(self, init_cond):
        """
        Removes an initial condition object from the model object.

        :param init_cond: initial condition object to be removed.
        :type init_cond: spatialpy.core.InitialCondition

        :raises ModelError: If the initial condition is not part of the model.
        """
        try:
            index = self.listOfInitialConditions.index(init_cond)
            self.listOfInitialConditions.pop(index)
        except ValueError as err:
            raise ModelError(
                f"{self.name} does not contain this initial condition."
            ) from err

    def delete_all_initial_conditions(self):
        """
        Removes all initial conditions from the model object.
        """
        self.listOfInitialConditions.clear()

    def get_all_initial_conditions(self):
        """
        Returns a list of all initial conditions in the model.

        :returns: A list of all initial conditions in the model.
        :rtype: list
        """
        return self.listOfInitialConditions

    def add_parameter(self, parameters):
        """
        Adds a parameter, or list of parameters to the model.

        :param parameters:  The parameter or list of parameters to be added to the model object.
        :type parameters: spatialpy.core.parameter.Parameter | list(spatialpy.core.parameter.Parameter)

        :returns: A parameter or list of Parameters that were added to the model.
        :rtype: spatialpy.core.parameter.Parameter | list(spatialpy.core.parameter.Parameter)

        :raises ModelError: If an invalid parameter is provided or if Parameter.validate fails.
        """
        if isinstance(parameters, list):
            for param in parameters:
                self.add_parameter(param)
        elif isinstance(parameters, Parameter) or type(parameters).__name__ == 'Parameter':
            self._problem_with_name(parameters.name)
            self._resolve_parameter(parameters)
            self.listOfParameters[parameters.name] = parameters
            self._listOfParameters[parameters.name] = f'P{len(self._listOfParameters)}'
        else:
            errmsg = f"parameters must be of type Parameter or list of Parameter not {type(parameters)}"
            raise ModelError(errmsg)
        return parameters

    def delete_parameter(self, name):
        """
        Removes a parameter object by name.

        :param name: Name of the parameter object to be removed.
        :type name: str
        """
        try:
            self.listOfParameters.pop(name)
            if name in self._listOfParameters:
                self._listOfParameters.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a parameter named {name}"
            ) from err

    def delete_all_parameters(self):
        """
        Removes all parameters from model object.
        """
        self.listOfParameters.clear()
        self._listOfParameters.clear()

    def get_parameter(self, name):
        """
        Returns a parameter object by name.

        :param name: Name of the parameter object to be returned
        :type name: str

        :returns: The specified parameter object.
        :rtype: Spatialpy.core.parameter.Parameter

        :raises ModelError: If the parameter is not part of the model.
        """
        if name not in self.listOfParameters:
            raise ModelError(f"Parameter {name} could not be found in the model.")
        return self.listOfParameters[name]

    def get_all_parameters(self):
        """
        Return a dictionary of all model parameters, indexed by name.

        :returns: A dict of all parameters in the model, in the form: {name : parameter object}
        :rtype: OrderedDict
        """
        return self.listOfParameters

    def sanitized_parameter_names(self):
        """
        Generate a dictionary mapping user chosen parameter names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user parameter names to their internal SpatialPy notation.
        :rtype: OrderedDict
        """
        parameter_name_mapping = OrderedDict()
        for i, name in enumerate(self.listOfParameters.keys()):
            if name not in parameter_name_mapping:
                parameter_name_mapping[name] = f'P{i}'
        return parameter_name_mapping

    def add_reaction(self, reactions):
        """
        Adds a reaction, or list of reactions to the model.

        :param reactions: The reaction or list of reactions to be added to the model object
        :type reactions: spatialpy.core.reaction.Reaction | list(spatialpy.core.reaction.Reaction)

        :returns: The reaction or list of reactions that were added to the model.
        :rtype: spatialpy.core.reaction.Reaction | list(spatialpy.core.reaction.Reaction)

        :raises ModelError: If an invalid reaction is provided or if Reaction.validate fails.
        """
        if isinstance(reactions, list):
            for reaction in reactions:
                self.add_reaction(reaction)
        elif isinstance(reactions, Reaction) or type(reactions).__name__ == "Reaction":
            self._problem_with_name(reactions.name)
            self._resolve_reaction(reactions)
            self.listOfReactions[reactions.name] = reactions
            # Build Sanitized reaction as well
            sanitized_reaction = reactions._create_sanitized_reaction(
                len(self.listOfReactions), self._listOfSpecies, self._listOfParameters
            )
            self._listOfReactions[reactions.name] = sanitized_reaction
        else:
            errmsg = f"reactions must be of type Reaction or list of Reaction not {type(reactions)}"
            raise ModelError(errmsg)
        return reactions

    def delete_reaction(self, name):
        """
        Removes a reaction object by name.

        :param name: Name of the reaction object to be removed.
        :type name: str
        """
        try:
            self.listOfReactions.pop(name)
            if name in self._listOfReactions:
                self._listOfReactions.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a reaction named {name}"
            ) from err

    def delete_all_reactions(self):
        """
        Removes all reactions from the model object.
        """
        self.listOfReactions.clear()
        self._listOfReactions.clear()

    def get_reaction(self, rname):
        """
        Returns a reaction object by name.

        :param name: Name of the reaction object to be returned
        :type name: str

        :returns: The specified reaction object.
        :rtype: spatialpy.core.reaction.Reaction

        :raises ModelError: If the reaction is not part of the model.
        """
        if name not in self.listOfReactions:
            raise ModelError(f"Reaction {name} could not be found in the model.")
        return self.listOfReactions[name]

    def get_all_reactions(self):
        """
        Returns a dictionary of all model reactions using names as keys.

        :returns: A dict of all reaction in the model, in the form: {name : reaction object}.
        :rtype: OrderedDict
        """
        return self.listOfReactions

    def get_num_reactions(self):
        """
        Returns the number of reactions in this model.

        :returns: The total number of different reactions in the model.
        :rtype: int
        """
        return len(self.listOfReactions)

    def add_boundary_condition(self, bound_cond):
        """
        Add an boundary condition object to the model.

        :param bound_cond: Boundary condition to be added
        :type bound_cond: spatialpy.core.boundarycondition.BoundaryCondition

        :returns: The boundary condition or list of boundary conditions that were added to the model.
        :rtype: spatialpy.core.boundarycondition.BoundaryCondition | \
                list(spatialpy.core.boundarycondition.BoundaryCondition)

        :raises ModelError: If an invalid boundary conidition is provided.
        """
        if isinstance(bound_cond, list):
            for boundary_condition in bound_cond:
                self.add_boundary_condition(boundary_condition)
        elif isinstance(bound_cond, BoundaryCondition) or type(bound_cond).__name__ in "BoundaryCondition":
            bound_cond.model = self
            self.listOfBoundaryConditions.append(bound_cond)
        else:
            errmsg = f"bound_cond must be of type BoundaryCondition or list of BoundaryCondition not {type(bound_cond)}"
            raise ModelError(errmsg)
        return bound_cond

    def delete_boundary_condition(self, bound_cond):
        """
        Removes an boundary condition object from the model object.

        :param bound_cond: boundary condition object to be removed.
        :type bound_cond: spatialpy.core.BoundaryCondition

        :raises ModelError: If the boundary condition is not part of the model.
        """
        try:
            index = self.listOfBoundaryConditions.index(bound_cond)
            self.listOfBoundaryConditions.pop(index)
        except ValueError as err:
            raise ModelError(
                f"{self.name} does not contain this boundary condition."
            ) from err

    def delete_all_boundary_conditions(self):
        """
        Removes all boundary conditions from the model object.
        """
        self.listOfBoundaryConditions.clear()

    def get_all_boundary_conditions(self):
        """
        Returns a list of all boundary conditions in the model.

        :returns: A list of all boundary conditions in the model.
        :rtype: list
        """
        return self.listOfBoundaryConditions

    def add_data_function(self, data_function):
        """
        Add a scalar spatial function to the simulation. This is useful if you have a
        spatially varying input to your model. Argument is a instances of subclass of the
        spatialpy.DataFunction class. It must implement a function 'map(point)' which takes a
        the spatial positon 'point' as an array, and it returns a float value.

        :param data_function: Data function to be added.
        :type data_function: spatialpy.DataFunction

        :returns: DataFunction object(s) added tothe model.
        :rtype: spatialpy.core.datafunction.DataFunction | list(spatialpy.core.datafunction.DataFunction)

        :raises ModelError: Invalid DataFunction
        """
        if isinstance(data_function, list):
            for data_fn in data_function:
                self.add_data_function(data_fn)
        elif isinstance(data_function, DataFunction) or type(data_function).__name__ == 'DataFunction':
            self._problem_with_name(data_function.name)
            self.listOfDataFunctions[data_function.name] = data_function
        else:
            errmsg = f"data_function must be of type DataFunction or list of DataFunction not {type(data_function)}"
            raise ModelError(errmsg)
        return data_function

    def delete_data_function(self, name):
        """
        Removes an data function object from the model object.

        :param name: data function object to be removed.
        :type name: spatialpy.core.DataFunction

        :raises ModelError: If the data function is not part of the model.
        """
        try:
            self.listOfDataFunctions.pop(name)
        except ValueError as err:
            raise ModelError(
                f"{self.name} does not contain a data function named {name}."
            ) from err

    def delete_all_data_functions(self):
        """
        Removes all data functions from the model object.
        """
        self.listOfDataFunctions.clear()

    def get_data_function(self, name):
        """
        Returns a data function object by name.

        :param name: Name of the data function object to be returned
        :type name: str

        :returns: The specified data function object.
        :rtype: spatialpy.core.datafunction.DataFunction

        :raises ModelError: If the data function is not part of the model.
        """
        if name not in self.listOfDataFunctions:
            raise ModelError(f"Data function {name} could not be found in the model.")
        return self.listOfDataFunctions[name]

    def get_all_data_functions(self):
        """
        Returns a dict of all data functions in the model.

        :returns: A dict of all data functions in the model.
        :rtype: OrderedDict
        """
        return self.listOfDataFunctions

    def sanitized_data_function_names(self):
        """
        Generate a dictionary mapping user chosen data function names to simplified formats which will be used
        later on by SpatialPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user data function names to their internal SpatialPy notation.
        :rtype: OrderedDict
        """
        data_fn_name_mapping = OrderedDict([])
        for i, name in enumerate(self.listOfDataFunctions.keys()):
            data_fn_name_mapping[name] = f'data_fn[{i}]'
        return data_fn_name_mapping

    def set_timesteps(self, output_interval, num_steps, timestep_size=None):
        """
        Set the simlation time span parameters.

        :param output_interval: size of each output timestep in seconds
        :type output_interval:  float

        :param num_steps: total number of steps to take. Note: the number of output times will be num_steps+1 \
        as the first output will be at time zero.
        :type num_steps: int

        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: float
        """
        self.tspan = TimeSpan.arange(
            output_interval, t=num_steps * output_interval, timestep_size=timestep_size
        )

    def timespan(self, time_span, timestep_size=None):
        """
        Set the time span of simulation. The SSA-SDPD engine does not support
        non-uniform timespans.

        :param tspan: Evenly-spaced list of times at which to sample the species populations during the simulation.
        :type tspan: numpy.ndarray

        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: float
        """
        if isinstance(time_span, TimeSpan) or type(time_span).__name__ == "TimeSpan":
            self.tspan = time_span
            if timestep_size is not None:
                self.tspan.timestep_size = timestep_size
                self.tspan.validate(coverage="all")
        else:
            self.tspan = TimeSpan(time_span, timestep_size=timestep_size)

    def compile_prep(self):
        """
        Make sure all paramters are evaluated to scalars, update the models diffusion restrictions,
        create the models expression utility, and generate the domain list of type ids in preperation
        of compiling the simulation files.

        :returns: The stoichiometric and dependency_graph
        :rtype: tuple

        :raises ModelError: Timestep size exceeds output frequency or Model is missing a domain
        """
        try:
            self.tspan.validate(coverage="all")
        except ModelError as err:
            raise ModelError(f"Failed to validate timespan. Reason given: {err}") from err

        if self.domain is None:
            raise ModelError("The model's domain is not set.  Use 'add_domain()'.")
        self.domain.compile_prep()
        
        self.__update_diffusion_restrictions()
        self.__apply_initial_conditions()
        self._resolve_all_parameters()
        self._resolve_all_reactions()

        sanitized_params = self.sanitized_parameter_names()
        for species in self.listOfSpecies.values():
            diff_coeff = species.diffusion_coefficient
            if isinstance(diff_coeff, str):
                if diff_coeff not in sanitized_params:
                    raise ModelError(f"Parameterm {diff_coeff} doesn't exist.")
                species.diffusion_coefficient = sanitized_params[diff_coeff]

        self.__get_expression_utility()
        stoich_matrix = self.__create_stoichiometric_matrix()
        dep_graph = self.__create_dependency_graph()

        return stoich_matrix, dep_graph

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
        :rtype: spatialpy.core.result.Result
        """
        from spatialpy.solvers.solver import Solver # pylint: disable=import-outside-toplevel

        sol = Solver(self, debug_level=debug_level)

        return sol.run(number_of_trajectories=number_of_trajectories, seed=seed, timeout=timeout,
                       number_of_threads=number_of_threads, debug=debug, profile=profile)
