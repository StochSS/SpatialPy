#This module defines a model that simulates a discrete, stoachastic, mixed biochemical reaction network in python.
    

from __future__ import division # is this still necessary?
import uuid
from collections import OrderedDict
from spatialpy.Solver import Solver
import numpy
import scipy


class Model():
    """ Representation of a spatial biochemical model. """
    reserved_names = ['vol']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']

    
    def __init__(self, name="", volume=1.0):
        """ Create an empty SpatialPy model. """
        
        # The name that the model is referenced by (should be a String)
        self.name = name
        
        # Dictionaries with Species, Reactions and Parameter objects.
        # Species,Reactio and Paramter names are used as keys.
        self.listOfParameters = OrderedDict()
        self.listOfSpecies    = OrderedDict()
        self.listOfReactions  = OrderedDict()
        
        # A well mixed model has an optional volume parameter
        self.volume = volume;
        
        # Dict that holds flattended parameters and species for
        # evaluation of expressions in the scope of the model.
        self.namespace = OrderedDict([])

        ######################
        self.mesh = None
        self.listOfSubdomainIDs = [1] # starts with subdomain '1'
        self.listOfDiffusionRestrictions = {}
        self.timestep_size = None
        self.num_timesteps = None
        self.listOfDataFunctions = []
        self.listOfInitialConditions = []
        self.species_map = {}
        self.tspan = None
        self.staticDomain = True;


    def run(self, number_of_trajectories=1, solver=None, seed=None, report_level=0):
        """ Simulate the model.
        Args:
            solver: A str or class type that is a subclass of SpatialPy.Solver.  Default: NSM solver.
            number_of_trajectories: How many trajectories should be run.
            seed: An int, the random seed given to the solver.
            report_level: An int, Level of output from the solver: 0, 1, or 2. Default: 0.
        Returns:
            A SpatialPY.Result object with the results of the simulation.
        """
        if solver is not None:
            if ((isinstance(solver, type)
                    and issubclass(solver, Solver))) or issubclass(type(solver), Solver):
                sol = solver(self, report_level=report_level)
        else:
            from spatialpy.nsmsolver import NSMSolver
            sol = NSMSolver(self, report_level=report_level)

        return sol.run(number_of_trajectories=number_of_trajectories, seed=seed)


    def set_timesteps(self, step_size, num_steps):
        """" Set the simlation time span parameters
        Args:
            step_size: float, size of each timestep in seconds
            num_steps: int, total number of steps to take

        Note: the number of output times will be num_steps+1 as the first
              output will be at time zero.
        """
        #TODO, add checking
        self.timestep_size = step_size
        self.num_timesteps = num_steps

    def timespan(self, time_span):
        """
        Set the time span of simulation. The SSA-SDPD engine does not support 
        non-uniform timespans.

        tspan : numpy ndarray
            Evenly-spaced list of times at which to sample the species
            populations during the simulation.
        """

        self.tspan = time_span

        items_diff = numpy.diff(time_span)
        items = map(lambda x: round(x, 10), items_diff)
        isuniform = (len(set(items)) == 1)

        if isuniform:
            self.timestep_size = items_diff[0]
            self.num_timesteps = len(items_diff)
        else:
            raise InvalidModelError("Only uniform timespans are supported")



    def add_subdomain(self, subdomain, domain_id, mass=1.0):
        """ Add a subdomain definition to the model.  By default, all regions are set to
        subdomain 1.
        Args:
            subdomain: an instance of a 'spatialpy.SubDomain' subclass.  The 'inside()' method
                       of this object will be used to assign domain_id to points.
            domain_id: the identifier for this subdomain (usually an int).
        Return values:
            Number of mesh points that were tagged with this domain_id
        """
        if self.mesh is None:
            raise Exception("SpatialPy models must have a mesh before subdomains can be attached");
        if domain_id not in self.listOfSubdomainIDs:
            # index is the "particle type", value is the "subdomain ID"
            self.listOfSubdomainIDs.append(domain_id)
        # apply the subdomain to all points, set sd for any points that match
        count =0
        on_boundary = self.mesh.find_boundary_points()
        for v_ndx in range(self.mesh.get_num_voxels()):
            if subdomain.inside( self.mesh.coordinates()[v_ndx,:], on_boundary[v_ndx]):
                self.mesh.sd[v_ndx] = domain_id
                self.mesh.mass[v_ndx] = mass
                count +=1
        return count

    def restrict(self, species, listOfSubDomains):
        """ Set the diffusion coefficient to zero for 'species' in all subdomains not in
            'listOfSubDomains'. This effectively restricts the movement of 'species' to
            the subdomains specified in 'listOfSubDomains'.
        Args:
            species: an instance of a 'spatialpy.Species'.
            listOfSubdomains: a list, each object in the list should be a 'domain_id'
        """
        #x = Species()
        #if not isinstance(species, Species):
        #if str(type(species)) != 'Species':
        #    raise ModelError("First argument to restrict() must be a Species object, not {0}".format(str(type(species))))
        if not isinstance(listOfSubDomains,list):
            self.listOfDiffusionRestrictions[species] = [listOfSubDomains]
        else:
            self.listOfDiffusionRestrictions[species] = listOfSubDomains

    def add_data_function(self, data_function):
        """ Add a scalar spatial function to the simulation.  This is useful if you have a 
            spatially varying in put to your model.  Argument is a instances of subclass of the 
            spatialpy.DataFunction class. It must implement a function 'map(x)' which takes a 
            the spatial positon 'x' as an array, and it returns a float value. 
        """
        #TODO validate input
        self.listOfDataFunctions.append(data_function)

    def add_initial_condition(self, ic):
        """ Add an initial condition object to the initialization of the model."""
        #TODO: validate model
        self.listOfInitialConditions.append(ic)
        

    def update_namespace(self):
        """ Create a dict with flattened parameter and species objects. """
        
        for param in self.listOfParameters:
            self.namespace[param]=self.listOfParameters[param].value
        # Dictionary of expressions that can be evaluated in the scope of this model.
        self.expressions = {}
    
    def get_species(self, sname):
        return self.listOfSpecies[sname]
    
    def get_num_species(self):
        return len(self.listOfSpecies)
    
    def get_all_species(self):
        return self.listOfSpecies

    def add_species(self, obj):
        """
        Adds a species, or list of species to the model.

        Attributes
        ----------
        obj : Species, or list of Species
            The species or list of species to be added to the model object.
        """


        if isinstance(obj, list):
            for S in obj:
                self.add_species(S)
        elif type(obj).__name__ == 'Species':
            problem = self.problem_with_name(obj.name)
            if problem is not None:
                raise problem
            self.species_map[obj] = len(self.listOfSpecies)
            self.listOfSpecies[obj.name] = obj
        else:
            raise ModelError("Unexpected parameter for add_species. Parameter must be Species or list of Species.")
        return obj


    def delete_species(self, obj):
        """ Remove a Species from model.listOfSpecies. """
        self.listOfSpecies.pop(obj)        
        
        
    def delete_all_species(self):
        self.listOfSpecies.clear()

    def get_parameter(self,pname):
        try:
            return self.listOfParameters[pname]
        except:
            raise ModelError("No parameter named "+pname)

    def get_all_parameters(self):
        return self.listOfParameters

    def problem_with_name(self, name):
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
        """
        if isinstance(params,list): 
            for p in params:
                self.add_parameter(p)
        else:
            #if isinstance(params, type(Parameter())):
            x = Parameter()
            if str(type(params)) == str(type(x)):
                problem = self.problem_with_name(params.name)
                if problem is not None:
                    raise problem
                # make sure that you don't overwrite an existing parameter??
                if params.name in self.listOfParameters.keys():
                    raise ParameterError("Parameter '{0}' has already been added to the model.".format(params.name))
                self.listOfParameters[params.name] = params
            else:
                #raise ParameterError("Could not resolve Parameter expression {} to a scalar value.".format(params))
                raise ParameterError("Parameter '{0}' needs to be of type '{2}', it is of type '{1}'".format(params.name,str(type(params)),str(type(x))))
        return params

    def delete_parameter(self, obj):
        self.listOfParameters.pop(obj)

    def set_parameter(self,pname,expression):
        """ Set the expression of an existing paramter. """
        p = self.listOfParameters[pname]
        p.expression = expression
        p.evaluate()
        
    def resolve_parameters(self):
        """ Attempt to resolve all parameter expressions to scalar floating point values. 
            Must be called prior to exporting the model.  """
        self.update_namespace()
        for param in self.listOfParameters:
            try:
                self.listOfParameters[param].evaluate(self.namespace)
            except:
                raise ParameterError("Could not resolve Parameter expression "+param + "to a scalar value.")
    
    def delete_all_parameters(self):
        self.listOfParameters.clear()

    def add_reaction(self,reacs):
        """ Add Reaction(s) to the model. Input can be single instance, a list of instances
            or a dict with name,instance pairs. """
        if isinstance(reacs, list): 
            for r in reacs:
                if r.name is None or r.name == "":
                    r.name = 'rxn' + str(uuid.uuid4()).replace('-', '_')
                self.listOfReactions[r.name] = r
        elif type(reacs).__name__ == "Reaction":
                if reacs.name is None or reacs.name == "":
                    reacs.name = 'rxn' + str(uuid.uuid4()).replace('-', '_')
                self.listOfReactions[reacs.name] = reacs
        else:
            raise ModelError("add_reaction() takes a spatialpy.Reaction object or list of objects")

    def get_reaction(self, rname):
        return self.listOfReactions[rname]

    def get_num_reactions(self):
        return len(self.listOfReactions)

    def get_all_reactions(self):
        return self.listOfReactions
    
    def delete_reaction(self, obj):
        self.listOfReactions.pop(obj)
        
    def delete_all_reactions(self):
        self.listOfReactions.clear()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return (self.listOfParameters == other.listOfParameters and \
            self.listOfSpecies == other.listOfSpecies and \
            self.listOfReactions == other.listOfReactions and \
            self.name == other.name)

    def create_stoichiometric_matrix(self):
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


    def create_dependency_graph(self):
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

    def apply_initial_conditions(self):
        """ Initalize the u0 matrix (zeros) and then apply each initial condition"""
        # initalize
        ns = self.get_num_species()
        nv = self.mesh.get_num_voxels()
        self.u0 = numpy.zeros((ns, nv))
        # apply initial condition functions
        for ic in self.listOfInitialConditions:
            ic.apply(self)
            




class Species():
    """ Model of a biochemical species. """
    
    def __init__(self,name=None,diffusion_constant=None,diffusion_coefficient=None,D=None):
        # A species has a name (string) and an initial value (positive integer)
        if name is None: 
            raise ModelError("Species must have a name")
        else:
            self.name = name
        if diffusion_constant is not None:
            self.diffusion_constant=diffusion_constant
        elif  diffusion_coefficient is not None:
            self.diffusion_constant=diffusion_coefficient
        elif D is not None:
            self.diffusion_constant=D
        else:
            raise ModelError("Species must have a diffusion_constant")


    def __str__(self):
        return self.name

class Parameter():
    """ 
        Model of a rate paramter. 
        A parameter can be given as a String expression (function) or directly as a scalar value.
        If given a String expression, it should be evaluable in the namespace of a parent Model.
        
    """
    
    def __init__(self,name="",expression=None,value=None):

        self.name = name        
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.expression = expression
        if expression != None:
            self.expression = str(expression)
        
        self.value = value
            
        # self.value is allowed to be None, but not self.expression. self.value
        # might not be evaluable in the namespace of this parameter, but defined
        # in the context of a model or reaction.
        if self.expression == None:
            #raise TypeError
            self.value = 0
    
        if self.value is None:
            self.evaluate()
    
    def evaluate(self,namespace={}):
        """ Evaluate the expression and return the (scalar) value """
        try:
            self.value = (float(eval(self.expression, namespace)))
        except:
            self.value = None
            
    def set_expression(self,expression):
        self.expression = expression
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        if expression is not None:
            self.expression = str(expression)
                    
        if self.expression is None:
            raise TypeError
    
        self.evaluate()
        
    def __str__(self):
        return str(self.value)
        

class Reaction():
    """ 
        Models a biochemical reaction. A reaction conatains dictinaries of species (reactants and products) and parameters.
        The reaction's propensity function needs to be evaluable and result in a non-negative scalar value
        in the namespace defined by the union of its Reactant, Product and Parameter dictionaries.
        
    """

    def __init__(self, name = "", reactants = {}, products = {}, propensity_function=None, massaction=None, rate=None, annotation=None,restrict_to=None):
        """ 
            Initializes the reaction using short-hand notation. 
            
            Input: 
                name:                       string that the model is referenced by
                parameters:                 a list of parameter instances
                propensity_function:        String with the expression for the reaction's propensity
                reactants:                  List of (species,stoichiometry) tuples
                product:                    List of (species,stoichiometry) tuples
                annotation:                 Description of the reaction (meta)
            
                massaction True,{False}     is the reaction of mass action type or not?
                rate                        if mass action, rate is a reference to a parameter instance.
            
            If massaction is set to true, propensity_function is not a valid argument. Instead, the 
            propensity function is constructed automatically. For mass-action, zeroth, first and second
            order reactions are supported, appemting to used higher orders will result in an error.
            
            Raises: ReactionError
            
        """
            
        # Metadata
        self.name = name
        self.annotation = ""
        
        self.massaction = massaction

        self.propensity_function = propensity_function

        if self.propensity_function is None: 
            if rate is None:
                errmsg = "Reaction "+self.name +": You must either set the reaction to be mass-action or specifiy a propensity function."
                raise ReactionError(errmsg)
            else:
                # If they don't give us a propensity function and do give a rate, assume mass-action.
                self.massaction = True
        else:
            if rate is not None:
                errmsg = "Reaction "+self.name +": You cannot set the propensity type to mass-action and simultaneously set a propensity function."
                raise ReactionError(errmsg)
            else:
                self.massaction = False
        
        self.reactants = {}
        if reactants is not None:
            for r in reactants:
                rtype = type(r).__name__
                if rtype=='instance':
                    self.reactants[r.name] = reactants[r]
                else:
                    self.reactants[r]=reactants[r]
    
        self.products = {}
        if products is not None:
            for p in products:
                rtype = type(p).__name__
                if rtype=='instance':
                    self.products[p.name] = products[p]
                else:
                    self.products[p]=products[p]

        if self.massaction:
            self.type = "mass-action"
            if rate is None:
                raise ReactionError("Reaction "+self.name +": A mass-action propensity has to have a rate.")
            self.marate = rate
            self.create_mass_action()
        else:
            self.type = "customized"

        self.restrict_to = restrict_to
                
    def create_mass_action(self):
        """ Create a mass action propensity function given self.reactants and a single parameter value.
        """
        # We support zeroth, first and second order propensities only.
        # There is no theoretical justification for higher order propensities.
        # Users can still create such propensities if they really want to,
        # but should then use a custom propensity.
        total_stoch = 0
        for r in self.reactants:
            total_stoch += self.reactants[r]
        if total_stoch > 2:
            raise ReactionError("Reaction: A mass-action reaction cannot involve more than two of one species or one "
                                "of two species.")
        # Case EmptySet -> Y

        propensity_function = self.marate.name

        # There are only three ways to get 'total_stoch==2':
        for r in self.reactants:
            # Case 1: 2X -> Y
            if self.reactants[r] == 2:
                propensity_function = ("0.5*" + propensity_function +
                                       "*" + str(r) + "*(" + str(r) + "-1)/vol")
            else:
                # Case 3: X1, X2 -> Y;
                propensity_function += "*" + str(r)

        # Set the volume dependency based on order.
        order = len(self.reactants)
        if order == 2:
            propensity_function += "/vol"
        elif order == 0:
            propensity_function += "*vol"

        self.propensity_function = propensity_function
            
    def set_type(self,type):
        if type not in {'mass-action','customized'}:
            raise ReactionError("Invalid reaction type.")
        self.type = type
    
    def add_reactant(self,S,stoichiometry):
        if stoichiometry <= 0:
            raise ReactionError("Reaction "+self.name+"Stoichiometry must be a positive integer.")
        self.reactants[S.name]=stoichiometry

    def add_product(self,S,stoichiometry):
        self.products[S.name]=stoichiometry

    def annotate(self,annotation):
        self.annotation = annotation


# Module exceptions
class ModelError(Exception):
    pass

class SpeciesError(ModelError):
    pass

class ReactionError(ModelError):
    pass

class ParameterError(ModelError):
    pass
