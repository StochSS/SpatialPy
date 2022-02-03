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
import uuid

from spatialpy.core.Species import Species
from spatialpy.core.Parameter import Parameter
from spatialpy.core.spatialpyError import ReactionError

class Reaction():
    """
    Models a biochemical reaction. A reaction conatains dictionaries of species (reactants and products) \
    and parameters. The reaction's propensity function needs to be evaluable and result in a \
    non-negative scalar value in the namespace defined by the union of its Reactant, Product and \
    Parameter dictionaries.  For mass-action, zeroth, first \
    and second order reactions are supported, attempting to used higher orders will result in an error.

    :param name: String that the model is referenced by
    :type name: str
    :param parameters: A list of parameter instances
    :type parameters: list(spatialpy.Model.Parameter)
    :param propensity_function: String with the expression for the reaction's propensity
    :type propensity_function: str
    :param reactants: Dictionary of {species:stoichiometry} of reaction reactants
    :type reactants: dict
    :param products: Dictionary of {species:stoichiometry} of reaction products
    :type products: dict
    :param annotation: Description of the reaction (meta)
    :type annotation: str
    :param rate: if mass action, rate is a reference to a parameter instance.
    :type rate: spatialpy.model.Parameter
    """

    def __init__(self, name="", reactants={}, products={}, propensity_function=None, rate=None, annotation=None, restrict_to=None):

        if not isinstance(name, str) and name is not None:
            raise ReactionError("Reaction name must be of type str.")
        if name is None or name == "":
            self.name = f'rxn{uuid.uuid4()}'.replace('-', '_')
        else:
            self.name = name
        
        if not isinstance(reactants, dict):
            raise ReactionError("Reaction reactants must be of type dict.")
        self.reactants = reactants
        
        if not isinstance(products, dict):
            raise ReactionError("Reaction products must be of type dict.")
        self.products = products

        if not (isinstance(propensity_function, (str, int, float)) or propensity_function is None):
            raise ReactionError("Reaction propensity_function must be one of the following types: str, int, float, or None.")
        if isinstance(propensity_function, (int, float)):
            self.propensity_function = str(propensity_function)
        else:
            self.propensity_function = propensity_function

        if not (isinstance(rate, (str, int, float)) or rate is None or \
                                isinstance(rate, Parameter) or type(rate).__name__ == 'Parameter'):
            raise ReactionError("Reaction rate must be one of the following types: spatialpy.Parameter, str, int, float, or None.")
        if isinstance(rate, (int, float)):
            self.rate = str(rate)
        else:
            self.rate = rate

        if isinstance(restrict_to, int):
            self.restrict_to = [restrict_to]
        elif isinstance(restrict_to, list) or restrict_to is None:
            self.restrict_to = restrict_to
        else:
            errmsg = f"Reaction {self.name}: restrict_to must be an integer or list of integers."
            raise ReactionError(errmsg)

        self.annotation = annotation
        

    def __str__(self):
        print_string = f"{self.name}, Active in: {str(self.restrict_to)}"
        if len(self.reactants):
            print_string += f"\n\tReactants"
            for species, stoichiometry in self.reactants.items():
                name = species if isinstance(species, str) else species.name
                print_string += f"\n\t\t{name}: {stoichiometry}"
        if len(self.products):
            print_string += f"\n\tProducts"
            for species, stoichiometry in self.products.items():
                name = species if isinstance(species, str) else species.name
                print_string += f"\n\t\t{name}: {stoichiometry}"
        print_string += f"\n\tPropensity Function: {self.propensity_function}"
        return print_string


    def _create_mass_action(self):
        """
        Create a mass action propensity function given self.reactants and a single parameter value.

        We support zeroth, first and second order propensities only.
        There is no theoretical justification for higher order propensities.
        Users can still create such propensities if they really want to,
        but should then use a custom propensity.
        """
        total_stoch = 0
        for r in self.reactants:
            total_stoch += self.reactants[r]
        if total_stoch > 2:
            raise ReactionError(
                      """Reaction: A mass-action reaction cannot involve more than two of one species or one "
                         of two species (no more than 2 total reactants).
                         SpatialPy support zeroth, first and second order propensities only.
                         There is no theoretical justification for higher order propensities.
                         Users can still create such propensities using a 'custom propensity'.""")
        # Case EmptySet -> Y

        self.type = "mass-action"
        rtype = type(self.rate).__name__
        if rtype == 'Parameter':
            self.marate = self.rate.name
        elif rtype == 'int' or rtype == 'float':
            self.marate = str(self.rate)
        else:
            self.marate = self.rate

        propensity_function = self.marate
        ode_propensity_function = self.marate

        # There are only three ways to get 'total_stoch==2':
        for r in self.reactants:
            # Case 1: 2X -> Y
            if self.reactants[r] == 2:
                propensity_function = "0.5 * {0} * {1} * ({1} - 1) / vol".format(propensity_function, r.name)
                ode_propensity_function += " * {0} * {0}".format(r.name)
            else:
                # Case 3: X1, X2 -> Y;
                propensity_function += f" * {r.name}"
                ode_propensity_function += f" * {r.name}"

        # Set the volume dependency based on order.
        order = len(self.reactants)
        if order == 2:
            propensity_function += " / vol"
        elif order == 0:
            propensity_function += " * vol"

        self.propensity_function = propensity_function
        self.ode_propensity_function = ode_propensity_function

    def add_product(self, S, stoichiometry):
        """
        Add a product to this reaction

        :param S: Species object to be produced by the reaction
        :type S: spatialpy.Model.Species
        :param stoichiometry: Stoichiometry of this product.
        :type stoichiometry: int
        """
        if not (isinstance(S, (str, Species)) or type(S).__name__ == 'Species'):
            raise ReactionError("S must be of type string or SpatialPy.Species. ")
        if not isinstance(stoichiometry, int) or stoichiometry <= 0:
            raise ReactionError("Stoichiometry must be a positive integer.")
        name = S if isinstance(S, str) else S.name
        self.products[name] = stoichiometry

    def add_reactant(self, S, stoichiometry):
        """
        Add a reactant to this reaction
            
        :param S: reactant Species object
        :type S: spatialpy.Model.Species
        :param stoichiometry: Stoichiometry of this participant reactant
        :type stoichiometry: int
        """
        if not (isinstance(S, (str, Species)) or type(S).__name__ == 'Species'):
            raise ReactionError(f"S must be of type string or SpatialPy.Species. ")
        if not isinstance(stoichiometry, int) or stoichiometry <= 0:
            raise ReactionError(f"stoichiometry must be a positive integer.")
        name = S if isinstance(S, str) else S.name
        self.reactants[name] = stoichiometry

    def annotate(self, annotation):
        """
        Add an annotation to this reaction.

        :param annotation: Annotation note to be added to reaction
        :type annotation: str
        """
        if annotation is None:
            raise ReactionError("Annotation cannot be None.")
        self.annotation = annotation
    
    def initialize(self, model):
        """
        Defered object initialization, called by model.add_reaction().
        """
        self.ode_propensity_function = self.propensity_function

        if self.propensity_function is None:
            if self.rate is None:
                errmsg = f"Reaction {self.name}: You must either set the reaction to be mass-action or specifiy a propensity function."
                raise ReactionError(errmsg)
            self.massaction = True
        else:
            if self.rate is not None:
                errmsg = f"Reaction {self.name}: You cannot set the propensity type to mass-action and simultaneously set a propensity function."
                raise ReactionError(errmsg)
            # If they don't give us a propensity function and do give a rate, assume mass-action.
            self.massaction = False
            self.marate = None


        reactants = self.reactants
        self.reactants = {}
        if reactants is not None:
            for r in reactants:
                rtype = type(r).__name__
                if rtype=='Species':
                    if r.name not in model.listOfSpecies:
                        raise ReactionError(f"Could not find species '{r.name}' in model.")
                    self.reactants[r] = reactants[r]
                elif rtype=='str':
                    if r not in model.listOfSpecies:
                        raise ReactionError(f"Could not find species '{r}' in model.")
                    self.reactants[model.listOfSpecies[r]] = reactants[r]

        products = self.products
        self.products = {}
        if products is not None:
            for p in products:
                rtype = type(p).__name__
                if rtype=='Species':
                    if p.name not in model.listOfSpecies:
                        raise ReactionError(f"Could not find species '{p.name}' in model.")
                    self.products[p] = products[p]
                else:
                    if p not in model.listOfSpecies:
                        raise ReactionError(f"Could not find species '{p}' in model.")
                    self.products[model.listOfSpecies[p]] = products[p]

        if self.massaction:
            self._create_mass_action()
        else:
            self.type = "customized"


    