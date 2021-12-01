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
import re
import unittest

import spatialpy
from spatialpy import Reaction, Parameter, Species
from spatialpy import ReactionError

class TestReaction(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for spatialpy.Reaction.
    ################################################################################################
    '''
    def test_constructor__mass_action(self):
        """ Test the Reaction constructor for a mass action reaction. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = "k1"
        reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate)
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, test_reactants)
        self.assertEqual(reaction.products, test_products)
        self.assertEqual(reaction.rate, test_rate)


    def test_constructor__custom_propensity(self):
        """ Test the Reaction constructor for a custom propensity reaction. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_propensity = "k1 * A * B"
        reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, propensity_function=test_propensity
        )
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, test_reactants)
        self.assertEqual(reaction.products, test_products)
        self.assertEqual(reaction.propensity_function, test_propensity)


    def test_constructor__no_name(self):
        """ Test the Reaction constructor with no name provided. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_propensity = "k1 * A * B"
        reaction = Reaction(
            reactants=test_reactants, products=test_products, propensity_function=test_propensity
        )
        self.assertIsNotNone(re.search("rxn.*", reaction.name))


    def test_constructor__name_is_none(self):
        """ Test the Reaction constructor with None as name. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_propensity = "k1 * A * B"
        reaction = Reaction(
            name=None, reactants=test_reactants, products=test_products, propensity_function=test_propensity
        )
        self.assertIsNotNone(re.search("rxn.*", reaction.name))


    def test_constructor__name_not_str_or_None(self):
        """ Test the Reaction constructor with non-str name. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = "k1"
        with self.assertRaises(ReactionError):
            reaction = Reaction(name=0, reactants=test_reactants, products=test_products, rate=test_rate)


    def test_constructor__reactants_not_dict(self):
        """ Test the Reaction constructor with non-dict reactants. """
        test_reactants = [["C", 1]]
        test_products = {"A": 1, "B": 1}
        test_rate = "k2"
        with self.assertRaises(ReactionError):
            reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate)

    def test_constructor__products_not_dict(self):
        """ Test the Reaction constructor with non-dict products. """
        test_reactants = {"A": 1, "B": 1}
        test_products = [["C", 1]]
        test_rate = "k1"
        with self.assertRaises(ReactionError):
            reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate)


    def test_constructor__propensity_function_not_accepted_type(self):
        """ Test the Reaction constructor with a propensity function that is not of the proper type. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_propensity = ["k1 * A * B"]
        with self.assertRaises(ReactionError):
            reaction = Reaction(
                name="test_reaction", reactants=test_reactants, products=test_products, propensity_function=test_propensity
            )


    def test_constructor__int_propensity_function(self):
        """ Test the Reaction constructor with an int propensity function. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_propensity = 20
        reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, propensity_function=test_propensity
        )
        self.assertIsInstance(reaction.propensity_function, str)
        self.assertEqual(reaction.propensity_function, "20")


    def test_constructor__float_propensity_function(self):
        """ Test the Reaction constructor with a float propensity function. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_propensity = 0.5
        reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, propensity_function=test_propensity
        )
        self.assertIsInstance(reaction.propensity_function, str)
        self.assertEqual(reaction.propensity_function, "0.5")


    def test_constructor__rate_not_accepted_type(self):
        """ Test the Reaction constructor with a rate that is not of the proper type. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = ["k1"]
        with self.assertRaises(ReactionError):
            reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate)


    def test_constructor__int_rate(self):
        """ Test the Reaction constructor with an int rate. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = 20
        reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate)
        self.assertIsInstance(reaction.rate, str)
        self.assertEqual(reaction.rate, "20")


    def test_constructor__float_rate(self):
        """ Test the Reaction constructor with a float rate. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = 0.5
        reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate)
        self.assertIsInstance(reaction.rate, str)
        self.assertEqual(reaction.rate, "0.5")


    def test_constructor__restrict_to_not_accepted_type(self):
        """ Test the Reaction constructor with a restrict_to that is not the proper type. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = "k1"
        test_restrict_to = "1"
        with self.assertRaises(ReactionError):
            reaction = Reaction(
                name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate, restrict_to=test_restrict_to
            )


    def test_constructor__int_restrict_to(self):
        """ Test the Reaction constructor with a int restrict_to. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = "k1"
        test_restrict_to = 1
        reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate, restrict_to=test_restrict_to
        )
        self.assertIsInstance(reaction.restrict_to, list)
        self.assertEqual(reaction.restrict_to, [1])


    def test___str__(self):
        """ Test Reaction.__str__ method. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_rate = "k1"
        test_restrict_to = 1
        reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate, restrict_to=test_restrict_to
        )
        self.assertIsInstance(str(reaction), str)


    def test__create_mass_action__total_stoch_3(self):
        """ Test Reaction._create_mass_action total stochiometry > 2. """
        test_reactants = {"A": 1, "B": 2}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction._create_mass_action()


    def test__create_mass_action__marate_type_as_string(self):
        """ Test Reaction._create_mass_action marate as string. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction._create_mass_action()
        self.assertEqual(test_reaction.propensity_function, "k1 * vol")
        self.assertEqual(test_reaction.ode_propensity_function, "k1")


    def test__create_mass_action__marate_type_as_int(self):
        """ Test Reaction._create_mass_action marate as int. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = 1
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction._create_mass_action()
        self.assertEqual(test_reaction.propensity_function, "1 * vol")
        self.assertEqual(test_reaction.ode_propensity_function, "1")


    def test__create_mass_action__marate_type_as_float(self):
        """ Test Reaction._create_mass_action marate as float. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = 0.5
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction._create_mass_action()
        self.assertEqual(test_reaction.propensity_function, "0.5 * vol")
        self.assertEqual(test_reaction.ode_propensity_function, "0.5")


    def test__create_mass_action__marate_type_as_parameter(self):
        """ Test Reaction._create_mass_action marate as parameter. """
        test_parameter = Parameter("k1", expression=0.1)
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = test_parameter
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction._create_mass_action()
        self.assertEqual(test_reaction.propensity_function, "k1 * vol")
        self.assertEqual(test_reaction.ode_propensity_function, "k1")


    def test__create_mass_action__X_to_Y(self):
        """ Test Reaction._create_mass_action X -> Y. """
        test_species_x = Species(name="X", diffusion_coefficient=0.1)
        test_reactants = {test_species_x: 1}
        test_products = {"Y": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction._create_mass_action()
        self.assertEqual(test_reaction.propensity_function, "k1 * X")
        self.assertEqual(test_reaction.ode_propensity_function, "k1 * X")


    def test__create_mass_action__X_plus_Y_to_Z(self):
        """ Test Reaction._create_mass_action X + Y -> Z. """
        test_species_x = Species(name="X", diffusion_coefficient=0.1)
        test_species_y = Species(name="Y", diffusion_coefficient=0.1)
        test_reactants = {test_species_x: 1, test_species_y: 1}
        test_products = {"Z": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction._create_mass_action()
        self.assertEqual(test_reaction.propensity_function, "k1 * X * Y / vol")
        self.assertEqual(test_reaction.ode_propensity_function, "k1 * X * Y")


    def test__create_mass_action__2X_to_Y(self):
        """ Test Reaction._create_mass_action 2X -> Y. """
        test_species_x = Species(name="X", diffusion_coefficient=0.1)
        test_reactants = {test_species_x: 2}
        test_products = {"Y": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction._create_mass_action()
        self.assertEqual(test_reaction.propensity_function, "0.5 * k1 * X * (X - 1) / vol")
        self.assertEqual(test_reaction.ode_propensity_function, "k1 * X * X")


    def test_add_product__species_object(self):
        """ Test Reaction.add_product species is SpatialPy.Species. """
        test_species = Species(name="X", diffusion_coefficient=0.1)
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction.add_product(test_species, 1)
        self.assertIn(test_species.name, test_reaction.products)
        self.assertEqual(test_reaction.products[test_species.name], 1)


    def test_add_product__species_string(self):
        """ Test Reaction.add_product species is string. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction.add_product("X", 1)
        self.assertIn("X", test_reaction.products)
        self.assertEqual(test_reaction.products["X"], 1)


    def test_add_product__species_not_accepted_type(self):
        """ Test Reaction.add_product species not string or SpatialPy.Species. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction.add_product(["X"], 1)


    def test_add_product__stochiometry_not_int(self):
        """ Test Reaction.add_product stochiometry not integer. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction.add_product("X", 0.1)


    def test_add_product__stochiometry_less_than_zero(self):
        """ Test Reaction.add_product stochiometry <= 0. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction.add_product("X", 0)


    def test_add_reactant__species_object(self):
        """ Test Reaction.add_reactant species is SpatialPy.Species. """
        test_species = Species(name="X", diffusion_coefficient=0.1)
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction.add_reactant(test_species, 1)
        self.assertIn(test_species.name, test_reaction.reactant)
        self.assertEqual(test_reaction.reactant[test_species.name], 1)


    def test_add_reactant__species_string(self):
        """ Test Reaction.add_reactant species is string. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_reaction.add_reactant("X", 1)
        self.assertIn("X", test_reaction.reactant)
        self.assertEqual(test_reaction.reactant["X"], 1)


    def test_add_reactant__species_not_accepted_type(self):
        """ Test Reaction.add_reactant species not string or SpatialPy.Species. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction.add_reactant(["X"], 1)


    def test_add_reactant__stochiometry_not_int(self):
        """ Test Reaction.add_reactant stochiometry not integer. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction.add_reactant("X", 0.1)


    def test_add_reactant__stochiometry_less_than_zero(self):
        """ Test Reaction.add_reactant stochiometry <= 0. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction.add_reactant("X", 0)


    def test_annotation(self):
        """ Test Reaction.annotation. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        test_annotation = "Testing the reaction annotation message."
        test_reaction.annotation(test_annotation)
        self.assertEqual(test_reaction.annotation, test_annotation)


    def test_annotation__annotation_is_None(self):
        """ Test Reaction.annotation annotation is None. """
        test_reactants = {}
        test_products = {"C": 1}
        test_rate = "k1"
        test_reaction = Reaction(
            name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate
        )
        with self.assertRaises(ReactionError):
            test_reaction.annotation(None)
