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
import unittest

import spatialpy
from spatialpy import Reaction
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
        rate_reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products, rate=test_rate)
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, test_reactants)
        self.assertEqual(reaction.products, test_products)
        self.assertEqual(reaction.rate, test_rate)


    def test_constructor__custom_propensity(self):
        """ Test the Reaction constructor for a custom propensity reaction. """
        test_reactants = {"A": 1, "B": 1}
        test_products = {"C": 1}
        test_propensity = "k1 * A * B"
        prop_reaction = Reaction(name="test_reaction", reactants=test_reactants, products=test_products,
                                 propensity_function=test_propensity)
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, test_reactants)
        self.assertEqual(reaction.products, test_products)
        self.assertEqual(reaction.propensity_function, test_propensity)
