# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2023 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import unittest

import spatialpy
from spatialpy import Species, Parameter
from spatialpy import SpeciesError

class TestSpecies(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for spatialpy.Species.
    ################################################################################################
    '''
    def setUp(self):
        self.species = Species(name="test_species", diffusion_coefficient=0)

    def test_constructor(self):
        """ Test the Species constructor. """
        species = Species(name="test_species", diffusion_coefficient=0)
        self.assertEqual(species.name, "test_species")
        self.assertEqual(species.diffusion_coefficient, 0)

    def test_constructor__no_name(self):
        """ Test the Species constructor without name. """
        with self.assertRaises(SpeciesError):
            species = Species(diffusion_coefficient=0)

    def test_constructor__invalid_name(self):
        """ Test the Species constructor with an invalid name. """
        test_names = ["", None, 0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(SpeciesError):
                    species = Species(name=test_name, diffusion_coefficient=0)

    def test_constructor__no_diffusion_coefficient(self):
        """ Test the Species constructor without diffusion_coefficient. """
        with self.assertRaises(SpeciesError):
            species = Species(name="test_species")

    def test_constructor__negative_diffusion_coefficient(self):
        """ Test the Species constructor with negative diffusion_coefficient. """
        with self.assertRaises(SpeciesError):
            species = Species(name="test_species", diffusion_coefficient=-1)

    def test_constructor__parameter_diffusion_coefficient(self):
        """ Test the Species constructor with a parameter for diffusion_coefficient. """
        test_parameter = Parameter(name="test_parameter", expression=0.5)
        species = Species(name="test_species", diffusion_coefficient=test_parameter)
        self.assertIsInstance(species.diffusion_coefficient, Parameter)

    def test_constructor__invalid_diffusion_coefficient(self):
        """ Test the Species constructor with an invalid diffusion_coefficient. """
        test_dcs = [None, [0]]
        for test_dc in test_dcs:
            with self.subTest(initial_value=test_dc):
                with self.assertRaises(SpeciesError):
                    species = Species(name="test_species", diffusion_coefficient=test_dc)

    def test_constructor__invalid_restrict_to(self):
        """ Test the Species constructor with an invalid restrict_to. """
        test_rts = [0.5, [], (5, 2), {'x': 5}]
        for test_rt in test_rts:
            with self.subTest(restrict_to=test_rt):
                with self.assertRaises(SpeciesError):
                    species = Species(name="test_species", diffusion_coefficient="0", restrict_to=test_rt)

    def test_constructor__int_restrict_to(self):
        """ Test the Species constructor with a int restrict_to. """
        species = Species(name="test_species", diffusion_coefficient=0, restrict_to=1)
        self.assertIsInstance(species.restrict_to, list)
        self.assertEqual(species.restrict_to, ["type_1"])

    def test_constructor__str_restrict_to(self):
        """ Test the Species constructor with a string restrict_to. """
        species = Species(name="test_species", diffusion_coefficient=0, restrict_to="Walls")
        self.assertIsInstance(species.restrict_to, list)
        self.assertEqual(species.restrict_to, ["type_Walls"])

    def test___str___(self):
        """ Test Species.__str__ method. """
        self.assertIsInstance(str(self.species), str)

    def test_set_diffusion_coefficient(self):
        """ Test Species.set_diffusion_coefficient method. """
        self.species.set_diffusion_coefficient(diffusion_coefficient=1)
        self.assertEqual(self.species.diffusion_coefficient, 1)

    def test_set_diffusion_coefficient__invalid_diffusion_coefficient(self):
        """ Test Species.set_diffusion_coefficient method with an invalid diffusion_coefficient. """
        test_dcs = [None, [1]]
        for test_dc in test_dcs:
            with self.subTest(diffusion_coefficient=test_dc):
                with self.assertRaises(SpeciesError):
                    self.species.set_diffusion_coefficient(diffusion_coefficient=test_dc)

    def test_set_diffusion_coefficient__negative_diffusion_coefficient(self):
        """ Test Species.set_diffusion_coefficient method with negative diffusion_coefficient. """
        with self.assertRaises(SpeciesError):
            self.species.set_diffusion_coefficient(diffusion_coefficient=-1)

    def test_validate__invalid_name(self):
        """ Test Species.validate with an invalid name. """
        test_names = ["", None, 0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(SpeciesError):
                    self.species.name = test_name
                    self.species.validate()

    def test_validate__invalid_diffusion_coefficient(self):
        """ Test Species.validate with an invalid diffusion_coefficient. """
        test_dcs = [None, [0]]
        for test_dc in test_dcs:
            with self.subTest(initial_value=test_dc):
                with self.assertRaises(SpeciesError):
                    self.species.diffusion_coefficient = test_dc
                    self.species.validate()

    def test_validate__invalid_diffusion_coefficient_value(self):
        """ Test Species.validate with an invalid diffusion_coefficient value. """
        test_dcs = [-0.5, -1, -3, -5]
        for test_dc in test_dcs:
            with self.subTest(initial_value=test_dc):
                with self.assertRaises(SpeciesError):
                    self.species.diffusion_coefficient = test_dc
                    self.species.validate()

    def test_validate__invalid_restrict_to(self):
        """ Test Species.validate with an invalid restrict_to. """
        test_rts = ["5", 1, 0.5, [], (1, 2), {'x': 1}]
        for test_rt in test_rts:
            with self.subTest(restrict_to=test_rt):
                with self.assertRaises(SpeciesError):
                    self.species.restrict_to = test_rt
                    self.species.validate()

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        import time
        from datetime import datetime
        start = time.time()
        self.species.validate()
        tic = datetime.utcfromtimestamp(time.time() - start)
        print(f"Total time to run validate: {tic.strftime('%M mins %S secs %f msecs')}")
