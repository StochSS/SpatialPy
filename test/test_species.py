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
from spatialpy import Species, Parameter
from spatialpy import SpeciesError

class TestSpecies(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for spatialpy.Species.
    ################################################################################################
    '''
    def test_constructor(self):
        """ Test the Species constructor. """
        species = Species(name="test_species", diffusion_coefficient=0)
        self.assertEqual(species.name, "test_species")
        self.assertEqual(species.diffusion_coefficient, 0)


    def test_constructor__no_name(self):
        """ Test the Species constructor without name. """
        with self.assertRaises(SpeciesError):
            species = Species(diffusion_coefficient=0)


    def test_constructor__name_not_str(self):
        """ Test the Species constructor with non-str name. """
        with self.assertRaises(SpeciesError):
            species = Species(name=0, diffusion_coefficient=0)


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
        self.assertEqual(species.diffusion_coefficient, test_parameter.expression)


    def test_constructor__diffusion_coefficient_not_int_or_float(self):
        """ Test the Species constructor with non-int or non-float diffusion_coefficient. """
        with self.assertRaises(SpeciesError):
            species = Species(name="test_species", diffusion_coefficient="0")


    def test___str___(self):
        """ Test Species.__str__ method. """
        species = Species(name="test_species", diffusion_coefficient=0)
        self.assertIsInstance(str(species), str)


    def test_set_diffusion_coefficient(self):
        """ Test Species.set_diffusion_coefficient method. """
        species = Species(name="test_species", diffusion_coefficient=0)
        species.set_diffusion_coefficient(diffusion_coefficient=1)
        self.assertEqual(species.diffusion_coefficient, 1)


    def test_set_diffusion_coefficient__negative_diffusion_coefficient(self):
        """ Test Species.set_diffusion_coefficient method with negative diffusion_coefficient. """
        species = Species(name="test_species", diffusion_coefficient=0)
        with self.assertRaises(SpeciesError):
            species.set_diffusion_coefficient(diffusion_coefficient=-1)


    def test_set_diffusion_coefficient__diffusion_coefficient_not_int_or_float(self):
        """ Test Species.set_diffusion_coefficient method with non-int or non-float diffusion_coefficient. """
        species = Species(name="test_species", diffusion_coefficient=0)
        with self.assertRaises(SpeciesError):
            species.set_diffusion_coefficient(diffusion_coefficient="1")
