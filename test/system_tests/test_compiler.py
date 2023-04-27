#!/usr/bin/env python3
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
''' Tests from the compiler. '''
import os
import sys
import unittest

sys.path.insert(1, os.path.abspath(os.getcwd()))
from test.models.birth_death import create_birth_death # pylint: disable=wrong-import-position

class TestCompiler(unittest.TestCase):
    '''
    ################################################################################################
    System test for compiler.
    ################################################################################################
    '''
    def setUp(self):
        self.model = create_birth_death()

    def test_constructor(self):
        """ Test the compiler. """
        _ = self.model.run()

if __name__ == '__main__':
    unittest.main()
