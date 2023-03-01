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
import math
import numpy
import unittest

from spatialpy.core.timespan import TimeSpan
from spatialpy.core.spatialpyerror import TimespanError

class TestTimeSpan(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for spatialpy.TimeSpan.
    ################################################################################################
    '''
    def setUp(self):
        """ Setup a clean valid timespan for testing. """
        self.tspan = TimeSpan.linspace(t=10, num_points=11, timestep_size=0.001)

    def set_timesteps(self, tspan, timestep_size):
        items_diff = numpy.diff(tspan)
        output_interval = items_diff[0]
        num_steps = len(items_diff)
        
        output_freq = output_interval / timestep_size
        num_timesteps = math.ceil(num_steps * output_freq)

        output_steps = numpy.arange(
            0, num_timesteps + timestep_size, output_freq
        )
        output_steps = numpy.unique(numpy.round(output_steps).astype(int))
        items = numpy.zeros((output_steps.size), dtype=float)
        for i, step in enumerate(output_steps):
            items[i] = step * timestep_size
        return items

    def test_constructor(self):
        """ Test the TimeSpan constructor. """
        test_tspan = self.set_timesteps(numpy.linspace(0, 20, 401), 0.001)
        tspan = TimeSpan(numpy.linspace(0, 20, 401), timestep_size=0.001)
        self.assertEqual(tspan, test_tspan)
        self.assertEqual(tspan.timestep_size, 0.001)

    def test_constructor__valid_data_structures(self):
        """ Test the TimeSpan constructor with valid data structures. """
        test_tspans = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
            range(11)
        ]
        for raw_tspan in test_tspans:
            with self.subTest(tspan=raw_tspan, tspan_type=type(raw_tspan)):
                test_tspan = self.set_timesteps(numpy.array(raw_tspan), 0.001)
                tspan = TimeSpan(raw_tspan, timestep_size=0.001)
                self.assertEqual(tspan, test_tspan)

    def test_constructor__invalid_items(self):
        """ Test the TimeSpan constructor with an invalid data structure type. """
        test_tspans = [
            None, "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]", 20, 50.5,
            set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        ]
        for test_tspan in test_tspans:
            with self.subTest(items=test_tspan):
                with self.assertRaises(TimespanError):
                    TimeSpan(test_tspan, timestep_size=0.001)

    def test_constructor__timestep_size_none(self):
        """ Test the TimeSpan constructor when timestep_size is omitted or set to None. """
        test_tspan = TimeSpan(numpy.linspace(0, 30, 301))
        self.assertEqual(test_tspan.timestep_size, 0.1)

    def test_constructor__invaild_timestep_size(self):
        """ Test the TimeSpan constructor when timestep_size is of an invalid type. """
        test_tsss = ["1", [0.001]]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    TimeSpan(numpy.linspace(0, 30, 301), timestep_size=test_tss)

    def test_constructor__invaild_timestep_size_value(self):
        """ Test the TimeSpan constructor when timestep_size is an invalid value. """
        test_tsss = [0, -0.5, -1, -5]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    TimeSpan(numpy.linspace(0, 30, 301), timestep_size=test_tss)

    def test_constructor__timestep_size_too_large(self):
        """ Test the TimeSpan constructor when the timestep_size > the diff of points. """
        with self.assertRaises(TimespanError):
            TimeSpan(numpy.arange(0, 20.5, 0.5), timestep_size=1)

    def test_linspace(self):
        """ Test TimeSpan.linspace. """
        tspan = TimeSpan.linspace(t=30, num_points=301, timestep_size=0.001)
        test_tspan = self.set_timesteps(numpy.linspace(0, 30, 301), 0.001)
        self.assertEqual(tspan, test_tspan)
        self.assertEqual(tspan.timestep_size, 0.001)

    def test_linspace__no_t(self):
        """ Test TimeSpan.linspace without passing t. """
        tspan = TimeSpan.linspace(num_points=201, timestep_size=0.001)
        test_tspan = self.set_timesteps(numpy.linspace(0, 20, 201), 0.001)
        self.assertEqual(tspan, test_tspan)

    def test_linspace__invalid_t(self):
        """ Test TimeSpan.linspace with invalid t. """
        test_values = [None, "5", 0, -0.5, -1, -2, -5, -10, [20.5]]
        for test_val in test_values:
            with self.subTest(t=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.linspace(t=test_val, num_points=301, timestep_size=0.001)

    def test_linspace__no_num_points(self):
        """ Test TimeSpan.linspace without passing num_points. """
        tspan = TimeSpan.linspace(t=30, timestep_size=0.001)
        test_tspan = self.set_timesteps(numpy.linspace(0, 30, int(30 / 0.05) + 1), 0.001)
        self.assertEqual(tspan, test_tspan)

    def test_linspace__invalid_num_points(self):
        """ Test TimeSpan.linspace with invalid num_points. """
        test_values = ["5", 0, -1, -2, -5, -10, 4.5, [40]]
        for test_val in test_values:
            with self.subTest(num_points=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.linspace(t=30, num_points=test_val, timestep_size=0.001)

    def test_linspace__timestep_size_none(self):
        """ Test TimeSpan.linspace when timestep_size is omitted or set to None. """
        test_tspan = TimeSpan.linspace(t=30, num_points=301)
        self.assertEqual(test_tspan.timestep_size, 0.1)

    def test_linspace__invaild_timestep_size(self):
        """ Test TimeSpan.linspace when timestep_size is of an invalid type. """
        test_tsss = ["1", [0.001]]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    TimeSpan.linspace(t=30, num_points=301, timestep_size=test_tss)

    def test_linspace__invaild_timestep_size_value(self):
        """ Test TimeSpan.linspace when timestep_size is an invalid value. """
        test_tsss = [0, -0.5, -1, -5]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    TimeSpan.linspace(t=30, num_points=301, timestep_size=test_tss)

    def test_linspace__timestep_size_too_large(self):
        """ Test TimeSpan.linspace when the timestep_size > the diff of points. """
        with self.assertRaises(TimespanError):
            TimeSpan.linspace(t=20, num_points=41, timestep_size=1)

    def test_linspace__no_args(self):
        """ Test TimeSpan.linspace without passing any args. """
        tspan = TimeSpan.linspace()
        test_tspan = self.set_timesteps(numpy.linspace(0, 20, 401), 0.05)
        self.assertEqual(tspan, test_tspan)
        self.assertEqual(tspan.timestep_size, 0.05)

    def test_arange(self):
        """ Test TimeSpan.arange. """
        tspan = TimeSpan.arange(0.1, t=30, timestep_size=0.001)
        test_tspan = self.set_timesteps(numpy.arange(0, 30.1, 0.1), 0.001)
        self.assertEqual(tspan, test_tspan)
        self.assertEqual(tspan.timestep_size, 0.001)

    def test_arange__no_t(self):
        """ Test TimeSpan.arange. """
        tspan = TimeSpan.arange(0.1, timestep_size=0.001)
        test_tspan = self.set_timesteps(numpy.arange(0, 20.1, 0.1), 0.001)
        self.assertEqual(tspan, test_tspan)

    def test_arange__invalid_t(self):
        """ Test TimeSpan.arange with invalid t. """
        test_values = [None, "5", 0, -0.5, -1, -2, -5, -10, [20.5]]
        for test_val in test_values:
            with self.subTest(t=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.arange(0.1, t=test_val, timestep_size=0.001)

    def test_arange__invalid_increment(self):
        """ Test TimeSpan.arange with invalid increment type. """
        test_values = [None, "0.05", 0, -1, -2, -5, -10, [0.05]]
        for test_val in test_values:
            with self.subTest(imcrement=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.arange(test_val, t=30, timestep_size=0.001)

    def test_arange__timestep_size_none(self):
        """ Test TimeSpan.arange when timestep_size is omitted or set to None. """
        test_tspan = TimeSpan.arange(0.1, t=30)
        self.assertEqual(test_tspan.timestep_size, 0.1)

    def test_arange__invaild_timestep_size(self):
        """ Test TimeSpan.arange when timestep_size is of an invalid type. """
        test_tsss = ["1", [0.001]]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    TimeSpan.arange(0.1, t=30, timestep_size=test_tss)

    def test_arange__invaild_timestep_size_values(self):
        """ Test TimeSpan.arange when timestep_size is an invalid values. """
        test_tsss = [0, -0.5, -1, -5]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    TimeSpan.arange(0.1, t=30, timestep_size=test_tss)

    def test_arange__timestep_size_too_large(self):
        """ Test TimeSpan.arange when the timestep_size > the diff of points. """
        with self.assertRaises(TimespanError):
            TimeSpan.arange(0.5, t=20, timestep_size=1)

    def test_validate__list(self):
        """ Test TimeSpan.validate with list data structure. """
        raw_tspan = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.tspan.items = raw_tspan
        self.tspan.validate(coverage="all")
        test_tspan = self.set_timesteps(numpy.array(raw_tspan), 0.001)
        self.assertIsInstance(self.tspan.items, numpy.ndarray)
        self.assertEqual(self.tspan, test_tspan)

    def test_validate__tuple(self):
        """ Test TimeSpan.validate with tuple data structure. """
        raw_tspan = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        self.tspan.items = raw_tspan
        self.tspan.validate(coverage="all")
        test_tspan = self.set_timesteps(numpy.array(raw_tspan), 0.001)
        self.assertIsInstance(self.tspan.items, numpy.ndarray)
        self.assertEqual(self.tspan, test_tspan)

    def test_validate__range(self):
        """ Test TimeSpan.validate with range data structure. """
        raw_tspan = range(11)
        self.tspan.items = raw_tspan
        self.tspan.validate(coverage="all")
        test_tspan = self.set_timesteps(numpy.array(raw_tspan), 0.001)
        self.assertIsInstance(self.tspan.items, numpy.ndarray)
        self.assertEqual(self.tspan, test_tspan)

    def test_validate__invalid_type(self):
        """ Test TimeSpan.validate with an invalid data structure type. """
        test_tspans = [
            None, "50", 20, 40.5, set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        ]
        for test_tspan in test_tspans:
            if test_tspan is not None:
                self.setUp()
            with self.subTest(items=test_tspan):
                with self.assertRaises(TimespanError):
                    self.tspan.items = test_tspan
                    self.tspan.validate(coverage="all")

    def test_validate__empty_timespan(self):
        """ Test TimeSpan.validate with an empty data structure. """
        test_tspans = [[], ()]
        for test_tspan in test_tspans:
            if test_tspan != []:
                self.setUp()
            with self.subTest(items=test_tspan):
                with self.assertRaises(TimespanError):
                    self.tspan.items = test_tspan
                    self.tspan.validate(coverage="all")

    def test_validate__all_same_values(self):
        """ Test TimeSpan.validate with an empty data structure. """
        with self.assertRaises(TimespanError):
            self.tspan.items = [2, 2, 2, 2, 2, 2, 2, 2, 2]
            self.tspan.validate(coverage="all")

    def test_validate__negative_start(self):
        """ Test TimeSpan.validate with an initial time < 0. """
        with self.assertRaises(TimespanError):
            self.tspan.items = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            self.tspan.validate(coverage="all")

    def test_validate__non_uniform_timespan(self):
        """ Test TimeSpan.validate with a non-uniform timespan. """
        with self.assertRaises(TimespanError):
            self.tspan.items = [2, 1, 3, 4, 5, 6, 7, 8, 9, 10]
            self.tspan.validate(coverage="all")

    def test_validate__invalid_timestep_size(self):
        """ Test TimeSpan.validate when timestep_size is of an invalid type. """
        test_tsss = [None, "1", [0.001]]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    self.tspan.timestep_size = test_tss
                    self.tspan.validate(coverage="all")

    def test_validate__invalid_timestep_size_values(self):
        """ Test TimeSpan.validate when timestep_size is an invalid values. """
        test_tsss = [0, -0.5, -1, -5]
        for test_tss in test_tsss:
            with self.subTest(timestep_size=test_tss):
                with self.assertRaises(TimespanError):
                    self.tspan.timestep_size = test_tss
                    self.tspan.validate(coverage="all")

    def test_validate__invalid_output_freq(self):
        """ Test TimeSpan.validate when output_freq is of an invalid type. """
        test_opfs = [None, "5", [0.5]]
        for test_opf in test_opfs:
            with self.subTest(output_freq=test_opf):
                with self.assertRaises(TimespanError):
                    self.tspan.output_freq = test_opf
                    self.tspan.validate(coverage="all")

    def test_validate__invalid_output_freq_values(self):
        """ Test TimeSpan.validate when output_freq is an invalid value. """
        test_opfs = [0, -0.5, -1, -5]
        for test_opf in test_opfs:
            with self.subTest(output_freq=test_opf):
                with self.assertRaises(TimespanError):
                    self.tspan.output_freq = test_opf
                    self.tspan.validate(coverage="all")

    def test_validate__output_freq_lessthan_timestep_size(self):
        """ Test TimeSpan.validate when output_freq < timestep_size. """
        with self.assertRaises(TimespanError):
            self.tspan.timestep_size = 0.5
            self.tspan.output_freq = 0.1
            self.tspan.validate(coverage="all")

    def test_validate__invalid_num_timesteps(self):
        """ Test TimeSpan.validate when num_timesteps is of an invalid type. """
        test_nsts = [None, "5", 0.5, [5]]
        for test_nst in test_nsts:
            with self.subTest(num_timestep=test_nst):
                with self.assertRaises(TimespanError):
                    self.tspan.num_timesteps = test_nst
                    self.tspan.validate(coverage="all")

    def test_validate__invalid_num_timesteps_value(self):
        """ Test TimeSpan.validate when num_timesteps is an invalid value. """
        test_nsts = [0, -1, -5]
        for test_nst in test_nsts:
            with self.subTest(num_timesteps=test_nst):
                with self.assertRaises(TimespanError):
                    self.tspan.num_timesteps = test_nst
                    self.tspan.validate(coverage="all")

    def test_validate__invalid_output_steps(self):
        """ Test TimeSpan.validate when output_steps is of an invalid type. """
        test_opss = [None, "5", 5, 0.5, [5], {"x":5}, (6,2)]
        for test_ops in test_opss:
            with self.subTest(output_steps=test_ops):
                with self.assertRaises(TimespanError):
                    self.tspan.output_steps = test_ops
                    self.tspan.validate(coverage="all")

    def test_validate__invalid_output_steps_value(self):
        """ Test TimeSpan.validate when output_steps is an invalid value. """
        with self.assertRaises(TimespanError):
            self.tspan.output_steps = numpy.array([6,2])
            self.tspan.validate(coverage="all")
