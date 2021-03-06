import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import spatialpy

class diffusion_debug(spatialpy.Model):

    def __init__(self, model_name='diffusion_debug_test'):
        spatialpy.Model.__init__(self, model_name)

        A = spatialpy.Species(name='A', diffusion_constant=0.01)
        self.add_species([A])

        self.mesh = spatialpy.Mesh.create_2D_domain(
            xlim=[-1, 1],
            ylim=[-1, 1],
            nx=50,
            ny=50,
            type_id=1.0,
            mass=1.0,
            nu=1.0,
            fixed=True,
            rho0=1.0,
            c0=1.0,
            P0=1.0,
            )

        self.add_initial_condition(spatialpy.PlaceInitialCondition(A,
                                   1000, [0, 0, 0]))

        self.timestep_size = 0.1
        self.num_timesteps = 10
        self.output_freq = 1


class test:

    def __init__(self):
        pass

    def test_different_seed(self):
        """ Test that the output is different if different seeds are given (default set on C level). """

        self.model = diffusion_debug()
        solver = spatialpy.Solver(self.model)
        resultOne = solver.run(seed=1)
        resultTwo = solver.run(seed=1)
        #resultOne = spatialpy.Result()
        #resultTwo = spatialpy.Result()
        resultOne.result_dir = None
        resultTwo.result_dir = None
        print('using == : {}'.format(resultOne == 1))
        #print('using __eq__() : {}'.format(resultOne.__eq__(resultTwo)))
        print('result 1 : {}'.format(type(resultOne)))
        print('result 2 : {}'.format(type(resultTwo)))


class init:

    test = test()
    test.test_different_seed()