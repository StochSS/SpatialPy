import unittest, sys, os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mode', default='develop', choices=['develop', 'release'],
                    help='Run tests in develop mode or release mode.')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.mode == 'develop':
        print('Running tests in develop mode. Appending repository directory to system path.')
        sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    import test_solver
    #import test_mincde

    modules = [
        test_solver,
        #test_mincde,
    ]

    for module in modules:
        suite = unittest.TestLoader().loadTestsFromModule(module)
        runner = unittest.TextTestRunner(failfast=args.mode == 'develop')

        print("Executing: {}".format(module))
        result = runner.run(suite)
        print('=' * 70)
        if not result.wasSuccessful():
            sys.exit(not result.wasSuccessful())
