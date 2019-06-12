import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
import warnings
import uuid


import numpy
import scipy.io
import scipy.sparse

from spatialpy.Model import *

import inspect

try:
    # This is only needed if we are running in an Ipython Notebook
    import IPython.display
except:
    pass


import pickle
import json
import functools

# module-level variable to for javascript export in IPython/Jupyter notebooks
__spatialpy_javascript_libraries_loaded = False
def load_javascript_libraries():
    global __spatialpy_javascript_libraries_loaded
    if not __spatialpy_javascript_libraries_loaded:
        __spatialpy_javascript_libraries_loaded = True
        import os.path
        import IPython.display
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/three.js')) as fd:
            bufa = fd.read()
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/render.js')) as fd:
            bufb = fd.read()
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/OrbitControls.js')) as fd:
            bufc = fd.read()
        IPython.display.display(IPython.display.HTML('<script>'+bufa+bufc+bufb+'</script>'))


def deprecated(func):
    '''This is a decorator which can be used to mark functions
     as deprecated. It will result in a warning being emitted
     when the function is used.'''

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
             "Call to deprecated function {}.".format(func.__name__),
             category=DeprecationWarning,
             filename=func.func_code.co_filename,
             lineno=func.func_code.co_firstlineno + 1
         )
        return func(*args, **kwargs)
    return new_func


class Solver:
    """ Abstract class for spatialpy solvers. """

    def __init__(self, model, solver_path=None, report_level=0, model_file=None, sopts=None):
        """ Constructor. """
        #TODO: fix class checking
        #if not isinstance(model, Model):
        #    raise SimulationError("Solver constructors must take a Model as an argument.")
        #if not issubclass(self.__class__, Solver):
        #    raise SimulationError("Solver classes must be a subclass of SpatialPy.Solver.")
        if not hasattr(self, 'NAME'):
            raise SimulationError("Solver classes must implement a NAME attribute.")

        self.model = model
        self.is_compiled = False
        self.report_level = report_level
        self.model_file = model_file
        self.infile_name = None
        self.delete_infile = False
        self.model_name = self.model.name
        self.solver_base_dir = None

        # For the remote execution
        self.temp_urdme_root = None

        self.SpatialPy_ROOT =  os.path.dirname(os.path.abspath(__file__))+"/spatialpy"

        #print("solver_path={0}".format(solver_path))
        if solver_path is None or solver_path == "":
            self.SpatialPy_BUILD = self.SpatialPy_ROOT + '/build/'
        else:
            self.SpatialPy_BUILD = solver_path + '/build/'
            os.environ['SOLVER_ROOT'] = solver_path


    def __del__(self):
        """ Deconstructor.  Removes the compiled solver."""
        try:
            if self.solver_base_dir is not None:
                try:
                    shutil.rmtree(self.solver_base_dir)
                except OSError as e:
                    print("Could not delete '{0}'".format(self.solver_base_dir))
            if self.temp_spatialpy_root is not None:
                try:
                    shutil.rmtree(self.temp_spatialpy_root)
                except OSError as e:
                    print("Could not delete '{0}'".format(self.temp_urdme_root))
        except Exception as e:
            print("__del__ failed: {0}".format(e))



    def compile(self):
        """ Compile the model."""

        # Create a unique directory each time call to compile.
        self.solver_base_dir = tempfile.mkdtemp(dir=os.environ.get('SPATIALPY_TMPDIR'))
        self.solver_dir = self.solver_base_dir + '/.spatialpy/'

        if self.report_level >= 1:
            print("Compiling Solver")

        if os.path.isdir(self.solver_dir):
            try:
                shutil.rmtree(self.solver_dir)
            except OSError as e:
                pass
        try:
            os.mkdir(self.solver_dir)
        except Exception as e:
            pass

        # Write the propensity file
        self.propfilename = self.model_name + '_generated_model'
        if self.model_file == None:
            prop_file_name = self.solver_dir + self.propfilename + '.c'
            if self.report_level > 1:
                print("Creating propensity file {0}".format(prop_file_name))
            self.create_propensity_file(file_name=prop_file_name)
        else:
            cmd = " ".join(['cp', self.model_file, self.solver_dir + self.propfilename + '.c'])
            if self.report_level > 1:
                print(cmd)
            subprocess.call(cmd, shell=True)

        # Build the solver
        makefile = 'Makefile.' + self.NAME
        cmd = " ".join([ 'cd', self.solver_base_dir , ';', 'make', '-f', self.SpatialPy_BUILD + makefile, 'SpatialPy_ROOT=' + self.SpatialPy_ROOT, 'SpatialPy_MODEL=' + self.propfilename])
        if self.report_level > 1:
            print("cmd: {0}\n".format(cmd))
        try:
            handle = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            return_code = handle.wait()
        except OSError as e:
            print("Error, execution of compilation raised an exception: {0}".format(e))
            print("cmd = {0}".format(cmd))
            raise SimulationError("Compilation of solver failed")

        if return_code != 0:
            try:
                print(handle.stdout.read())
                print(handle.stderr.read())
            except Exception as e:
                pass
            raise SimulationError("Compilation of solver failed, return_code={0}".format(return_code))

        if self.report_level > 1:
            print(handle.stdout.read())
            print(handle.stderr.read())

        self.is_compiled = True


    def run(self, number_of_trajectories=1, seed=None, input_file=None, loaddata=False):
        """ Run one simulation of the model.
        number_of_trajectories: How many trajectories should be run.
        seed: the random number seed (incremented by one for multiple runs).
        input_file: the filename of the solver input data file .
        loaddata: boolean, should the result object load the data into memory on creation.
        Returns:
            Result object.
                or, if number_of_trajectories > 1
            a list of Result objects
        """
        if number_of_trajectories > 1:
            result_list = []
        # Check if compiled, call compile() if not.
        if not self.is_compiled:
            self.compile()

        # Execute the solver
        for run_ndx in range(number_of_trajectories):
            outfile = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('SpatialPy_TMPDIR'))
            outfile.close()
            urdme_solver_cmd = [self.solver_dir + self.propfilename + '.' + self.NAME, self.infile_name, outfile.name]

            if seed is not None:
                urdme_solver_cmd.append(str(seed+run_ndx))
            if self.report_level > 1:
                print('cmd: {0}\n'.format(urdme_solver_cmd))
            stdout = ''
            stderr = ''
            try:
                if self.report_level >= 1:  #stderr & stdout to the terminal
                    handle = subprocess.Popen(urdme_solver_cmd)
                else:
                    handle = subprocess.Popen(urdme_solver_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                    stdout, stderr = handle.communicate()
                return_code = handle.wait()
            except OSError as e:
                print("Error, execution of solver raised an exception: {0}".format(e))
                print("urdme_solver_cmd = {0}".format(urdme_solver_cmd))

            if return_code != 0:
                if self.report_level >= 1:
                    try:
                        print(stderr)
                        print(stdout)
                    except Exception as e:
                        pass
                print("urdme_solver_cmd = {0}".format(urdme_solver_cmd))
                raise SimulationError("Solver execution failed, return code = {0}".format(return_code))


            #Load the result from the hdf5 output file.
            try:
                result = Result(self.model, outfile.name, loaddata=loaddata)
                result["Status"] = "Sucess"
                result.stderr = stderr
                result.stdout = stdout
                if number_of_trajectories > 1:
                    result_list.append(result)
                else:
                    return result
            except Exception as e:
                exc_info = sys.exc_info()
                os.remove(outfile.name)
                raise(exc_info[1], None, exc_info[2])

        return result_list


    def create_propensity_file(self, file_name=None):
        """ Generate the C propensity file that is used to compile the solvers.
            Only mass action propensities are supported.
        """

        template = open(os.path.abspath(os.path.dirname(__file__)) + '/data/propensity_file_template.c', 'r')
        propfile = open(file_name, "w")
        propfilestr = template.read()

        speciesdef = ""
        i = 0
        for S in self.model.listOfSpecies:
            speciesdef += "#define " + S + " " + "x[" + str(i) + "]" + "\n"
            i += 1

        propfilestr = propfilestr.replace("__DEFINE_SPECIES__", speciesdef)

        propfilestr = propfilestr.replace("__NUMBER_OF_REACTIONS__", str(self.model.get_num_reactions()))
        propfilestr = propfilestr.replace("__NUMBER_OF_SPECIES__", str(self.model.get_num_species()))
        propfilestr = propfilestr.replace("__NUMBER_OF_VOXELS__", str(self.model.mesh.get_num_voxels()))

        # Create defines for the DataFunctions.
        data_fn_str = ""
        i = 0
        for d in self.model.listOfDataFunctions:
            if d.name is None:
                raise SimulationError("DataFunction {0} does not have a name attributed defined.".format(i))
            data_fn_str += "#define " + d.name + " data[" + str(i) + "]\n"
            i += 1
        propfilestr = propfilestr.replace("__DEFINE_DATA_FUNCTIONS__", str(data_fn_str))

        # Make sure all paramters are evaluated to scalars before we write them to the file.
        self.model.resolve_parameters()
        parameters = ""
        for p in self.model.listOfParameters:
            parameters += "const double " + p + " = " + str(self.model.listOfParameters[p].value) + ";\n"
        propfilestr = propfilestr.replace("__DEFINE_PARAMETERS__", str(parameters))


        # Reactions
        funheader = "double __NAME__(const int *x, double t, const double vol, const double *data, int sd)"
        #funheader = "double __NAME__(const int *x, double t, const double vol, const double *data, int sd, int voxel, int *xx, const size_t *irK, const size_t *jcK, const double *prK)"

        funcs = ""
        funcinits = ""
        i = 0
        for R in self.model.listOfReactions:
            func = ""
            rname = self.model.listOfReactions[R].name
            func += funheader.replace("__NAME__", rname) + "\n{\n"
            if self.model.listOfReactions[R].restrict_to == None or (isinstance(self.model.listOfReactions[R].restrict_to, list) and len(self.model.listOfReactions[R].restrict_to) == 0):
                func += self.model.listOfReactions[R].propensity_function
            else:
                func += "if("
                if isinstance(self.model.listOfReactions[R].restrict_to, list) and len(self.model.listOfReactions[R].restrict_to) > 0:
                    for sd in self.model.listOfReactions[R].restrict_to:
                        func += "sd == " + str(sd) + "||"
                    func = func[:-2]
                elif isinstance(self.model.listOfReactions[R].restrict_to, int):
                    func += "sd == " +  str(self.model.listOfReactions[R].restrict_to)
                else:
                    raise SimulationError("When restricting reaction to subdomains, you must specify either a list or an int")
                func += "){\n"
                func += self.model.listOfReactions[R].propensity_function

                func += "\n}else{"
                func += "\n\treturn 0.0;}"


            func += "\n}"
            funcs += func + "\n\n"
            funcinits += "    ptr[" + str(i) + "] = " + rname + ";\n"
            i += 1

        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__", funcs)
        propfilestr = propfilestr.replace("__DEFINE_PROPFUNS__", funcinits)
        propfile.write(propfilestr)
        propfile.close()





class DataFunction():
    """ Abstract class used to constuct the data vector. """
    name = None
    def __init__(self, name=None):
        if name is not None:
            self.name = name
        if self.name is None:
            raise Exception("DataFunction must have a 'name'")

    def map(self, x):
        """ map() takes the coordinate 'x' and returns a list of doubles.
        Args:
            x: a list of 3 ints.
        Returns:
            a list of floats.
        """
        raise Exception("DataFunction.map() not implemented.")





