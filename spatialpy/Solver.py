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

from model import *

import inspect

try:
    # This is only needed if we are running in an Ipython Notebook
    import IPython.display
except:
    pass

try:
    import h5py
except:
    raise Exception("SpatialPy requires h5py.")

import pickle
import json
import functools

# module-level variable to for javascript export in IPython/Jupyter notebooks
__spatialpy_javascript_libraries_loaded = False
def load_pyurdme_javascript_libraries():
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


class SpatialPySolver:
    """ Abstract class for spatialpy solvers. """

    def __init__(self, model, solver_path=None, report_level=0, model_file=None, sopts=None):
        """ Constructor. """
        if not isinstance(model, SpatialPyModel):
            raise SpatialPyError("URDMEsolver constructors must take a SpatialPyModel as an argument.")
        if not issubclass(self.__class__, URDMESolver):
            raise SpatialPyError("Solver classes must be a subclass of SpatialPySolver.")
        if not hasattr(self, 'NAME'):
            raise SpatialPyError("Solver classes must implement a NAME attribute.")

        self.model = model
        self.is_compiled = False
        self.report_level = report_level
        self.model_file = model_file
        self.infile_name = None
        self.delete_infile = False
        self.model_name = self.model.name
        self.solver_base_dir = None
        if sopts is None:
            self.sopts = [0,0,0]
        else:
            self.sopts = sopts

        # For the remote execution
        self.temp_urdme_root = None

        self.SpatialPy_ROOT =  os.path.dirname(os.path.abspath(__file__))+"/spatialpy"

        #print "solver_path={0}".format(solver_path)
        if solver_path is None or solver_path == "":
            self.SpatialPy_BUILD = self.SpatialPy_ROOT + '/build/'
        else:
            self.SpatialPy_BUILD = solver_path + '/build/'
            os.environ['SOLVER_ROOT'] = solver_path

    def __getstate__(self):
        """ Save the state of the solver, saves all instance variables
            and reads all the files necessary to compile the solver off
            of the file system and stores it in a separate state variable.
            If the solver model files is specified, it saves that too.
            This is used by Pickle.
        """
        ret = {}
        # Save the instance variables
        ret['vars'] = self.__dict__.copy()
        # The model object is not picklabe due to the Swig-objects from Dolfin
        #ret['vars']['model'] = None
        ret['vars']['is_compiled'] = False
        # Create temp root
        tmproot = tempfile.mkdtemp(dir=os.environ.get('PYURDME_TMPDIR'))
        # Get the propensity file
        model_file = tmproot+'/'+self.model_name + '_pyurdme_generated_model'+ '.c'
        ret['model_file'] = os.path.basename(model_file)
        if self.model_file == None:
            self.create_propensity_file(file_name=model_file)
        else:
            subprocess.call('cp '+self.model_file+' '+model_file, shell=True)
        # Get the solver source files
        os.mkdir(tmproot+'/include')
        os.mkdir(tmproot+'/src')
        os.mkdir(tmproot+'/src/'+self.NAME)
        #TODO: what if solverdir is not the same as URDME_ROOT ?
        subprocess.call('cp '+self.URDME_ROOT+'/src/*.c '+tmproot+'/src/', shell=True)
        subprocess.call('cp '+self.URDME_ROOT+'/src/'+self.NAME+'/*.* '+tmproot+'/src/'+self.NAME+'/', shell=True)
        subprocess.call('cp '+self.URDME_ROOT+'/include/*.h '+tmproot+'/include/', shell=True)
        #TODO: get the include files from solvers not in the default path (none currently implement this).
        # Get the Makefile
        os.mkdir(tmproot+'/build')
        subprocess.call('cp '+self.URDME_BUILD+'Makefile.'+self.NAME+' '+tmproot+'/build/Makefile.'+self.NAME, shell=True)
        # Get the input file
        input_file = tmproot+'/model_input.mat'
        ret['input_file'] = os.path.basename(input_file)
        self.serialize(filename=input_file, report_level=self.report_level)
        ##
        origwd = os.getcwd()
        os.chdir(tmproot)
        tarname = tmproot+'/'+self.NAME+'.tar.gz'
        subprocess.call('tar -czf '+tarname+' src include build '+os.path.basename(input_file)+' '+os.path.basename(model_file), shell=True)
        with open(tarname, 'r') as f:
            ret['SolverFiles'] = f.read()
        os.chdir(origwd)
        shutil.rmtree(tmproot)
        # return the state
        return ret

    def __setstate__(self, state):
        """ Set all instance variables for the object, and create a unique temporary
            directory to store all the solver files.  URDME_BUILD is set to this dir,
            and is_compiled is always set to false.  This is used by Pickle.
        """
        # 0. restore the instance variables
        for key, val in state['vars'].iteritems():
            self.__dict__[key] = val
        # 1. create temporary directory = SPATIALPY_ROOT
        self.temp_spatialpy_root = tempfile.mkdtemp(dir=os.environ.get('SPATIALPY_TMPDIR'))
        self.SPATIALPY_ROOT = self.temp_spatialpy_root
        self.SPATIALPY_BUILD = self.temp_spatialpy_root+'/build/'
        origwd = os.getcwd()
        os.chdir(self.temp_spatialpy_root)
        tarname = self.temp_spatialpy_root+'/'+self.NAME+'.tar.gz'
        with open(tarname, 'wd') as f:
            f.write(state['SolverFiles'])
        subprocess.call('tar -zxf '+tarname, shell=True)
        os.chdir(origwd)
        # Model File
        self.model_file = self.temp_spatialpy_root+'/'+state['model_file']
        # Input File
        self.infile_name = self.temp_spatialpy_root+'/'+state['input_file']


    def __del__(self):
        """ Deconstructor.  Removes the compiled solver."""
        if self.delete_infile:
            try:
                os.remove(self.infile_name)
            except OSError as e:
                print "Could not delete '{0}'".format(self.infile_name)
        if self.solver_base_dir is not None:
            try:
                shutil.rmtree(self.solver_base_dir)
            except OSError as e:
                print "Could not delete '{0}'".format(self.solver_base_dir)
        if self.temp_spatialpy_root is not None:
            try:
                shutil.rmtree(self.temp_spatialpy_root)
            except OSError as e:
                print "Could not delete '{0}'".format(self.temp_urdme_root)


    def serialize(self, filename=None, report_level=0, sopts=None):
        """ Write the datastructures needed by the the core URDME solvers to a .mat input file. """
        spatialpy_solver_data = self.model.get_solver_datastructure()
        spatialpy_solver_data['report'] = report_level
        if sopts is None:
            spatialpy_solver_data['sopts'] = self.sopts
        else:
            spatialpy_solver_data['sopts'] = sopts

        self.model.validate(spatialpy_solver_data)
        scipy.io.savemat(filename, spatialpy_solver_data, oned_as='column')


    def compile(self):
        """ Compile the model."""

        # Create a unique directory each time call to compile.
        self.solver_base_dir = tempfile.mkdtemp(dir=os.environ.get('SPATIALPY_TMPDIR'))
        self.solver_dir = self.solver_base_dir + '/.spatialpy/'
        #print "URDMESolver.compile()  self.solver_dir={0}".format(self.solver_dir)

        if self.report_level >= 1:
            print "Compiling Solver"

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
        self.propfilename = self.model_name + '_pyurdme_generated_model'
        if self.model_file == None:
            prop_file_name = self.solver_dir + self.propfilename + '.c'
            if self.report_level > 1:
                print "Creating propensity file {0}".format(prop_file_name)
            self.create_propensity_file(file_name=prop_file_name)
        else:
            cmd = " ".join(['cp', self.model_file, self.solver_dir + self.propfilename + '.c'])
            if self.report_level > 1:
                print cmd
            subprocess.call(cmd, shell=True)

        # Build the solver
        makefile = 'Makefile.' + self.NAME
        cmd = " ".join([ 'cd', self.solver_base_dir , ';', 'make', '-f', self.URDME_BUILD + makefile, 'URDME_ROOT=' + self.URDME_ROOT, 'URDME_MODEL=' + self.propfilename])
        if self.report_level > 1:
            print "cmd: {0}\n".format(cmd)
        try:
            handle = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            return_code = handle.wait()
        except OSError as e:
            print "Error, execution of compilation raised an exception: {0}".format(e)
            print "cmd = {0}".format(cmd)
            raise URDMEError("Compilation of solver failed")

        if return_code != 0:
            try:
                print handle.stdout.read()
                print handle.stderr.read()
            except Exception as e:
                pass
            raise URDMEError("Compilation of solver failed, return_code={0}".format(return_code))

        if self.report_level > 1:
            print handle.stdout.read()
            print handle.stderr.read()

        self.is_compiled = True


    def run(self, number_of_trajectories=1, seed=None, input_file=None, loaddata=False):
        """ Run one simulation of the model.
        number_of_trajectories: How many trajectories should be run.
        seed: the random number seed (incremented by one for multiple runs).
        input_file: the filename of the solver input data file .
        loaddata: boolean, should the result object load the data into memory on creation.
        Returns:
            URDMEResult object.
                or, if number_of_trajectories > 1
            a list of URDMEResult objects
        """
        if number_of_trajectories > 1:
            result_list = []
        # Check if compiled, call compile() if not.
        if not self.is_compiled:
            self.compile()

        if input_file is None:
            if self.infile_name is None or not os.path.exists(self.infile_name):
                # Get temporary input and output files
                infile = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))

                # Write the model to an input file in .mat format
                self.serialize(filename=infile, report_level=self.report_level)
                infile.close()
                self.infile_name = infile.name
                self.delete_infile = True
        else:
            self.infile_name = input_file
            self.delete_infile = False

        if not os.path.exists(self.infile_name):
            raise URDMEError("input file not found.")

        # Execute the solver
        for run_ndx in range(number_of_trajectories):
            outfile = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))
            outfile.close()
            urdme_solver_cmd = [self.solver_dir + self.propfilename + '.' + self.NAME, self.infile_name, outfile.name]

            if seed is not None:
                urdme_solver_cmd.append(str(seed+run_ndx))
            if self.report_level > 1:
                print 'cmd: {0}\n'.format(urdme_solver_cmd)
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
                print "Error, execution of solver raised an exception: {0}".format(e)
                print "urdme_solver_cmd = {0}".format(urdme_solver_cmd)

            if return_code != 0:
                if self.report_level >= 1:
                    try:
                        print stderr, stdout
                    except Exception as e:
                        pass
                print "urdme_solver_cmd = {0}".format(urdme_solver_cmd)
                raise URDMEError("Solver execution failed, return code = {0}".format(return_code))


            #Load the result from the hdf5 output file.
            try:
                result = URDMEResult(self.model, outfile.name, loaddata=loaddata)
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
                raise exc_info[1], None, exc_info[2]

        return result_list


    def create_propensity_file(self, file_name=None):
        """ Generate the C propensity file that is used to compile the URDME solvers.
            Only mass action propensities are supported.
        """

        template = open(os.path.abspath(os.path.dirname(__file__)) + '/data/propensity_file_template.c', 'r')
        propfile = open(file_name, "w")
        propfilestr = template.read()

        speciesdef = ""
        i = 0
        for S in self.model.listOfSpecies:
            speciesdef += "#define " + S + " " + "x[" + str(i) + "]" + "\n"
            speciesdef += "#define " + S + "_INDEX " +  str(i) + "\n"
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
                raise URDMEError("DataFunction {0} does not have a name attributed defined.".format(i))
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
                    raise URDMEError("When restricting reaction to subdomains, you must specify either a list or an int")
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



def spatialpy(model=None, solver='nsm', solver_path="", model_file=None, input_file=None, seed=None, report_level=0):
    """ SPATIALPY solver interface.
        Similar to model.run() the urdme() function provides an interface that is backwards compatiable with the
        previous URDME implementation.
        After sucessful execution, urdme returns a URDMEResults object with the following members:
        U:         the raw copy number output in a matrix with dimension (Ndofs, num_time_points)
        tspan:     the time span vector containing the time points that corresponds to the columns in U
        status:    Sucess if the solver executed without error
        stdout:    the standard ouput stream from the call to the core solver
        stderr:    the standard error stream from the call to the core solver
    """


    #If solver is a subclass of URDMESolver, use it directly.
    if isinstance(solver, (type, types.ClassType)) and  issubclass(solver, URDMESolver):
        sol = solver(model, solver_path, report_level, model_file=model_file)
    elif type(solver) is str:
        if solver == 'nsm':
            from nsmsolver import NSMSolver
            sol = NSMSolver(model, solver_path, report_level, model_file=model_file)
        elif solver == 'nem':
            from nemsolver import NEMSolver
            sol = NEMSolver(model, solver_path, report_level, model_file=model_file)
        else:
            raise URDMEError("Unknown solver: {0}".format(solver_name))
    else:
        raise URDMEError("solver argument to spatialpy() must be a string or a spatialpysolver class object.")

    sol.compile()
    return sol.run(seed=seed, input_file=input_file)


class SPATIALPYDataFunction():
    """ Abstract class used to constuct the URDME data vector. """
    name = None
    def __init__(self, name=None):
        if name is not None:
            self.name = name
        if self.name is None:
            raise Exception("URDMEDataFunction must have a 'name'")

    def map(self, x):
        """ map() takes the coordinate 'x' and returns a list of doubles.
        Args:
            x: a list of 3 ints.
        Returns:
            a list of floats.
        """
        raise Exception("URDMEDataFunction.map() not implemented.")





if __name__ == '__main__':
    """ Command line interface to URDME. Execute URDME given a model file. """
