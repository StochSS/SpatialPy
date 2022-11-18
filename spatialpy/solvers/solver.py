# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2022 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import shutil
import signal
import subprocess
import tempfile
import threading
import time
import getpass
import re

import numpy

from spatialpy.core.spatialpyerror import ModelError, SimulationError, SimulationTimeout

def _read_from_stdout(stdout ,verbose=True):
    try:
        while True:
            line = stdout.readline()
            if line != b'':
                if verbose:
                    print(line.decode(),end='')
            else:
                #got empty line, ending
                return
    except Exception as err:
        print(f"_read_from_stdout(): {err}")


class Solver:
    """
    SpatialPy solver object.

    :param model: Target model of solver simulation.
    :type model: spatialpy.core.model.Model

    :param debug_level: Target level of debugging.
    :type debug_level: int
    """
    def __init__(self, model, debug_level=0):
        from spatialpy.core.model import Model # pylint: disable=import-outside-toplevel
        if not (isinstance(model, Model) or type(model).__name__ == 'Model'):
            raise SimulationError("Model must be of type spatialpy.Model.")
        if not issubclass(self.__class__, Solver):
            raise SimulationError("Solver classes must be a subclass of spatialpy.Solver.")

        self.model = model
        self.is_compiled = False
        self.debug_level = debug_level
        self.model_name = self.model.name
        self.build_dir = None
        self.propfilename = None
        self.prop_file_name = None
        self.executable_name = 'ssa_sdpd'
        self.h = None  # basis function width

        self.spatialpy_root = os.path.dirname(
            os.path.abspath(__file__)) + "/c_base/ssa_sdpd-c-simulation-engine"
        self.spatialpy_rootdir =  self.spatialpy_root.replace(" ","\\ ")
        self.spatialpy_rootinc =  self.spatialpy_root.replace(" ","\\\\ ")
        self.spatialpy_rootparam =  self.spatialpy_root.replace(" ","?")

        tmpdir = tempfile.gettempdir()
        self.core_dir = os.path.join(os.path.join(tmpdir, 'spatialpy_core'), getpass.getuser())

        if not os.path.isdir(os.path.join(tmpdir, 'spatialpy_core')):
            os.mkdir(os.path.join(tmpdir, 'spatialpy_core'))
        if not os.path.isdir(self.core_dir):
            os.mkdir(self.core_dir)

        self.debugger_url = None
        self.debugger_process = None

    def __del__(self):
        try:
            if self.build_dir is not None:
                try:
                    shutil.rmtree(self.build_dir, ignore_errors=True)
                except OSError:
                    print(f"Could not delete '{self.build_dir}'")
        except Exception:
            pass

    def __create_propensity_file(self, stoich_matrix, dep_graph, file_name=None):
        # process initial conditions here
        nspecies = self.model.u0.shape[0]
        ncells = self.model.u0.shape[1]

        num_types = len(self.model.domain.listOfTypeIDs)
        if self.model.enable_pde:
            num_chem_species = len(self.model.listOfSpecies)
            num_chem_rxns = len(self.model.listOfReactions)
        else:
            num_chem_species = 0
            num_chem_rxns = 0
        if self.model.enable_rdme:
            num_stoch_species = len(self.model.listOfSpecies)
            num_stoch_rxns = len(self.model.listOfReactions)
        else:
            num_stoch_species = 0
            num_stoch_rxns = 0
        num_data_fn = len(self.model.listOfDataFunctions)

        # Reactions
        funcs, funcinits = self.__get_reaction_prop()
        deterministic_chem_rxn_functions, deterministic_chem_rxn_function_init = self.__get_chem_rxn_prop()

        # End of pyurdme replacements
        # SSA-SDPD values here
        if(self.model.enable_rdme and len(self.model.listOfSpecies) > 0):
            init_rdme = "initialize_rdme(system, input_irN, input_jcN, input_prN, input_irG, input_jcG, input_u0);"
        else:
            init_rdme=''

        init_bc = ""
        for bound_cond in self.model.listOfBoundaryConditions:
            init_bc += bound_cond.expression()

        replacements = {
            "__NUMBER_OF_REACTIONS__": str(self.model.get_num_reactions()),
            "__NUMBER_OF_SPECIES__": str(self.model.get_num_species()),
            "__NUMBER_OF_VOXELS__": str(self.model.domain.get_num_voxels()),
            "__DEFINE_PARAMETERS__": self.__get_param_defs(),
            "__DEFINE_REACTIONS__": funcs,
            "__DEFINE_PROPFUNS__": funcinits,
            "__DEFINE_CHEM_FUNS__": deterministic_chem_rxn_functions,
            "__DEFINE_CHEM_FUN_INITS__": deterministic_chem_rxn_function_init,
            "__INIT_PARTICLES__": self.__get_particle_inits(num_chem_species),
            "__DATA_FUNCTION_ASSIGN__": self.__get_data_fn_assign(ncells),
            "__INPUT_CONSTANTS__": self.__get_input_constants(nspecies, ncells, stoich_matrix, dep_graph),
            "__SYSTEM_CONFIG__": self.__get_system_config(num_types, num_chem_species, num_chem_rxns,
                                                          num_stoch_species, num_stoch_rxns, num_data_fn),
            "__INIT_RDME__": init_rdme,
            "__BOUNDARY_CONDITIONS__": init_bc,
            "__DEFINE_GET_NEXT_OUTPUT__": self.__get_next_output()

        }

        propfilestr = self.__get_propfile_str(replacements)

        #### Write the data to the file ####
        self.__write_prop_file(file_name, propfilestr)

    def __get_chem_rxn_prop(self):
        funheader = "double det__NAME__(const double *x, double t, const double vol, const double *data_fn, int sd)"

        deterministic_chem_rxn_functions = ""
        deterministic_chem_rxn_function_init = ""
        for i, (rname, reac) in enumerate(self.model.listOfReactions.items()):
            ode_propensity_function = self.model.expr.getexpr_cpp(reac.ode_propensity_function)
            func = funheader.replace("__NAME__", rname)
            func += "\n{\n"
            if reac.restrict_to is None or (isinstance(reac.restrict_to, list) and len(reac.restrict_to) == 0):
                func += f"return {ode_propensity_function};"
            else:
                func += "if("
                if isinstance(reac.restrict_to, list) and len(reac.restrict_to) > 0:
                    conds = []
                    for type_id in reac.restrict_to:
                        conds.append(f"sd == {type_id}")
                    func += "||".join(conds)
                else:
                    errmsg = "When restricting reaction to types, you must specify either a list or an int"
                    raise SimulationError(errmsg)
                func += "){\n"
                func += f"return {ode_propensity_function};"
                func += "\n}else{"
                func += "\n\treturn 0.0;}"

            func += "\n}"
            deterministic_chem_rxn_functions += f"{func}\n\n"
            deterministic_chem_rxn_function_init += f"    ptr[{i}] = (ChemRxnFun) det{rname};\n"

        return deterministic_chem_rxn_functions, deterministic_chem_rxn_function_init

    def __get_data_fn_assign(self, ncells):
        data_fn_assign = ""
        if len(self.model.listOfSpecies) > 0:
            for ndf, _ in enumerate(self.model.listOfDataFunctions):
                data_fn_assign += f"this_particle->data_fn[{ndf}] = input_data_fn[{ndf}*{ncells}+id];"

        return data_fn_assign

    def __get_input_constant(self, header, data, to_int=False):
        outstr = f"{header}[{len(data)}] = "
        outstr += "{"
        outstr += ",".join([str(int(val) if to_int else val) for val in data])
        outstr += "};\n"

        return outstr

    def __get_input_constants(self, nspecies, ncells, stoich_matrix, dep_graph):
        input_constants = ""

        outstr = f"unsigned int input_u0[{nspecies * ncells}] = "
        outstr += "{"
        if len(self.model.listOfSpecies) > 0:
            for i in range(ncells):
                for j in range(nspecies):
                    if i + j > 0:
                        outstr += ','
                    outstr += str(int(self.model.u0[j, i]))
        outstr += "};\n"
        input_constants += outstr

        if len(self.model.listOfSpecies) > 0:
            if min(stoich_matrix.shape) > 0:
                dense_matrix = stoich_matrix.todense() # this will not work if Nrxn or Nspecies is zero
                outstr = f"static int input_N_dense[{dense_matrix.shape[0] * dense_matrix.shape[1]}] = "
                outstr += "{"
                for i in range(dense_matrix.shape[0]):
                    for j in range(dense_matrix.shape[1]):
                        if j + i > 0:
                            outstr += ','
                        outstr += f"{int(dense_matrix[i, j])}"
                outstr += "};\n"
                input_constants += outstr

                input_constants += self.__get_input_constant("static size_t input_irN", stoich_matrix.indices)
                input_constants += self.__get_input_constant("static size_t input_jcN", stoich_matrix.indptr)
                input_constants += self.__get_input_constant("static int input_prN", stoich_matrix.data, to_int=True)

                if len(self.model.listOfDataFunctions) > 0:
                    outstr = f"static double input_data_fn[{len(self.model.listOfDataFunctions) * ncells}] = "
                    outstr += "{"
                    coords = self.model.domain.coordinates()
                    for ndf, data_fn in enumerate(self.model.listOfDataFunctions.values()):
                        for i in range(ncells):
                            if i > 0 and ndf == 0:
                                outstr += ','
                            outstr += str(data_fn.map([coords[i, 0], coords[i, 1], coords[i, 2]]))
                    outstr += "};\n"
                    input_constants += outstr
            else:
                input_constants += "static int input_N_dense[0] = {};\n"
                input_constants += "static size_t input_irN[0] = {};\n"
                input_constants += "static size_t input_jcN[0] = {};\n"
                input_constants += "static int input_prN[0] = {};\n"

            input_constants += self.__get_input_constant("static size_t input_irG", dep_graph.indices)
            input_constants += self.__get_input_constant("static size_t input_jcG", dep_graph.indptr)

        if len(self.model.listOfSpecies) > 0:
            outstr = "const char* const input_species_names[] = {"
            outstr += ",".join(f'"{sname}"' for sname in self.model.listOfSpecies)
            outstr += ", 0};\n"
            input_constants += outstr

        num_types = len(self.model.domain.listOfTypeIDs)
        outstr = f"const int input_num_subdomain = {num_types-1};\n" # the backend ignores type_id=0 (a.k.a. unassigned)
        input_constants += outstr

        outstr = f"const double input_subdomain_diffusion_matrix[{len(self.model.listOfSpecies) * (num_types-1)}] = "
        outvec = []
        for i, species in enumerate(self.model.listOfSpecies.values()):
            for j, type_id in enumerate(self.model.domain.typeNdxMapping.keys()):
                if j==0: continue  # do not process type_id=0 in diffusion
                try:
                    if species not in self.model.listOfDiffusionRestrictions or \
                       type_id in self.model.listOfDiffusionRestrictions[species]:
                        outvec.append(f"{species.diffusion_coefficient}")
                    else:
                        outvec.append("0.0")
                except KeyError as err:
                    print(f"error: {err}")
                    print(self.model.listOfDiffusionRestrictions)
                    raise SimulationError(f"error: {self.model.listOfDiffusionRestrictions}") from err

        outstr += "{"+",".join(outvec)+"};\n"
        input_constants += outstr

        return input_constants

    def __get_next_output(self):
        output_step = "unsigned int get_next_output(ParticleSystem* system)\n{\n"
        output_step += "static int index = 0;\n"
        output_step += "const std::vector<unsigned int> output_steps = {"
        output_step += f"{', '.join(self.model.tspan.output_steps.astype(str).tolist())}"
        output_step += "};\nunsigned int next_step = output_steps[index];\n"
        output_step += "index++;\n"
        output_step += "return next_step;\n}\n"

        return output_step

    def __get_param_defs(self):
        sanitized_parameters = self.model.sanitized_parameter_names()
        parameters = ""
        for pname in self.model.listOfParameters:
            param = sanitized_parameters[pname]
            parameters += f"const double {param} = {self.model.listOfParameters[pname].value};\n"

        for name, ndx in self.model.domain.typeNdxMapping.items():
            parameters += f"const size_t {name} = {ndx};\n"
        return parameters

    def __get_particle_inits(self, num_chem_species):
        init_particles = ""
        if self.model.domain.type_id is None:
            self.model.domain.type_id = ["type_1"] * self.model.domain.get_num_voxels()
        for i, type_id in enumerate(self.model.domain.type_id):
            if "UnAssigned" in type_id:
                errmsg = "Not all particles have been defined in a type. Mass and other properties must be defined"
                raise SimulationError(errmsg)
            x = self.model.domain.coordinates()[i, 0]
            y = self.model.domain.coordinates()[i, 1]
            z = self.model.domain.coordinates()[i, 2]
            nu = self.model.domain.nu[i]
            mass = self.model.domain.mass[i]
            c = self.model.domain.c[i]
            rho = self.model.domain.rho[i]
            fixed = int(self.model.domain.fixed[i])
            init_particles += "init_create_particle(sys,id++,"
            init_particles += f"{x},{y},{z},{type_id},{nu},{mass},{c},{rho},{fixed},{num_chem_species});\n"

        return init_particles

    def __get_propfile_str(self, replacements):
        temp_path = os.path.abspath(os.path.dirname(__file__))
        temp_path += '/c_base/ssa_sdpd-c-simulation-engine/propensity_file_template.cpp'
        with open(temp_path, 'r', encoding="utf-8") as template:
            propfilestr = template.read()

        for holder, data in replacements.items():
            propfilestr = propfilestr.replace(holder, data)

        return propfilestr

    def __get_reaction_prop(self):
        funheader = "double __NAME__(const int *x, double t, const double vol, const double *data_fn, int sd)"

        funcs = ""
        funcinits = ""
        for i, (rname, reac) in enumerate(self.model.listOfReactions.items()):
            propensity_function = self.model.expr.getexpr_cpp(reac.propensity_function)
            func = funheader.replace("__NAME__", rname)
            func +=  "\n{\n"
            if reac.restrict_to is None or (isinstance(reac.restrict_to, list) and len(reac.restrict_to) == 0):
                func += f"return {propensity_function};"
            else:
                func += "if("
                if isinstance(reac.restrict_to, list) and len(reac.restrict_to) > 0:
                    conds = []
                    for type_id in reac.restrict_to:
                        conds.append(f"sd == {type_id}")
                    func += "||".join(conds)
                else:
                    errmsg = "When restricting reaction to types, you must specify either a list or an int"
                    raise SimulationError(errmsg)
                func += "){\n"
                func += f"return {propensity_function};"
                func += "\n}else{"
                func += "\n\treturn 0.0;}"
            func += "\n}"
            funcs += f"{func}\n\n"
            funcinits += f"    ptr[{i}] = (PropensityFun) {rname};\n"

        return funcs, funcinits

    def __get_system_config(self, num_types, num_chem_species, num_chem_rxns,
                            num_stoch_species, num_stoch_rxns, num_data_fn):
        system_config = f"debug_flag = {self.debug_level};\n"
        system_config += "ParticleSystem *system = new ParticleSystem("
        system_config += f"{num_types-1},{num_chem_species},{num_chem_rxns},"
        system_config += f"{num_stoch_species},{num_stoch_rxns},{num_data_fn});\n"
        system_config += f"system->static_domain = {int(self.model.staticDomain)};\n"
        if len(self.model.listOfSpecies) > 0:
            system_config += "system->subdomain_diffusion_matrix = input_subdomain_diffusion_matrix;\n"
            system_config += "system->stoichiometric_matrix = input_N_dense;\n"
            system_config += "system->chem_rxn_rhs_functions = ALLOC_ChemRxnFun();\n"
            system_config += "system->stoch_rxn_propensity_functions = ALLOC_propensities();\n"
            system_config += "system->species_names = input_species_names;\n"

        system_config += f"system->dt = {self.model.tspan.timestep_size};\n"
        system_config += f"system->nt = {self.model.tspan.num_timesteps};\n"
        if self.h is None:
            self.h = self.model.domain.find_h()
        if self.h == 0.0:
            raise ModelError('h (basis function width) can not be zero.')
        system_config += f"system->h = {self.h};\n"
        system_config += f"system->rho0 = {self.model.domain.rho0};\n"
        system_config += f"system->c0 = {self.model.domain.c0};\n"
        system_config += f"system->P0 = {self.model.domain.P0};\n"
        #// bounding box
        system_config += f"system->xlo = {self.model.domain.xlim[0]};\n"
        system_config += f"system->xhi = {self.model.domain.xlim[1]};\n"
        system_config += f"system->ylo = {self.model.domain.ylim[0]};\n"
        system_config += f"system->yhi = {self.model.domain.ylim[1]};\n"
        system_config += f"system->zlo = {self.model.domain.zlim[0]};\n"
        system_config += f"system->zhi = {self.model.domain.zlim[1]};\n"

        if not numpy.count_nonzero(self.model.domain.vertices[:,1]):
            self.model.domain.dimensions = 1
        elif not numpy.count_nonzero(self.model.domain.vertices[:,2]):
            self.model.domain.dimensions = 2
        else:
            self.model.domain.dimensions = 3
        system_config += f"system->dimension = {self.model.domain.dimensions};\n"

        if self.model.domain.gravity is not None:
            for i, val in enumerate(self.model.domain.gravity):
                system_config += f"system->gravity[{i}] = {val};\n"

        return system_config

    def __read_profile_info(self, result):
        profile_data_path = os.path.join(result.result_dir, 'gmon.out')
        exe_path = os.path.join(self.build_dir, self.executable_name)
        cmd = f'gprof {exe_path} {profile_data_path}'
        print(cmd)
        with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) as proc:
            stdout, _ = proc.communicate()
            proc.wait()
            print(f'Gprof report for {result.result_dir}')
            print(stdout.decode('utf-8'))

    def __run_debugger(self):
        self.debugger_url = 'http://127.0.0.1:5000'
        if not hasattr(self, 'debugger_process'):
            with subprocess.Popen('gdbgui -r', shell=True) as debugger_process:
                self.debugger_process = debugger_process
        print(f'Your debugger is running at {self.debugger_url}')

    def __write_prop_file(self, file_name, propfilestr):
        with open(file_name, "w", encoding="utf-8") as propfile:
            propfile.write(propfilestr)

    def compile(self, debug=False, profile=False):
        """
        Compile the model.

        :param debug: If True, will print additional build debugging
        :type debug: bool

        :param profile: If True, will print additional profiling information
        :type profile: bool

        :raises SimulationError: Failed to compile
        """
        stoich_matrix, dep_graph = self.model.compile_prep()

        # Create a unique directory each time call to compile.
        self.build_dir = tempfile.mkdtemp(
            prefix='spatialpy_build_', dir=os.environ.get('SPATIALPY_TMPDIR'))
        # Write the propensity file
        # Match except word characters \w = ([a-zA-Z0-9_]) and \_ = _ replace with ''
        propfilename = re.sub(r'[^\w\_]', '', self.model_name)
        self.propfilename = f"{propfilename}_generated_model"
        self.prop_file_name = os.path.join(self.build_dir, f'{self.propfilename}.c')

        if self.debug_level >= 1:
            print(f"Compiling Solver.  Build dir: {self.build_dir}")

        if self.debug_level >= 1:
            print(f"Creating propensity file {self.prop_file_name}")
        self.__create_propensity_file(stoich_matrix, dep_graph, file_name=self.prop_file_name)

        # Build the solver
        for make_exe_location in ["make", "mingw32-make"]:
            make_exe = shutil.which("make")
            if make_exe is not None:
                break
        if make_exe is None:
            raise SimulationError("Make executable could not be found")
        makefile = self.spatialpy_rootdir+'/build/Makefile'
        makefile_ann = self.spatialpy_rootdir+'/external/ANN/src/Makefile.spatialpy'

        cmd_ann = [
            make_exe, '-d', '-C', self.core_dir, '-f', makefile_ann,
            f'ROOTINC={self.spatialpy_rootinc}',
        ]
        cmd_core = [
            make_exe, '-d', '-C', self.core_dir, 'CORE', '-f',  makefile,
            f'ROOT={self.spatialpy_rootparam}',
            f'ROOTINC={self.spatialpy_rootinc}',
            f'BUILD={self.core_dir}',
        ]
        cmd_build = [
            make_exe, '-d', '-C', self.build_dir, '-I', self.core_dir, '-f', makefile,
            'ROOT=' + self.spatialpy_rootparam,
            'ROOTINC=' + self.spatialpy_rootinc,
            'COREDIR=' + self.core_dir,
            'MODEL=' + self.prop_file_name, 'BUILD='+self.build_dir
        ]
        if profile:
            cmd_build.append('GPROFFLAG=-pg')
        if profile or debug:
            cmd_build.append('GDB_FLAG=-g')
        if self.debug_level > 1:
            cmd = " && ".join([*cmd_ann, *cmd_core, *cmd_build])
            print(f"cmd: {cmd}\n")
        try:
            for cmd_target in [cmd_ann, cmd_core, cmd_build]:
                result = subprocess.check_output(cmd_target, stderr=subprocess.STDOUT)
                if self.debug_level > 1:
                    print(result.stdout.decode("utf-8"))
        except subprocess.CalledProcessError as err:
            try:
                print(err.stdout.decode("utf-8"))
            except Exception:
                pass
            raise SimulationError(f"Compilation of solver failed, return_code={result.return_code}")
        except OSError as err:
            print(f"Error, execution of compilation raised an exception: {err}")
            print(f"cmd = {cmd}")
            raise SimulationError("Compilation of solver failed") from err

        self.is_compiled = True


    def run(self, number_of_trajectories=1, seed=None, timeout=None,
                number_of_threads=None, debug=False, profile=False, verbose=True):
        """
        Run one simulation of the model.

        :param number_of_trajectories: How many trajectories should be simulated.
        :type number_of_trajectories: int

        :param seed: the random number seed (incremented by one for multiple runs).
        :type seed: int

        :param timeout: maximum number of seconds the solver can run.
        :type timeout: int

        :param number_of_threads: the number threads the solver will use.
        :type number_of_threads: int

        :param debug: start a gdbgui debugger (also compiles with debug symbols if compilation hasn't happened)
        :type debug: bool

        :param profile: Output gprof profiling data if available
        :type profile: bool

        :param verbose: If true, prints addtional data to console
        :type verbose: bool

        :returns: A SpatialPy Result object containing spatial and time series data from simulation.
        :rtype: spatialpy.Result.Result | list(spatialpy.Result.Result)

        :raises SimulationTimeout: Simulation exceeded timeout.
        :raises SimulationError: Simulation execution failed.
        """
        from spatialpy.core.result import Result # pylint: disable=import-outside-toplevel
        if number_of_trajectories > 1:
            result_list = []
        # Check if compiled, call compile() if not.
        if not self.is_compiled:
            self.compile(debug=debug, profile=profile)

        # Execute the solver
        for run_ndx in range(number_of_trajectories):
            outfile = tempfile.mkdtemp(
                prefix='spatialpy_result_', dir=os.environ.get('SPATIALPY_TMPDIR'))
            result = Result(self.model, outfile)
            if self.debug_level >= 1:
                print(f"Running simulation. Result dir: {outfile}")
            solver_cmd = os.path.join(self.build_dir, self.executable_name)

            if number_of_threads is not None:
                solver_cmd += " -t " + str(number_of_threads)

            if seed is not None:
                solver_cmd += " -s "+str(seed+run_ndx)

            if self.debug_level >= 1:
                print(f'cmd: {solver_cmd}')

            start = time.monotonic()
            return_code = None
            try:
                with subprocess.Popen(solver_cmd, cwd=outfile, shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        start_new_session=True) as process:
                    try:
                        # start thread to read process stdout to stdout
                        thread = threading.Thread(target=_read_from_stdout, args=(process.stdout,verbose))
                        thread.start()
                        if timeout is not None:
                            return_code = process.wait(timeout=timeout)
                        else:
                            return_code = process.wait()
                        thread.join()
                    except KeyboardInterrupt:
                        # send signal to the process group
                        os.killpg(process.pid, signal.SIGINT)
                        print('Terminated by user after seconds: {:.2f}'.format(time.monotonic() - start))
                    except subprocess.TimeoutExpired as err:
                        result.timeout = True
                        # send signal to the process group
                        os.killpg(process.pid, signal.SIGINT)
                        raise SimulationTimeout("SpatialPy solver timeout exceded.") from err
            except OSError as err:
                print(f"Error, execution of solver raised an exception: {err}")
                print(f"cmd = {solver_cmd}")

            if self.debug_level >= 1:  # output time
                print('Elapsed seconds: {:.2f}'.format(time.monotonic() - start))

            if return_code is not None and return_code != 0:
                print(f"solver_cmd = {solver_cmd}")
                raise SimulationError(f"Solver execution failed, return code = {return_code}")

            result.success = True
            if profile:
                self.__read_profile_info(result)
            if number_of_trajectories > 1:
                result_list.append(result)
            else:
                return result

        return result_list
