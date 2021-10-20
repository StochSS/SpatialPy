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

import os
import shutil
import tempfile

def cleanup_tempfiles():
	'''
	Cleanup all tempfiles in spatialpy core, build, and results.
	'''
	cleanup_core_files()
	tempdir = tempfile.gettempdir()
	for file_obj in os.listdir(tempdir):
		if file_obj.startswith("spatialpy_build"):
			cleanup_build_files(build_dir=os.path.join(tempdir, file_obj))
		elif file_obj.startswith("spatialpy_result"):
			cleanup_result_files(result_dir=os.path.join(tempdir, file_obj))

def cleanup_core_files():
	'''
	Cleanup all tempfiles in spatialpy core.
	'''
	tempdir = tempfile.gettempdir()
	core_dir = os.path.join(tempdir, "spatialpy_core")
	if os.path.isdir(core_dir):
		shutil.rmtree(core_dir)
	print(f"Spatialpy core directory was removed")

def cleanup_build_files(build_dir=None):
	'''
	Cleanup all spatialpy_build directories.

	:param build_dir: Path to the build directory to be removed. (optional)
	:type build_dir: string
	'''

	if build_dir is not None:
		shutil.rmtree(build_dir)
		print(f"Build directory'{build_dir}' was removed")
	else:
		count = 0
		tempdir = tempfile.gettempdir()
		for file_obj in os.listdir(tempdir):
			if file_obj.startswith("spatialpy_build"):
				build_dir = os.path.join(tempdir, file_obj)
				shutil.rmtree(build_dir)
				count += 1
		print(f"{count} build directories were removed")

def cleanup_result_files(result_dir=None):
	'''
	Cleanup all spatialpy_result directories.

	:param result_dir: Path to the result directory to be removed. (optional)
	:type result_dir: string
	'''
	if result_dir is not None:
		shutil.rmtree(result_dir)
		print(f"Result directory '{result_dir}' was removed")
	else:
		count = 0
		tempdir = tempfile.gettempdir()
		for file_obj in os.listdir(tempdir):
			if file_obj.startswith("spatialpy_result"):
				result_dir = os.path.join(tempdir, file_obj)
				shutil.rmtree(result_dir)
				count += 1
		print(f"{count} result directories were removed")