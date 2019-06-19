from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.easy_install import easy_install
import os

SETUP_DIR = os.path.dirname(os.path.abspath(__file__))


def stoch_path(command_subclass):
    """
    A decorator for classes subclassing one of the setuptools commands.
    It modifies the run() method.
    """
    orig_run = command_subclass.run

    def modified_run(self):
        success = False
        orig_run(self)

    command_subclass.run = modified_run
    return command_subclass


# update all install classes with our new class
@stoch_path
class develop_new(develop):
    pass


@stoch_path
class install_new(install):
    pass


@stoch_path
class bdist_egg_new(bdist_egg):
    pass


@stoch_path
class easy_install_new(easy_install):
    pass


with open('README.md', 'r') as fh:
    full_description = fh.read()



setup(name="spatialpy",
      version="0.1.0",
      packages=['spatialpy'],

      include_package_data = True,
      #package_data={'spatialpy':['data/*.c','data/three.js_templates/js/*','data/three.js_templates/*.html','spatialpy/AUTHORS','spatialpy/LICENCE','spatialpy/bin/*','spatialpy/build/*','spatialpy/include/*','spatialpy/src/*.c','spatialpy/src/nsm/*']},
      
       description='Python Interface for Spatial Stochastic Biochemical Simulations', 

      install_requires = ["numpy",
                          "scipy"],
      
      author="Brian Drawert, Evie Hilton",
      author_email="briandrawert@gmail.com",
      license = "GPL",
      keywords = "spatialpy, spatial stochastic simulation, RDME",
      url="https://spatialpy.github.io/SpatialPy/",
      download_url="https://github.com/SpatialPy/SpatialPy/tarball/master/",


      cmdclass={'bdist_egg': bdist_egg_new,
                'install': install_new,
                'develop': develop_new,
                'easy_install': easy_install_new},
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research'
      ],

      )
      
