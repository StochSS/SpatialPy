from setuptools import setup

setup(name="spatialpy",
      version="1.1.1",
      packages=['spatialpy'],
      
      #include_package_data = True,
      package_data={'spatialpy':['data/*.c','data/three.js_templates/js/*','data/three.js_templates/*.html','spatialpy/AUTHORS','spatialpy/LICENCE','spatialpy/bin/*','spatialpy/build/*','spatialpy/include/*','spatialpy/src/*.c','spatialpy/src/nsm/*']},
      
      install_requires = ["numpy",
                          "scipy",
                          "h5py"],
      
      author="Brian Drawert",
      author_email=["briandrawert@gmail.com"],
      license = "GPL",
      keywords = "spatialpy, spatial stochastic simulation, RDME",
     
      
      )
      
