# SurfaceVolumeMeshingToolKit

Surface Volume Meshing Toolkit (SVMTK) is a Python frontend for specific 2D and 3D mesh generation methods in the Computational Geometry Algorithms Library [(CGAL)](https://www.cgal.org/) . 

SVMTK was designed to create volume meshes of soft organic tissue surfaces, like the human brain, with the option to repair unphysical errors on the surfaces.


## Installation

First, ensure that all of the prerequisites needed by SVMTK are installed.  Installation instructions for all of the SVMTK requirements can be found in [REQUIREMENTS.md](REQUIREMENTS.md)

Next, clone the SVMT repository by 

`git clone --recursive https://github.com/SVMTK/SVMTK`

If you have already installed SVMTK you can update your installation with the following command executed in the /SVMTK directory

`git pull origin master`
`git submodule update --init`

Third, you will need to setup SVMTK if you have not already done so.  To do this you can execute 

`python3 setup.py install`

from within the SVMT directory.  If you have, instead, updated an existing SVMT installation 
you should first clean any residual files by issuing the following command

`python3 setup.py clean`

It is advised, but not required, that you use a virtual environment for python.  (see [https://docs.python.org/3/tutorial/venv.html] for instructions on using a python virtual environment)

You may see warnings as the setup process runs; these warnings arise from the CGAL source code, which is used by SVMTK, and are no cause for concern.

Finally, you can test your SVMTK installation by running 

`python3 setup.py test`

Note that, if you do not use a virtual environment, you may need to preface the above commands, for setting up and testing SVMTK, with the `sudo' keyword.  Thus, 

`sudo python3 setup.py install`
`sudo python3 setup.py clean`

If the test is successful you should see the following output 

`----------------------------------------------------------------------`

`Ran 43 tests in XXXX s`

`OK`

Also check the installation with any of the examples in `examples/`

## Docker

There is a Dockerfile in the directory `docker/`.

For more information on how to use docker, take a look at the docker tutorial:
[https://docs.docker.com/get-started/]

## Build instructions on SAGA HPC Cluster

`export PROJECT_HOME=path-to-base-install-dir`

`module load Python/3.7.4-GCCcore-8.3.0 CMake/3.13.3-GCCcore-8.2.0 GCC/8.3.0 Boost/1.71.0-GCC-8.3.0 Eigen/3.3.7 MPFR/4.0.2-GCCcore-8.3.0.lua GMP/6.1.2-GCCcore-8.3.0.lua`

`export CMAKE_PREFIX_PATH=$PROJECT_HOME/src/SVMTK/local:/cluster/software/Eigen/3.3.7:/cluster/software/Boost/1.71.0-GCC-8.3.0:/cluster/software/MPFR/4.0.2-GCCcore-8.3.0:/cluster/software/GMP/6.1.2-GCCcore-8.3.0`

`cd SVMTK/external/pybind11/`

`mkdir build`

`cd build/`

`cmake -DPYBIND11_TEST=off -DCMAKE_INSTALL_PREFIX=$PROJECT_HOME/src/SVMTK/local ..`

`make install`

`cd ../..`

`curl -sL https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.3/CGAL-5.0.3.tar.xz | tar -xJf -`

`cd CGAL-5.0.3/`

`mkdir build`

`cd build/`

`cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PROJECT_HOME/src/SVMTK/local -DWITH_Eigen3:BOOL=ON .. && make -j4`

`make install`

`cd ../../../`

`mkdir -p $PROJECT_HOME/src/SVMTK/local/lib/python3.7/site-packages/`

`export PYTHONPATH=$PROJECT_HOME/src/SVMTK/local/lib/python3.7/site-packages:$PYTHONPATH`

`python setup.py install --prefix=$PROJECT_HOME/src/SVMTK/local`

### Running SVMTK on SAGA HPC Cluster

`export PROJECT_HOME=path-to-base-install-dir`

`module load Python/3.7.4-GCCcore-8.3.0 CMake/3.12.1 GCC/8.3.0 Boost/1.71.0-GCC-8.3.0 Eigen/3.3.7 MPFR/4.0.2-GCCcore-8.3.0.lua GMP/6.1.2-GCCcore-8.3.0.lua`

`export PYTHONPATH=$PROJECT_HOME/src/SVMTK/local/lib/python3.7/site-packages:$PYTHONPATH`

`python -c "import SVMTK"`

## About

SVMTK is written by Lars Magnus Valnes and Jakob Schreiner
