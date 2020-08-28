# SurfaceVolumeMeshingToolKit

Clone SVMTK with the flag '--recursive' to also clone the submodules Catch2, Pybind11 and CGAL.
Alternatively run `git submodule update --init`.

## Requirements

 - CGAL-5.0.3 with EIGEN >= 3.2
 - Python>=3.6

## Pybind11

Install pybind and place it in `external`

```
mkdir external && cd external
git clone https://github.com/pybind/pybind11.git --branch=v2.4.3
```

## Installation

It is highly recommended to use a virtual environment for python.

run `python3 setup.py install`

Check the installation with any of the examples in `examples/`

## Docker

There is a Dockerfile in `docker/`

For more information on how to use docker, take a look at the docker tutorial:
[https://docs.docker.com/get-started/]

### Install CGAL with Eigen3

Either download the eigen source code or install with `sudo apt-get install libeigen3-dev`

Download and install CGAL 4.13 with

`curl -sL https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13/CGAL-4.13.tar.xz | tar -xJf -`

or download the source code manually and

`cd CGAL-4.13 && cmake -DWITH_Eigen3:BOOL=ON . && make``

## Build instructions on SAGA

`export PROJECT_HOME=path-to-base-install-dir`
`module load Python/3.7.4-GCCcore-8.3.0 CMake/3.12.1 GCC/8.3.0 Boost/1.71.0-GCC-8.3.0 Eigen/3.3.7 MPFR/4.0.2-GCCcore-8.3.0.lua GMP/6.1.2-GCCcore-8.3.0.lua`
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
`mkdir -p /cluster/home/johannr/src/SVMTK/local/lib/python3.7/site-packages/`
`export PYTHONPATH=$PROJECT_HOME/src/SVMTK/local/lib/python3.7/site-packages:$PYTHONPATH`
`python setup.py install --prefix=$PWD/local`

### Running SVMTK on SAGA

`export PROJECT_HOME=path-to-base-install-dir`
`module load Python/3.7.4-GCCcore-8.3.0 CMake/3.12.1 GCC/8.3.0 Boost/1.71.0-GCC-8.3.0 Eigen/3.3.7 MPFR/4.0.2-GCCcore-8.3.0.lua GMP/6.1.2-GCCcore-8.3.0.lua`
`export PYTHONPATH=$PROJECT_HOME/src/SVMTK/local/lib/python3.7/site-packages:$PYTHONPATH`
`python -c "import SVMTK"`

## About

SVMTK is written by Lars Magnus Valnes and Jakob Schreiner
