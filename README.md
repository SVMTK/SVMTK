# SurfaceVolumeMeshingToolKit

## Requirements

 - CGAL-4.13 with EIGEN >= 3.2
 - Python3.6

## Pybind11

If you did not clone this project using `git clone --recursive`, pull pybind11 by
`git submodule update --recursive`

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
