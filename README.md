# SurfaceVolumeMeshingToolKit

## Requirements

 - CGAL-5.0.2 with EIGEN >= 3.2
 - Python3.6

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

Download and install CGAL 5.0.2 with

`curl -sL https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.2/CGAL-5.0.2.tar.xz | tar -xJf -`

or download the source code manually and

`cd CGAL-5.0.2 && cmake -DWITH_Eigen3:BOOL=ON . && make``


## About

SVMTK is written by Lars Magnus Valnes and Jakob Schreiner
