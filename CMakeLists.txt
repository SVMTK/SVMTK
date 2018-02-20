cmake_minimum_required(VERSION 2.8.12)
project(BrainMesh)

# CGAL and its components
find_package(CGAL QUIET COMPONENTS)
include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost REQUIRED)
find_package(OpenMesh QUIET)    # TODO: What is this for and do we need it?

# Eigen3.2 and its components
find_package(Eigen3 3.2.0)  # Requires 3.2.0 or greater
include(${EIGEN3_USE_FILE})

# pybin11
add_subdirectory(external/pybind11)
include_directories(${source_dir}/include)
set(PYBIND11_PYTHON_VERSION 3.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# headers
include_directories(include)
pybind11_add_module(brainmesh python/src/brainmesh.cpp)
