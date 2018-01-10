cmake_minimum_required(VERSION 2.8.12)
project(BrainMesh)

# CGAL and its components
find_package(CGAL QUIET COMPONENTS)
include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost REQUIRED)

find_package(OpenMesh QUIET)

# Eigen3.2 and its components
find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
include(${EIGEN3_USE_FILE})

# pybind11 
add_subdirectory(lib/pybind11)
set(PYBIND11_PYTHON_VERSION 3.6)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")

# Header files
include_directories(include)

pybind11_add_module(surface_mesh python/src/surface_mesh.cpp)
