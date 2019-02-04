cmake_minimum_required(VERSION 2.8.12)
project(BrainMesh)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE RELEASE)

# CGAL and its components
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost REQUIRED)
include_directories(${BOOST_USE_FILES})

# Eigen3.2 and its components
if(EXISTS "$ENV{EIGEN3_INCLUDE_DIR}")
    include("$ENV{EIGEN3_INCLUDE_DIR}")
else()
    find_package(Eigen3 3.2.0 REQUIRED)  # Requires 3.2.0 or greater     # REQUIRED
    include(${EIGEN3_USE_FILE})
endif()

# pybin11
add_subdirectory(external/pybind11)

set(SOURCE_DIR "include")
include_directories(${SOURCE_DIR})

set(PYBIND11_PYTHON_VERSION 3.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
pybind11_add_module(brainmesh SHARED python/src/brainmesh.cpp)
