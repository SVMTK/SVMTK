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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include(CGAL_CreateSingleSourceCGALProgram)
create_single_source_cgal_program("CGALSurface.cpp")
#create_single_source_cgal_program("CGALMeshCreator.cpp")
#target_link_libraries(CGALSurface ${Boost_LIBRARIES})
#target_link_libraries( ${BOOST_LIBRARIES} )
