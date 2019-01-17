cmake_minimum_required(VERSION 2.8.12)
project(BrainMesh)

# CGAL and its components
find_package(CGAL 4.13 REQUIRED)
include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost REQUIRED)
include_directories(${BOOST_USE_FILES})

# Eigen3.2 and its components
find_package(Eigen3 3.2.0 REQUIRED)  # Requires 3.2.0 or greater
include(${EIGEN3_USE_FILE})

find_package(GMP REQUIRED)
include_directories(${GOBJECT_INCLUDE_DIR})

find_package(MPFR REQUIRED)
include_directories(${MPFR_INCLUDE_DIR})

# pybin11
add_subdirectory(external/pybind11)
include_directories(${source_dir}/include)
# set(PYBIND11_PYTHON_VERSION 3.7)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -frounding-math -fPIC")

add_library(CGALMeshCreator include/CGALMeshCreator.cpp)
add_library(CGALSurface include/CGALSurface.cpp)
set(BMLIB CGALMeshCreator CGALSurface)

include_directories(include)
pybind11_add_module(brainmesh python/src/brainmesh.cpp)
target_link_libraries(brainmesh PRIVATE ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GOBJECT_LIBRARIES} ${OpenMesh} ${Eigen3} ${Boost} ${MPFR} ${BMLIB})

#    target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

#if (CGAL_AUTO_LINK_ENABLED)
#    target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} )
#else()
#    target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
#endif()
