#cmake_minimum_required(VERSION 2.8.12)
cmake_minimum_required(VERSION 3.5)
project(SVMTK)
set(SVMTK_VERSION_RELEASE 0)
set(SVMKT_VERSION_MAJOR "2020")
set(SVMTK_VERSION_MINOR "1")
set(SVMTK_VERSION_MICRO "0")
set(SVMTK_VERSION "${SVMTK_VERSION_MAJOR}.${SVMTK_VERSION_MINOR}.${SVMTK_VERSION_MICRO}")
set(CMAKE_BUILD_TYPE Release)


#add_subdirectory(CGAL)

#include_directories(BEFORE CGAL)
# CGAL and its components
find_package(CGAL 5 REQUIRED COMPONENTS Core GMP MPFR)
include(${CGAL_USE_FILE})
#add_subdirectory(cgal/Installation)
# Boost and its components
find_package(Boost REQUIRED  )
include_directories(${BOOST_USE_FILES})
message( ${CGAL_MAJOR_VERSION})
# Eigen3.2 and its components
find_package(Eigen3 REQUIRED)  # Requires 3.2.0 or greater
include(${EIGEN3_USE_FILE})

#find_package(GMP 4.2 REQUIRED)
#include_directories(${GOBJECT_INCLUDE_DIR})

#find_package(MPFR REQUIRED)
#include_directories(${MPFR_INCLUDE_DIR})

# pybin11
add_subdirectory(external/pybind11)
include_directories(${source_dir}/include)
# set(PYBIND11_PYTHON_VERSION 3.7)

set( CMAKE_CXX_FLAGS "-fPIC -shared ${CMAKE_CXX_FLAGS} -std=c++17 -frounding-math -lgmpxx -lgmp ")



include_directories(include)
pybind11_add_module(SVMTK python/src/SVMTK.cpp)
target_link_libraries(SVMTK PRIVATE ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES} ${OpenMesh} ${Eigen3} ${Boost})
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)

#    target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

#if (CGAL_AUTO_LINK_ENABLED)
#    target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} )
#else()
#    target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
#endif()


#enable_testing()

#ENABLE_TESTING()
#ADD_SUBDIRECTORY( tests/cpp )

#add_test (NAME python-tests
#  COMMAND ${PYTHON_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/mypyspac
