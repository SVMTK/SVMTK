cmake_minimum_required(VERSION 3.13)

# Policy CMP0077 from cmake 3.13 
cmake_policy(SET CMP0077 NEW) 

# Policy CMP0127 fromcmake 3.22 
#if(POLICY CMP0127)
#  cmake_policy(SET CMP0127 OLD)
#endif()


#set(Python3_FIND_REGISTRY "LAST")
project(SVMTK)
include(FeatureSummary)

set(SVMTK_VERSION_RELEASE 2)
set(SVMKT_VERSION_MAJOR "2")
set(SVMTK_VERSION_MINOR "0")
set(SVMTK_VERSION_MICRO "0")
set(SVMTK_VERSION "${SVMTK_VERSION_MAJOR}.${SVMTK_VERSION_MINOR}.${SVMTK_VERSION_MICRO}")
set(SVMTK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include) 

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


option(BUILD_TESTING OFF "Build SVMTK tests")
add_feature_info(BUILD_TESTING BUILD_TESTING  "Build SVMTK tests")

option(DOWNLOAD_PYBIND11 "Download Pybind11. Requires that git-repo has been cloned in recursive mode" ON)
add_feature_info(DOWNLOAD_PYBIND11 DOWNLOAD_PYBIND11 "Download Pybind11. Requires that git-repo has been cloned in recursive mode")

option(DOWNLOAD_CGAL "Download CGAL. Requires that git-repo has been cloned in recursive mode" ON)
add_feature_info(DOWNLOAD_CGAL DOWNLOAD_CGAL "Download CGAL. Requires that git-repo has been cloned in recursive mode")

if (DOWNLOAD_PYBIND11)
  set(PYBIND11_FINDPYTHON OFF) 
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7)
    message(FATAL_ERROR "GCC version must be at least 7 (for sufficient C++1z  support and to avoid some compiler bugs). You have version ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
endif()


if (${BUILD_TESTING})
   add_subdirectory(tests)
endif()

find_package (Eigen3 3.2 REQUIRED NO_MODULE)
include(${EIGEN3_USE_FILE})
include_directories(${EIGEN_INCLUDE_DIRS})
set_package_properties(Eigen3 PROPERTIES TYPE REQUIRED
  DESCRIPTION "Lightweight C++ template library for linear algebra"
  URL "http://eigen.tuxfamily.org")

if (DOWNLOAD_CGAL)
  set(CGAL_DIR external/cgal)
  find_package(CGAL REQUIRED COMPONENTS Core)
else()
    find_package(CGAL REQUIRED COMPONENTS Core)
endif()

set_package_properties(CGAL PROPERTIES TYPE REQUIRED
  DESCRIPTION "C++ template library for efficient and reliable algorithms in computational geometry"
  URL "http://www.cgal.org")
set_package_properties(MPFR PROPERTIES TYPE REQUIRED
  DESCRIPTION "C library for multiple-precision floating-point computations with correct rounding"
  URL "http://www.mpfr.org")
set_package_properties(GMP PROPERTIES TYPE REQUIRED
  DESCRIPTION "The GNU Multiple Precision Arithmetic Library"
  URL "http://gmplib.org")

find_package(Boost REQUIRED) 
include_directories(${Boost_INCLUDE_DIRS})
set_package_properties(Boost PROPERTIES TYPE REQUIRED
  DESCRIPTION "Boost C++ libraries"
  URL "http://www.boost.org")

if (DOWNLOAD_PYBIND11)
  add_subdirectory(external/pybind11)
else()
  find_package(pybind11 REQUIRED)
endif()

pybind11_add_module(SVMTK SHARED python/src/SVMTK.cpp)
set( CMAKE_CXX_FLAGS "-fPIC -flto=auto -std=gnu++1z -frounding-math -lgmpxx -lgmp -Wall -I${Boost_INCLUDE_DIRS}")   
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(SVMTK PRIVATE CGAL::CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES} ${Eigen3} ${Boost}) 

get_target_property(OUT SVMTK LINK_LIBRARIES)
message(STATUS ${OUT})
