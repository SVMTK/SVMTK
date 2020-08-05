#cmake_minimum_required(VERSION 2.8.12)
cmake_minimum_required(VERSION 3.5)
project(SVMTK)
set(SVMTK_VERSION_RELEASE 1)
set(SVMKT_VERSION_MAJOR "2020")
set(SVMTK_VERSION_MINOR "0")
set(SVMTK_VERSION_MICRO "0")
set(SVMTK_VERSION "${SVMTK_VERSION_MAJOR}.${SVMTK_VERSION_MINOR}.${SVMTK_VERSION_MICRO}")
set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")




find_package(Eigen3 REQUIRED)  # Requires 3.2.0 or greater
include(${EIGEN3_USE_FILE})

find_package(CGAL REQUIRED COMPONENTS Core GMP MPFR )
include(${CGAL_USE_FILE})

find_package(Boost REQUIRED  )
include_directories(${BOOST_USE_FILES})


# 
if (${CGAL_MAJOR_VERSION} GREATER 4)
  message(STATUS "CGAL version ${CGAL_MAJOR_VERSION}")
  add_definitions(-DSVMTK_INSTALL_CGAL_5)
else()
  message(STATUS "CGAL version ${CGAL_MAJOR_VERSION}")
endif()


add_subdirectory(external/pybind11)
include_directories(${source_dir}/include)

set( CMAKE_CXX_FLAGS "-fPIC  -shared ${CMAKE_CXX_FLAGS} -std=c++17 -frounding-math -lgmpxx -lgmp ")

include_directories(include)
pybind11_add_module(SVMTK python/src/SVMTK.cpp)
target_link_libraries(SVMTK INTERFACE ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES} ${OpenMesh} ${Eigen3} ${Boost})
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)


#SET(TEST_DIR "tests")
#SET(TESTS
    "${TEST_DIR}/test_main.cpp"
    "${TEST_DIR}/test_Domain.cpp")

#include_directories(lib/catch/include)
#add_executable("${PROJECT_NAME}_test" ${TESTS})
#option(ENABLE_TESTS "Enable testing" On)
#message(STATUS "Tests are ${ENABLE_TESTS}")

#enable_testing()
#add_subdirectory(tests EXCLUDE_FROM_ALL)


#  add_subdirectory(tests)
#  message(STATUS "Testing Enabled")




