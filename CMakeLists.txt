cmake_minimum_required(VERSION 3.5)
project(SVMTK)
set(SVMTK_VERSION_RELEASE 1)
set(SVMKT_VERSION_MAJOR "2")
set(SVMTK_VERSION_MINOR "0")
set(SVMTK_VERSION_MICRO "0")
set(SVMTK_VERSION "${SVMTK_VERSION_MAJOR}.${SVMTK_VERSION_MINOR}.${SVMTK_VERSION_MICRO}")
set(CMAKE_BUILD_TYPE Release)
message("Compiler ID ${CMAKE_CXX_COMPILER_ID}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Catch2/single_include/catch2 ) 
set(SVMTK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

SET(TEST_DIR "tests")
SET(TESTS 
    ${TEST_DIR}/test_main.cpp
    ${TEST_DIR}/test_Surface.cpp
    ${TEST_DIR}/test_Domain.cpp
)



add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable("${PROJECT_NAME}_test" ${TESTS}  )
target_link_libraries("${PROJECT_NAME}_test" PRIVATE SVMTK Catch) 
target_include_directories("${PROJECT_NAME}_test" PRIVATE  ${CATCH_INCLUDE_DIR} ${SVMTK_INCLUDE_DIR} )

find_package(Eigen3 REQUIRED) 
include(${EIGEN3_USE_FILE})

find_package(CGAL REQUIRED COMPONENTS Core )
include(${CGAL_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${BOOST_USE_FILES})

if (${CGAL_MAJOR_VERSION} GREATER 4)
  message(STATUS "CGAL version ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
  add_definitions(-DSVMTK_INSTALL_CGAL_5)
else()
  message(STATUS "CGAL version ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
endif()

add_subdirectory(external/pybind11)

set( CMAKE_CXX_FLAGS "-fPIC -std=gnu++1z -frounding-math -lgmpxx -lgmp -Wall")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
pybind11_add_module(SVMTK SHARED python/src/SVMTK.cpp)
target_link_libraries(SVMTK PRIVATE ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES} ${OpenMesh} ${Eigen3} ${Boost})
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)





