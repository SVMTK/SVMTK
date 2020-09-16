cmake_minimum_required(VERSION 3.5)
project(SVMTK-tests)

set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)

set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_Surface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_Domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_Slice.cpp

)

set( SVMTK PARENT_SCOPE) 
set( EIGEN3_INCLUDE_DIR PARENT_SCOPE)

find_package (Eigen3 3.2 REQUIRED NO_MODULE)
include(${EIGEN3_USE_FILE})

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable(SVMTK_test ${TEST_SOURCES})
target_link_libraries(SVMTK_test PRIVATE SVMTK Catch ${CGAL_LIBRARIES} ${EIGEN3} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost}) 
target_include_directories(SVMTK_test PRIVATE ${CATCH_INCLUDE_DIR}  ${SVMTK_INCLUDE_DIR})
set_target_properties( SVMTK_test
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)
