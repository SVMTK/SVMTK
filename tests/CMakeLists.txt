cmake_minimum_required(VERSION 3.5)
project(SVMTK-tests)


set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)

set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_Surface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_Domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_Slice.cpp

)
message("AFSFAFAGEAGEFA") 
find_package (Eigen3 3.2 REQUIRED)
include(${EIGEN3_USE_FILE})

add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable(SVMTK_test ${TEST_SOURCES})

target_link_libraries(SVMTK_test PRIVATE CGAL::CGAL SVMTK Catch ) 
target_include_directories(SVMTK_test PRIVATE CGAL::Eigen3_support ${CATCH_INCLUDE_DIR} ${SVMTK_INCLUDE_DIR} )

set_target_properties(SVMTK_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)



