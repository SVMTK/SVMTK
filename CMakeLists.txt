#cmake_minimum_required(VERSION 2.8.12)
cmake_minimum_required(VERSION 3.5)
project(SVMTK)
set(SVMTK_VERSION_RELEASE 1)
set(SVMKT_VERSION_MAJOR "2020")
set(SVMTK_VERSION_MINOR "0")
set(SVMTK_VERSION_MICRO "0")
set(SVMTK_VERSION "${SVMTK_VERSION_MAJOR}.${SVMTK_VERSION_MINOR}.${SVMTK_VERSION_MICRO}")
set(CMAKE_BUILD_TYPE Release)
message("Compiler ID ${CMAKE_CXX_COMPILER_ID}")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")
#configure_file(SVMTK.h.in inlcude/SVMTK.h)

SET(TEST_DIR "tests")
SET(TESTS 
    ${TEST_DIR}/test_main.cpp
#    ${TEST_DIR}/test_Slice.cpp
#    ${TEST_DIR}/test_Domain.cpp
   ${TEST_DIR}/test_Surface.cpp
)

#set(CATCH_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/Catch2/include")
#set(SVMTK_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
#add_library(Catch INTERFACE)
#target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})
#message(${CATCH_INCLUDE_DIR})

# Generate a test executable

include_directories(external/Catch2/include)
add_executable("${PROJECT_NAME}_test" ${TESTS} )
set_source_files_properties(${TESTS} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}" )
target_link_libraries(${PROJECT_NAME}_test PRIVATE SVMTK)
#target_include_directories(${PROJECT_NAME}_test PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>/include $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)


find_package(Eigen3 REQUIRED) 
include(${EIGEN3_USE_FILE})

find_package(CGAL REQUIRED COMPONENTS Core )
include(${CGAL_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${BOOST_USE_FILES})

if (${CGAL_MAJOR_VERSION} GREATER 4)
  message(STATUS "CGAL version ${CGAL_MAJOR_VERSION}")
  add_definitions(-DSVMTK_INSTALL_CGAL_5)
else()
  message(STATUS "CGAL version ${CGAL_MAJOR_VERSION}")
endif()

add_subdirectory(external/pybind11)

set( CMAKE_CXX_FLAGS "-fPIC -std=c++17 -frounding-math -lgmpxx -lgmp")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
pybind11_add_module(SVMTK SHARED python/src/SVMTK.cpp)
target_link_libraries(SVMTK PRIVATE ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES} ${OpenMesh} ${Eigen3} ${Boost}) #INTERFACE ?? OPENMESH??
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)






#option(ENABLE_TESTS "Enable testing" On)
#message(STATUS "Tests are ${ENABLE_TESTS}")

#enable_testing()
#add_subdirectory(tests EXCLUDE_FROM_ALL)


#  add_subdirectory(tests)
#  message(STATUS "Testing Enabled")




