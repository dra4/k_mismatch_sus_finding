cmake_minimum_required(VERSION 2.8.8)

set(PROJECT_NAME_STR k_mismatch_sus_finding)
PROJECT(${PROJECT_NAME_STR} C CXX)
#TODO CHANGE THIS OPTION
set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_BUILD_TYPE Release)

find_package(Threads REQUIRED)

# SDSL directory
set(SDSL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/ext/sdsl-lite")

# flags
set(CMAKE_CXX_FLAGS                "-Wall -Wextra -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# include the sdsl files and googletest library
add_subdirectory("${PROJECT_SOURCE_DIR}/ext/libdivsufsort")

# add libdivsufsort includes path
include_directories("${CMAKE_BINARY_DIR}/ext/libdivsufsort/include")
include_directories("${PROJECT_SOURCE_DIR}/src/sdsl")

# Compile project
ADD_EXECUTABLE(${PROJECT_NAME_STR}.x
               src/sdsl/bits.cpp src/sdsl/construct_sa.cpp
               src/ReadsDB.cpp src/AppConfig.cpp
               src/ExactLCPk.cpp
               src/compute_klcp.cpp
               src/compute_kacs.cpp src/main.cpp)

#Link SDSL library
find_package (Threads)
target_link_libraries(${PROJECT_NAME_STR}.x divsufsort divsufsort64 ${CMAKE_THREAD_LIBS_INIT})

# Testing
enable_testing()
# add google test subdir
add_subdirectory("${PROJECT_SOURCE_DIR}/ext/googletest")
set(GTEST_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/ext/googletest/include")
set(GTEST_LIBS_DIR "${CMAKE_BINARY_DIR}/ext/googletest")
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)

add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES}
               src/sdsl/bits.cpp src/sdsl/construct_sa.cpp
               src/ReadsDB.cpp src/AppConfig.cpp
               src/ExactLCPk.cpp
               src/compute_klcp.cpp)
add_dependencies(${PROJECT_TEST_NAME} gtest gtest_main)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/src")
target_link_libraries(${PROJECT_TEST_NAME} ${CMAKE_THREAD_LIBS_INIT}
                      gtest gtest_main divsufsort divsufsort64)
add_test(alfred_test ${PROJECT_TEST_NAME})

# test files
configure_file(${PROJECT_SOURCE_DIR}/test/banana.fa . COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/t1.fa . COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/t2.fa . COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/primates/C.aethiops.fa . COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/primates/C.albifrons.fa . COPYONLY)
# configure_file(${PROJECT_SOURCE_DIR}/test/t2.fq . COPYONLY)
