cmake_minimum_required(VERSION 3.24)
project(particles)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # .dll, .so (shlib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # .a (static)

# Add/enable GoogleTest library
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG b796f7d44681514f58a683a3a71ff17c94edb0c1 # v1.13.0
)
FetchContent_MakeAvailable(googletest)
enable_testing()

add_library(part
        compute_boundary.cc
        droplets_factory.cc
        droplet.cc
        system.cc
        system_evolution.cc
        csv_reader.cc
        csv_writer.cc
        particle.cc
        particles_factory_interface.cc)

add_executable(particles main.cc)
target_link_libraries(particles part)

# Testing
# need to link the test_XXX.cc to "libpart" AND to GTest
# Examples:
#add_executable(test_XXX test_XXX.cc)
#target_link_libraries(test_XXX GTest::gtest_main part)

# Enable tests discovery
#include(GoogleTest)
#gtest_discover_tests(test_XXX) # the target must be an executable known to CMake