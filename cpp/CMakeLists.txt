cmake_minimum_required(VERSION 3.24)
project(particles)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # .dll, .so (shlib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # .a (static)

# Build a (static) library
add_library(part
        compute_boundary.cc
        droplets_factory.cc
        droplet.cc
        system.cc
        system_evolution.cc
        csv_reader.cc
        csv_writer.cc
        particle.cc
        particles_factory_interface.cc
        compute_verlet_integration.cc
        compute_evaporation.cc)

# Build/link main executable
add_executable(particles main.cc)
target_link_libraries(particles part)

# Tests
enable_testing()
add_subdirectory(tests)