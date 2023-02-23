cmake_minimum_required(VERSION 3.24)
project(particles)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(particles
        main.cc
        compute_boundary.cc
        droplets_factory.cc
        droplet.cc
        system.cc
        system_evolution.cc
        compute_interaction.cc
        csv_reader.cc
        csv_writer.cc
        particle.cc
        particles_factory_interface.cc)
