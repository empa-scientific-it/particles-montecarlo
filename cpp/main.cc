#include "compute_verlet_integration.hh"
#include "csv_reader.hh"
#include "csv_writer.hh"
#include "my_types.hh"
#include "system.hh"
#include "droplets_factory.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <iostream>
#include <sstream>
/* -------------------------------------------------------------------------- */

int main(int argc, const char** argv) {
  if (argc != 6) {
    std::cout << "Usage: " << argv[0]
              << " nsteps dump_freq input.csv particle_type timestep"
              << std::endl;
    std::cout << "\tparticle type can be: droplet" << std::endl;
    std::exit(EXIT_FAILURE);
  }


  // the number of steps to perform
  UInt num_steps;
  std::stringstream(argv[1]) >> num_steps;
  // freq to dump
  int dump_freq;
  std::stringstream(argv[2]) >> dump_freq;
  // init file
  std::string filename = argv[3];
  // particle type
  std::string type = argv[4];
  // timestep
  Real timestep;
  std::stringstream(argv[5]) >> timestep;

  if (type == "droplet")
    DropletsFactory::getInstance();
  else {
    std::cout << "Unknown particle type: " << type << std::endl;
    std::exit(EXIT_FAILURE);
  }

  ParticlesFactoryInterface& factory = ParticlesFactoryInterface::getInstance();

  SystemEvolution& system_evolution = factory.createSimulation(filename, timestep);

  system_evolution.setNSteps(num_steps);
  system_evolution.setDumpFreq(dump_freq);

    try {
        system_evolution.evolve();
        return EXIT_SUCCESS;
    } catch (...) {
        return EXIT_FAILURE;
    }
}
