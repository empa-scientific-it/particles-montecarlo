#include "droplets_factory.hh"
#include <iostream>
#include "csv_reader.hh"

std::unique_ptr<Particle> DropletsFactory::createParticle() {
    return std::make_unique<Droplet>();
}

SystemEvolution& DropletsFactory::createSimulation(const std::string &fname, Real timestep) {

    this->system_evolution = std::make_unique<SystemEvolution>(std::make_unique<System>());

    CsvReader reader(fname);
    reader.read(this->system_evolution->getSystem());

    return *system_evolution;
}

ParticlesFactoryInterface& DropletsFactory::getInstance() {
    if (not ParticlesFactoryInterface::factory)
        ParticlesFactoryInterface::factory = new DropletsFactory;

    return *factory;
}