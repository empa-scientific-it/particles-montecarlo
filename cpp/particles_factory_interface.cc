#include "particles_factory_interface.hh"
/* -------------------------------------------------------------------------- */
ParticlesFactoryInterface& ParticlesFactoryInterface::getInstance() {

  return *factory;

}

/* -------------------------------------------------------------------------- */
ParticlesFactoryInterface* ParticlesFactoryInterface::factory = nullptr;
