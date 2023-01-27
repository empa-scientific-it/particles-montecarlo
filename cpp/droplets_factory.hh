#ifndef PARTICLES_DROPLETS_FACTORY_HH
#define PARTICLES_DROPLETS_FACTORY_HH

#include "particles_factory_interface.hh"
#include "droplet.hh"

class DropletsFactory : public ParticlesFactoryInterface {

private:
    DropletsFactory() = default;

public:
    SystemEvolution& createSimulation(const std::string& fname, Real timestep) override;

    std::unique_ptr<Particle> createParticle() override;

    static ParticlesFactoryInterface& getInstance();
};


#endif //PARTICLES_DROPLETS_FACTORY_HH
