#include "compute_evaporation.hh"
#include "droplet.hh"
#include <cmath>

ComputeEvaporation::ComputeEvaporation(Real temperature, Real relative_h) : temperature(temperature), relative_h(relative_h) {
    // TODO: Compute the wet-bulb temperature
    // d0 = 0.151977
    // d1 = 8.313659
    // d2 = 1.676331
    // d3 = 0.00391838
    // d4 = 0.023101
    // d5 = 4.686035
    wet_bulb_temp = temperature * atan(0.151977 * 0.0);
    vapor_diff = 0.0;
    vapor_press = 0.0;
    film_temp = 0.0;

    alpha = 0.0;
    beta = 0.0;
}

void ComputeEvaporation::setParams(Real temperature, Real relative_h) {
    this->temperature = temperature;
    this->relative_h = relative_h;
}

void ComputeEvaporation::compute(System &system) {
    UInt num_particles = system.getNbParticles();
    for (UInt np = 0; np < num_particles; ++np) {
        auto &par = dynamic_cast<Droplet &>(system.getParticle(np));
        auto &pos = par.getPosition();
        auto &vel = par.getVelocity();
        auto &is_active = par.getState();

        if (!is_active)
            break;

        // TODO: implement the `evap` function
    }
}

