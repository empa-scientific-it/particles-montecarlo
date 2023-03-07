#include "compute_evaporation.hh"
#include "droplet.hh"
#include <cmath>
#include <constants.hh>

ComputeEvaporation::ComputeEvaporation(Real temperature, Real relative_h) : temperature(temperature), relative_h(relative_h) {

    wet_bulb_temp = temperature*atan(0.151977*pow((relative_h+8.313659),0.5)) + atan(temperature + relative_h) - atan(relative_h - 1.676331) + 0.00391838*pow(relative_h,1.5)*atan(0.023101*relative_h) - 4.686035;

    vapor_diff = 21.2e-6 * (1 + 0.0071 * wet_bulb_temp);
    vapor_press = 67 * (temperature - wet_bulb_temp);
    film_temp = wetBulb + 1/3*(temperature-wet_bulb_temp);

    alpha = 4 * WATER_MOL_W * vapor_diff * vapor_press / (DENSITY_L * R * (film_temp + 273.15));
    beta = 0.276 * pow ((DENSITY_G / (VISC_G * pow(vapor_diff,2))),(1 / 6));

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
        auto &diam = par.getDiameter();
        auto &is_active = par.getState();

        if (!is_active)
            break;

        Real sqrt_arg = pow(diam,2) - 2 * dt * alpha_evap * (1 + beta_evap * sqrt(diam * vel));
        Real temp_diam = sqrt_arg < 0 ? D_MIN : sqrt(sqrt_arg);

        //TODO: add calculation of rate of change of mass, Mp, needed for blowing velocity

        diam = temp_diam;
    }
}

