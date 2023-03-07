#ifndef PARTICLES_COMPUTE_EVAPORATION_HH
#define PARTICLES_COMPUTE_EVAPORATION_HH

#include "compute.hh"

class ComputeEvaporation : public Compute {
public:
    // Constructor
    ComputeEvaporation(Real temperature, Real relative_h);
public:
    // Methods
    void compute(System& system) override;
    void setParams(Real temperature, Real relative_h);
private:
    // User inputs
    Real temperature, relative_h;
    // Calculated params
    Real wet_bulb_temp, vapor_diff, vapor_press, film_temp;
    Real alpha, beta;
    // Minimum attainable diameter (hard-coded)
    const Real d_min = 1.0e-7;
};


#endif //PARTICLES_COMPUTE_EVAPORATION_HH
