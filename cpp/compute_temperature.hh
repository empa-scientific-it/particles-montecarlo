#ifndef __COMPUTE_TEMPERATURE__HH__
#define __COMPUTE_TEMPERATURE__HH__

/* -------------------------------------------------------------------------- */
#include "compute.hh"
#include "fft.hh"

//! Compute contact interaction between ping-pong balls
class ComputeTemperature : public Compute {
  // Virtual implementation
  public:
    ComputeTemperature(Real timestep, int side);
  public:
  // override of compute method
    void compute(System& system) override;
    void setDt(Real dt);
    void setParams(Real kappa, Real heat_c);
// Private members
  private:
    Matrix<Real> q; // FFT frequencies
    Real dt;
    Real kappa;
    Real heat_c;
};

/* -------------------------------------------------------------------------- */
#endif  //__COMPUTE_TEMPERATURE__HH__
