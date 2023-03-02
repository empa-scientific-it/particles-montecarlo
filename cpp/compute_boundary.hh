#ifndef __COMPUTE_BOUNDARY__HH__
#define __COMPUTE_BOUNDARY__HH__

/* -------------------------------------------------------------------------- */
#include "compute.hh"
/* -------------------------------------------------------------------------- */

//! Compute interaction with simulation box
class ComputeBoundary : public Compute {
  // Constructors/Destructors
public:
  ComputeBoundary(const Vector& box_min, const Vector& box_max);

  // Methods
public:
  void compute(System& system) override;

  // Members
protected:
  Vector box_min, box_max;
};

/* -------------------------------------------------------------------------- */
#endif  //__COMPUTE_BOUNDARY__HH__
