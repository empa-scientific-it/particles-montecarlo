#ifndef __MATERIAL_POINT__HH__
#define __MATERIAL_POINT__HH__

/* -------------------------------------------------------------------------- */
#include "particle.hh"

//! Class for MaterialPoint
class MaterialPoint : public Particle {
/* ------------------------------------------------------------------------ */
/* Methods                                                                  */
/* ------------------------------------------------------------------------ */

public:

  void printself(std::ostream& stream) const override;
  void initself(std::istream& sstr) override;

  // temperature
  const Real& getTemperature() const {return temperature;};
  void setTemperature(const Real& new_temperature);

  // heat rate
  const Real& getHeatRate() const {return heat_rate;};
  void setHeatRate(const Real& new_heat_rate);

  // get "boundary" points
  bool& isBoundary() {return is_boundary;};
  
/* ------------------------------------------------------------------------ */
/* Attributes                                                               */
/* ------------------------------------------------------------------------ */
private:
  Real temperature;
  Real heat_rate;
  bool is_boundary;
};

/* -------------------------------------------------------------------------- */
#endif  //__MATERIAL_POINT__HH__
