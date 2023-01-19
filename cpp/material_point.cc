#include "material_point.hh"

/* -------------------------------------------------------------------------- */
void MaterialPoint::printself(std::ostream& stream) const {
  Particle::printself(stream);
  stream << " " << temperature << " " << heat_rate << " " << is_boundary;
}

/* -------------------------------------------------------------------------- */

void MaterialPoint::initself(std::istream& sstr) {
  Particle::initself(sstr);
  sstr >> temperature  >> heat_rate >> is_boundary;
}

/* -------------------------------------------------------------------------- */

void MaterialPoint::setTemperature(const Real& new_temperature) {
    temperature = new_temperature;
}

/* -------------------------------------------------------------------------- */

void MaterialPoint::setHeatRate(const Real& new_heat_rate) {
    heat_rate = new_heat_rate;
}
