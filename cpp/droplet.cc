#include "droplet.hh"

void Droplet::printself(std::ostream &stream) const {
    Particle::printself(stream);
    stream << " " << diameter << " " << is_active;
}

/*------------------------------------------------------------*/

void Droplet::initself(std::istream &sstr) {
    Particle::initself(sstr);
    // droplet default state is active (i.e. not settled)
    is_active = true;
    sstr >> diameter;
}

/*------------------------------------------------------------*/

void Droplet::setState(const bool &new_state) {
    is_active = new_state;
}