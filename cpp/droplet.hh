#ifndef PARTICLES_DROPLET_H
#define PARTICLES_DROPLET_H

#include "particle.hh"

class Droplet : public Particle {

public:
    //! Get droplet radius
    Real& getDiameter() { return diameter; }

    //! Get droplet state
    bool& getState() { return is_active; }

    void printself(std::ostream& stream) const override;
    void initself(std::istream& sstr) override;

private:
    Real diameter;
    bool is_active;
};

#endif //PARTICLES_DROPLET_H
