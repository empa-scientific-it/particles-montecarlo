#include <cmath>
#include <stdexcept>
#include "compute_boundary.hh"
#include "droplet.hh"

/* -------------------------------------------------------------------------- */

ComputeBoundary::ComputeBoundary(const Vector &box_min, const Vector &box_max)
    : box_min(box_min), box_max(box_max)
{
    Vector d = box_max - box_min;
    for (UInt i = 0; i < Vector::dim; ++i)
        if (d[i] < 0)
        {
            std::cout << "Box extents do not form a domain range. Min values exceed Max values." << std::endl;
            std::exit(1);
        }
}

/* -------------------------------------------------------------------------- */

void ComputeBoundary::compute(System &system) {
    UInt num_particles = system.getNbParticles();
    for (UInt np = 0; np < num_particles; ++np) {
        auto &par = dynamic_cast<Droplet &>(system.getParticle(np));
        auto &pos = par.getPosition();
        auto &vel = par.getVelocity();
        auto &is_active = par.getState();

        Vector delta_max = box_max - pos, delta_min = pos - box_min;

        // check if the particle has gone past the min Z boundary
        if (delta_min[2] < 0)
            par.setState(false);

        // if particle is inactive, skip
        if (!is_active)
            break;

        // updated particle's position and velocity
        for (UInt i = 0; i < Vector::dim; ++i) {
            if (delta_max[i] < 0) {
                vel[i] *= -1;
                pos[i] -= 2 * abs(delta_max[i]);
            } else if (delta_min[i] < 0) {
                vel[i] *= -1;
                pos[i] += 2 * abs(delta_min[i]);
            }
            // check that particle's new position is inside the box
            if (box_max[i] - pos[i] < 0 || pos[i] - box_min[i] < 0)
                throw std::runtime_error("Particle position outside the bounding box!");
        }
    }
}
/* -------------------------------------------------------------------------- */