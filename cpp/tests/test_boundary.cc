#include "../droplet.hh"
#include "../system.hh"
#include "../droplets_factory.hh"
#include "../compute_boundary.hh"
#include <gtest/gtest.h>
/*******************************************************************/

// Fixture class: two droplets arbitrarily close to the boundary
class TwoDroplets : public ::testing::Test {
protected:
    void SetUp() override {
        DropletsFactory::getInstance();
        Droplet dp1, dp2;
        // Set mass (fake)
        dp1.getMass() = 1.;
        dp2.getMass() = 1.;
        // Set diameter (fake)
        dp1.getDiameter() = 1.;
        dp2.getDiameter() = 1.;
        // Set zero velocities
        dp1.getVelocity() = 0.;
        dp2.getVelocity() = 0.;
        // Set positions
        dp1.getPosition() = 0.;
        dp2.getPosition() = 0.;

        // Add the particles to the system
        system.addParticle(std::make_shared<Droplet>(dp1));
        system.addParticle(std::make_shared<Droplet>(dp2));

        // Set the boundary
        box_min = 0.0;
        box_max = 3.0;

        // Boundary interaction
        boundary = std::make_shared<ComputeBoundary>(box_min, box_max);
    }

    System system;
    Vector box_min, box_max;
    std::shared_ptr<ComputeBoundary> boundary;
};
/*******************************************************************/

TEST_F(TwoDroplets, BoundaryExit) {
    auto dp1 = system.getParticle(0);
    auto dp2 = system.getParticle(1);

    // TODO: test whether a droplet bounces back correctly
    //boundary->compute(system);

    // Dummy test
    EXPECT_EQ(dp1.getMass(), dp2.getMass());
    EXPECT_EQ(system.getNbParticles(), 3);
}
