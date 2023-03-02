#include "../droplet.hh"
#include "../system.hh"
#include "../droplets_factory.hh"
#include "../compute_boundary.hh"
#include <gtest/gtest.h>
#include <stdexcept>
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
        // Set state
        dp1.setState(true);
        dp2.setState(true);
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

// Inactive particle
TEST_F(TwoDroplets, InactiveParticle) {
    auto &dp1 = dynamic_cast<Droplet &>(system.getParticle(0));
    dp1.getPosition()[2] = -1.0;
    dp1.getPosition()[0] = -1.0;

    boundary->compute(system);

    EXPECT_FALSE(dp1.getState());
    EXPECT_NEAR(dp1.getPosition()[0], -1.0, 1e-10);
}

// Boundary exit
TEST_F(TwoDroplets, BoundaryExitX) {
    auto &dp1 = system.getParticle(0);
    auto &dp2 = system.getParticle(1);

    dp1.getPosition()[0] = 3.2;
    dp1.getVelocity()[0] = 0.3;
    dp2.getPosition()[0] = -0.1;
    dp2.getVelocity()[0] = -0.2;

    boundary->compute(system);

    EXPECT_LT(dp1.getVelocity()[0], 0);
    EXPECT_GT(dp2.getVelocity()[0], 0);

    EXPECT_NEAR(dp1.getPosition()[0], 2.8, 1e-5);
    EXPECT_NEAR(dp2.getPosition()[0], 0.1, 1e-5);
}

// Boundary exit along XY
TEST_F(TwoDroplets, BoundaryExitXY) {
    auto &dp1 = system.getParticle(0);
    dp1.getPosition()[0] = 3.6;
    dp1.getPosition()[1] = -0.2;
    dp1.getVelocity()[0] = 0.5;
    dp1.getVelocity()[1] = -0.5;

    boundary->compute(system);

    EXPECT_NEAR(dp1.getPosition()[0], 2.4, 1e-5);
    EXPECT_NEAR(dp1.getPosition()[1], 0.2, 1e-5);
}

// Error: box outside boundary
TEST_F(TwoDroplets, OutOfBound) {
    auto &dp1 = system.getParticle(0);
    dp1.getPosition()[0] = 6.2;
    EXPECT_THROW(boundary->compute(system), std::runtime_error);
}
