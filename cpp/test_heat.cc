#include "my_types.hh"
#include "material_point.hh"
#include "material_points_factory.hh"
#include "system.hh"
#include "csv_reader.hh"
#include "csv_writer.hh"
#include "fft.hh"
#include "compute_temperature.hh"
#include <gtest/gtest.h>

/*****************************************************************************/
class InitialConditions : public ::testing::Test {
    protected:
        void SetUp() override {

            // System size and timestep
            N = 144;
            dt = 0.1;
            
            // Initialize the points 
            std::vector<MaterialPoint> points;
            for (UInt i = 0; i < N; ++i) {
                MaterialPoint p;
                p.setTemperature(1.0);
                p.getMass() = 1.0;
                p.isBoundary() = 0;
                points.push_back(p);
            }

            // Add particles to system
            for (auto& p : points) {
                system.addParticle(std::make_shared<MaterialPoint>(p));
            }

            // Set the compute
            temperature = std::make_shared<ComputeTemperature>(dt, N);
            temperature->setParams(1.0, 1.0);
        }

        System system;
        UInt N;
        double dt;

        std::shared_ptr<ComputeTemperature> temperature;
};
/*****************************************************************************/
TEST_F(InitialConditions, homogeneous) {

    // No heat flux, uniform temperature
    for (UInt i = 0; i < N; ++i) {
        MaterialPoint &p = dynamic_cast<MaterialPoint&>(system.getParticle(i));
        p.setHeatRate(0.0);
    }

    temperature->compute(system);

    // Temperature should not change after a time step
    for (UInt i = 0; i < N; ++i) {
        MaterialPoint &p = dynamic_cast<MaterialPoint&>(system.getParticle(i));
        ASSERT_NEAR(p.getTemperature(), 1.0, 1.0e-15);
    }
}

TEST_F(InitialConditions, step) {

    // Delta function heat flux, uniform temperature
    UInt np = 0;
    for (UInt i = 0; i < sqrt(N); ++i) {
        for (UInt j = 0; j < sqrt(N); ++j) {
            MaterialPoint &p = dynamic_cast<MaterialPoint&>(system.getParticle(np));

            // Set the heat rates
            if (j == sqrt(N)/4) {
                p.setHeatRate(-sqrt(N));
            } else if (j == 3*sqrt(N)/4) {
                p.setHeatRate(sqrt(N));
            } else {
                p.setHeatRate(0.0);
            }
            np += 1;
        }
    }

    temperature->compute(system);

    // Check temperatures after one time step
    np = 0;
    for (UInt i = 0; i < sqrt(N); ++i) {
        for (UInt j = 0; j < sqrt(N); ++j) {
            MaterialPoint &p = dynamic_cast<MaterialPoint&>(system.getParticle(np));
            if (j == sqrt(N)/4) {
                ASSERT_NEAR(p.getTemperature(), -0.2, 1e-10);
            } else if (j == 3*sqrt(N)/4) {
                ASSERT_NEAR(p.getTemperature(), 2.2, 1e-10);
            } else {
                ASSERT_NEAR(p.getTemperature(), 1.0, 1e-10);
            }
            np += 1;
        }
    }
}

TEST_F(InitialConditions, radial) {

    // Set parameters for radial heat distribution
    UInt np;
    Real xmin = -1.0;
    Real ymin = -1.0;
    Real xmax = 1.0;
    Real ymax = 1.0;
    Real r = 0.5;
    Real dx = (xmax-xmin)/(sqrt(N)-1.0);
    Real dy = (ymax-ymin)/(sqrt(N)-1.0);
    Real x;
    Real y;
    
    // Set the heat rates
    np = 0;
    y = ymin;
    for (UInt i = 0; i < sqrt(N); ++i) {
        x = xmin;
        for (UInt j = 0; j < sqrt(N); ++j) {
            MaterialPoint &p = dynamic_cast<MaterialPoint&>(system.getParticle(np));
            np += 1;
            if ((x*x + y*y) < r) {
                p.setHeatRate(1.0);
            } else {
                p.setHeatRate(0.0);
            }
            x += dx;
        }
        y += dy;
    }

    temperature->compute(system);

    // Check temperatures after one time step
    np = 0;
    y = ymin;
    for (UInt i = 0; i < sqrt(N); ++i) {
        x = xmin;
        for (UInt j = 0; j < sqrt(N); ++j) {
            MaterialPoint &p = dynamic_cast<MaterialPoint&>(system.getParticle(np));
            np += 1;
            if ((x*x + y*y) < r) {
                ASSERT_NEAR(p.getTemperature(), 1.1, 1e-10);
            } else {
                ASSERT_NEAR(p.getTemperature(), 1.0, 1e-10);
            }
            x += dx;
        }
        y += dy;
    }

}
