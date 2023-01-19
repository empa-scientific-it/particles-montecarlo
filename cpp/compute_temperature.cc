#include "compute_temperature.hh"
#include "material_point.hh"
#include <cmath>

/* -------------------------------------------------------------------------- */
ComputeTemperature::ComputeTemperature(Real timestep, int side) : dt(timestep) {
    // once we have the size of the problem, we can pre-compute the frequencies
    Matrix<complex> qq = FFT::computeFrequencies(side);
    this->q.resize(side);
    for (UInt i = 0; i < side; ++i) {
        for (UInt j = 0; j < side; ++j) {
            this->q(i, j) = pow(qq(i, j).real(), 2) + pow(qq(i, j).imag(), 2);
        }
    }
    // normalize frequencies
    this->q /= (side*side);
}

void ComputeTemperature::setDt(Real dt) {
    // set time step
    this->dt = dt;
}

/*!
 \param kappa : heat conductivity (Real)
 \param heat_c : specific heat capacity (Real)

 method to initialize the parameters for the heat equation solver
 */
void ComputeTemperature::setParams(Real kappa, Real heat_c) {
    // set solver parameters
    this->kappa = kappa;
    this->heat_c = heat_c;
}

void ComputeTemperature::compute(System& system) {
    // Build matrix from our list of particles
    UInt N = sqrt(system.getNbParticles());

    // Matrices for FFT
    Matrix<complex> m(N);
    Matrix<complex> m_fft(N);
    Matrix<complex> hv(N);
    Matrix<complex> hv_fft(N);

    // Loop over particles and extract temperature and 
    // heat flux into matrices
    UInt np = 0;
    for (UInt i = 0; i < N; ++i) {
        for (UInt j = 0; j < N; ++j) {
            auto& p = static_cast<MaterialPoint&>( system.getParticle(np) );
            m(i, j) = p.getTemperature();
            hv(i, j) = p.getHeatRate();
            np += 1;
        }
    }


    // Perform FFT
    m_fft = FFT::transform(m);
    hv_fft = FFT::transform(hv);

    // Loop over particles and compute derivative
    np = 0;
    for (UInt i = 0; i < N; ++i) {
        for (UInt j = 0; j < N; ++j) {
            auto& p = static_cast<MaterialPoint&>( system.getParticle(np) );

            // Here we consider the input particle mass as the mass density, for simplicity
            // (also allows for spatially varying mass density)
            m_fft(i, j) = 1.0/(p.getMass()*heat_c)*(hv_fft(i, j)-kappa*m_fft(i, j)*q(i, j));
            np += 1;
        }
    }

    // Perform inverse FFT
    m = FFT::itransform(m_fft);

    // Loop over particles and perform the Euler integration
    // Don't update the temperature of "boundary" points
    np = 0;
    for (UInt i = 0; i < N; ++i) {
        for (UInt j = 0; j < N; ++j) {
            auto& p = static_cast<MaterialPoint&>( system.getParticle(np) );
            if ( p.isBoundary() ) p.setTemperature(0.0);
            else p.setTemperature( p.getTemperature() + dt * m(i,j).real() );
            np += 1;
        }
    }

}

/* -------------------------------------------------------------------------- */
