#ifndef FFT_HH
#define FFT_HH
/* ------------------------------------------------------ */
#include "matrix.hh"
#include "my_types.hh"
#include <fftw3.h>
/* ------------------------------------------------------ */

struct FFT {

  static Matrix<complex> transform(Matrix<complex>& m);
  static Matrix<complex> itransform(Matrix<complex>& m);

  static Matrix<complex> computeFrequencies(int size);

};

/* ------------------------------------------------------ */

inline Matrix<complex> FFT::transform(Matrix<complex>& m_in) {
    int n = m_in.size();
    Matrix<complex> m_out(n);
    fftw_plan p;

    // Do the 2D DFT with FFTW, casting the complex into fftw_complex
    p = fftw_plan_dft_2d(n, n,
                         reinterpret_cast<fftw_complex*>(m_in.data()),
                         reinterpret_cast<fftw_complex*>(m_out.data()),
                         FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    fftw_destroy_plan(p);

    return m_out;
}

/* ------------------------------------------------------ */

inline Matrix<complex> FFT::itransform(Matrix<complex>& m_in) {
    int n = m_in.size();
    Matrix<complex> m_out(n);
    fftw_plan p;

    // Do the 2D DFT with FFTW, casting the complex into fftw_complex
    p = fftw_plan_dft_2d(n, n,
                        reinterpret_cast<fftw_complex*>(m_in.data()),
                        reinterpret_cast<fftw_complex*>(m_out.data()),
                        FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    fftw_destroy_plan(p);

    // normalization factor for the backward FFT
    m_out /= (n*n);

    return m_out;
}

/* ------------------------------------------------------ */


/* ------------------------------------------------------ */

inline Matrix<complex> FFT::computeFrequencies(int size) {
    // Computes non-normalized FFT frequencies
    // in the same manner as np.fft.fftfreq
    //
    // The x frequencies are stored as the real part of
    // complex number; the y frequencies are stored
    // as the imaginary part
    
    Matrix<complex> m_out(size);
    int ii;
    int jj;
    
    // Loop over the four quadrants of the matrix
    // to assign frequencies.
    
    // for even side length
    if (size % 2 == 0) {
        for (int i = 0; i < size/2; ++i) {
            for (int j = 0; j < size/2; ++j) {
                m_out(i, j).real(-i);
                m_out(i, j).imag(-j);
            }
        }
        ii = 0;
        for (int i = size/2; i < size; ++i) {
            jj = 0;
            for (int j = size/2; j < size; ++j) {
                m_out(i, j).real(-i+ii);
                m_out(i, j).imag(-j+jj);
                jj += 2;
            }
            ii += 2;
        }
        for (int i = 0; i < size/2; ++i) {
            jj = 0;
            for (int j = size/2; j < size; ++j) {
                m_out(i, j).real(i);
                m_out(i, j).imag(-j+jj);
                jj += 2;
            }
        }
        ii = 0;
        for (int i = size/2; i < size; ++i) {
            for (int j = 0; j < size/2; ++j) {
                m_out(i, j).real(-i+ii);
                m_out(i, j).imag(j);
            }
            ii += 2;
        }

    // for odd side length
    } else {
        for (int i = 0; i <= (size-1)/2; ++i) {
            for (int j = 0; j <= (size-1)/2; ++j) {
                m_out(i, j).real(i);
                m_out(i, j).imag(j);
            }
        }
        ii = 1;
        for (int i = (size+1)/2; i < size; ++i) {
            jj = 1;
            for (int j = (size+1)/2; j < size; ++j) {
                m_out(i, j).real(-i+ii);
                m_out(i, j).imag(-j+jj);
                jj += 2;
            }
            ii += 2;
        }
        for (int i = 0; i <= (size-1)/2; ++i) {
            jj = 1;
            for (int j = (size+1)/2; j < size; ++j) {
                m_out(i, j).real(i);
                m_out(i, j).imag(-j+jj);
                jj += 2;
            }
        }
        ii = 1;
        for (int i = (size+1)/2; i < size; ++i) {
            for (int j = 0; j <= (size-1)/2; ++j) {
                m_out(i, j).real(-i+ii);
                m_out(i, j).imag(j);
            }
            ii += 2;
        }
    }

    // return non-normalized frequencies
    return m_out;
}

#endif  // FFT_HH
