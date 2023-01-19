#include "my_types.hh"
#include "fft.hh"
#include <gtest/gtest.h>

/*****************************************************************/
TEST(FFT, transform) {
  UInt N = 512;
  Matrix<complex> m(N);

  Real k = 2 * M_PI / N;
  for (auto&& entry : index(m)) {
    int i = std::get<0>(entry);
    int j = std::get<1>(entry);
    auto& val = std::get<2>(entry);
    val = cos(k * i);
  }

  Matrix<complex> res = FFT::transform(m);

  for (auto&& entry : index(res)) {
    int i = std::get<0>(entry);
    int j = std::get<1>(entry);
    auto& val = std::get<2>(entry);
    if (std::abs(val) > 1e-10)
      std::cout << i << "," << j << " = " << val << std::endl;

    if (i == 1 && j == 0)
      ASSERT_NEAR(std::abs(val), N * N / 2, 1e-10);
    else if (i == N - 1 && j == 0)
      ASSERT_NEAR(std::abs(val), N * N / 2, 1e-10);
    else
      ASSERT_NEAR(std::abs(val), 0, 1e-10);
  }
}
/*****************************************************************/

TEST(FFT, inverse_transform) {
  UInt N = 512;
  Matrix<complex> m(N);

  Real k = 2 * M_PI / N;
  for (auto&& entry : index(m)) {
    int i = std::get<0>(entry);
    int j = std::get<1>(entry);
    auto& val = std::get<2>(entry);

    if (i == 1 && j == 0)
        val = N * N / 2;
    else if (i == N - 1 && j == 0)
        val = N * N / 2;
    else
        val = 0.0;
  }

  Matrix<complex> res = FFT::itransform(m);

  //res /= N;

  for (auto&& entry : index(res)) {
    int i = std::get<0>(entry);
    int j = std::get<1>(entry);
    auto& val = std::get<2>(entry);

    if (val.imag() > 1e-10)
        std::cout << val << std::endl;
    
    ASSERT_NEAR(val.real(), cos(k * i), 1e-10);
  }
}
/*****************************************************************/


TEST(FFT, frequencies) {
  UInt N = 3;
  UInt ii;
  UInt jj;

  Matrix<complex> res = FFT::computeFrequencies(N);
  Matrix<Real> q(N);

  for (UInt i = 0; i < N; ++i) {
      for (UInt j = 0; j < N; ++j) {
          q(i, j) = pow(res(i, j).real(), 2) + pow(res(i, j).imag(), 2);
      }
  }

  ASSERT_NEAR(q(0, 0), 0.0, 1e-10);
  ASSERT_NEAR(q(0, 1), 1.0, 1e-10);
  ASSERT_NEAR(q(0, 2), 1.0, 1e-10);
  ASSERT_NEAR(q(1, 0), 1.0, 1e-10);
  ASSERT_NEAR(q(1, 1), 2.0, 1e-10);
  ASSERT_NEAR(q(1, 2), 2.0, 1e-10);
  ASSERT_NEAR(q(2, 0), 1.0, 1e-10);
  ASSERT_NEAR(q(2, 1), 2.0, 1e-10);
  ASSERT_NEAR(q(2, 2), 2.0, 1e-10);

}

