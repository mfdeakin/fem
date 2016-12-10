
#include <iomanip>
#include <iostream>

#include <random>

#include <array.hpp>
#include <ctmath.hpp>
#include <polynomial.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace Numerical;

TEST_CASE("Polynomial Evaluation", "[Polynomial]") {
  constexpr const int degree_range = 2;
  constexpr const int dim = 3;
  using CoeffT = double;
  using P = Polynomial<CoeffT, degree_range, dim>;
  P p;
  p.coeff(0, 0, 0) = 2.0;

  p.coeff(0, 0, 1) = 3.0;
  p.coeff(0, 1, 0) = 5.0;
  p.coeff(1, 0, 0) = 7.0;

  p.coeff(0, 0, 2) = 0.0;
  p.coeff(0, 1, 1) = 0.0;
  p.coeff(0, 2, 0) = 0.0;
  p.coeff(1, 0, 1) = 0.0;
  p.coeff(1, 1, 0) = 0.0;
  p.coeff(2, 0, 0) = 0.0;

  REQUIRE(p.eval(0.0, 0.0, 0.0) == 2.0);

  REQUIRE(p.eval(0.0, 0.0, 1.0) == 5.0);
  REQUIRE(p.eval(0.0, 1.0, 0.0) == 7.0);
  REQUIRE(p.eval(1.0, 0.0, 0.0) == 9.0);
}

TEST_CASE("Coefficient Indices", "[Polynomial]") {
  std::random_device rd;
  std::mt19937_64 engine(rd());
  using CoeffT = double;
  using pdf_uniform =
      std::uniform_real_distribution<CoeffT>;

  constexpr const int num_rand_tests = 10;

  SECTION("degree 2, dimension 2") {
    constexpr const int degree_range = 2;
    constexpr const int dim = 2;
    using P = Polynomial<CoeffT, degree_range, dim>;
    using A = Array<int, dim>;
    Array<CoeffT,
          CTMath::poly_num_coeffs<int>(degree_range, dim)>
        coeffs;
    P p;
    for(int i = 0; i < num_rand_tests; ++i) {
      coeffs[0] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 0) = coeffs[0];

      coeffs[1] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 1) = coeffs[1];
      coeffs[2] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 0) = coeffs[2];

      coeffs[3] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 2) = coeffs[3];
      coeffs[4] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 1) = coeffs[4];
      coeffs[5] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(2, 0) = coeffs[5];

      REQUIRE(p.coeff(0, 0) == coeffs[0]);
      REQUIRE(p.coeff(A(0, 0)) == coeffs[0]);

      REQUIRE(p.coeff(0, 1) == coeffs[1]);
      REQUIRE(p.coeff(A(0, 1)) == coeffs[1]);
      REQUIRE(p.coeff(1, 0) == coeffs[2]);
      REQUIRE(p.coeff(A(1, 0)) == coeffs[2]);

      REQUIRE(p.coeff(0, 2) == coeffs[3]);
      REQUIRE(p.coeff(A(0, 2)) == coeffs[3]);
      REQUIRE(p.coeff(1, 1) == coeffs[4]);
      REQUIRE(p.coeff(A(1, 1)) == coeffs[4]);
      REQUIRE(p.coeff(2, 0) == coeffs[5]);
      REQUIRE(p.coeff(A(2, 0)) == coeffs[5]);
    }
  }

  SECTION("degree 3, dimension 3") {
    constexpr const int degree_range = 3;
    constexpr const int dim = 3;
    using P = Polynomial<CoeffT, degree_range, dim>;
    using A = Array<int, dim>;
    Array<CoeffT,
          CTMath::poly_num_coeffs<int>(degree_range, dim)>
        coeffs;
    P p;
    for(int i = 0; i < num_rand_tests; ++i) {
      coeffs[0] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 0, 0) = coeffs[0];

      coeffs[1] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 0, 1) = coeffs[1];
      coeffs[2] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 1, 0) = coeffs[2];
      coeffs[3] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 0, 0) = coeffs[3];

      coeffs[4] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 0, 2) = coeffs[4];
      coeffs[5] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 1, 1) = coeffs[5];
      coeffs[6] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 2, 0) = coeffs[6];
      coeffs[7] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 0, 1) = coeffs[7];
      coeffs[8] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 1, 0) = coeffs[8];
      coeffs[9] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(2, 0, 0) = coeffs[9];

      coeffs[10] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 0, 3) = coeffs[10];
      coeffs[11] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 1, 2) = coeffs[11];
      coeffs[12] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 2, 1) = coeffs[12];
      coeffs[13] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(0, 3, 0) = coeffs[13];
      coeffs[14] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 0, 2) = coeffs[14];
      coeffs[15] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 1, 1) = coeffs[15];
      coeffs[16] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(1, 2, 0) = coeffs[16];
      coeffs[17] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(2, 0, 1) = coeffs[17];
      coeffs[18] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(2, 1, 0) = coeffs[18];
      coeffs[19] = pdf_uniform(-1.0, 1.0)(engine);
      p.coeff(3, 0, 0) = coeffs[19];

      REQUIRE(p.coeff(0, 0, 0) == coeffs[0]);
      REQUIRE(p.coeff(A(0, 0, 0)) == coeffs[0]);

      REQUIRE(p.coeff(0, 0, 1) == coeffs[1]);
      REQUIRE(p.coeff(A(0, 0, 1)) == coeffs[1]);
      REQUIRE(p.coeff(0, 1, 0) == coeffs[2]);
      REQUIRE(p.coeff(A(0, 1, 0)) == coeffs[2]);
      REQUIRE(p.coeff(1, 0, 0) == coeffs[3]);
      REQUIRE(p.coeff(A(1, 0, 0)) == coeffs[3]);

      REQUIRE(p.coeff(0, 0, 2) == coeffs[4]);
      REQUIRE(p.coeff(A(0, 0, 2)) == coeffs[4]);
      REQUIRE(p.coeff(0, 1, 1) == coeffs[5]);
      REQUIRE(p.coeff(A(0, 1, 1)) == coeffs[5]);
      REQUIRE(p.coeff(0, 2, 0) == coeffs[6]);
      REQUIRE(p.coeff(A(0, 2, 0)) == coeffs[6]);
      REQUIRE(p.coeff(1, 0, 1) == coeffs[7]);
      REQUIRE(p.coeff(A(1, 0, 1)) == coeffs[7]);
      REQUIRE(p.coeff(1, 1, 0) == coeffs[8]);
      REQUIRE(p.coeff(A(1, 1, 0)) == coeffs[8]);
      REQUIRE(p.coeff(2, 0, 0) == coeffs[9]);
      REQUIRE(p.coeff(A(2, 0, 0)) == coeffs[9]);

      REQUIRE(p.coeff(0, 0, 3) == coeffs[10]);
      REQUIRE(p.coeff(A(0, 0, 3)) == coeffs[10]);
      REQUIRE(p.coeff(0, 1, 2) == coeffs[11]);
      REQUIRE(p.coeff(A(0, 1, 2)) == coeffs[11]);
      REQUIRE(p.coeff(0, 2, 1) == coeffs[12]);
      REQUIRE(p.coeff(A(0, 2, 1)) == coeffs[12]);
      REQUIRE(p.coeff(0, 3, 0) == coeffs[13]);
      REQUIRE(p.coeff(A(0, 3, 0)) == coeffs[13]);
      REQUIRE(p.coeff(1, 0, 2) == coeffs[14]);
      REQUIRE(p.coeff(A(1, 0, 2)) == coeffs[14]);
      REQUIRE(p.coeff(1, 1, 1) == coeffs[15]);
      REQUIRE(p.coeff(A(1, 1, 1)) == coeffs[15]);
      REQUIRE(p.coeff(1, 2, 0) == coeffs[16]);
      REQUIRE(p.coeff(A(1, 2, 0)) == coeffs[16]);
      REQUIRE(p.coeff(2, 0, 1) == coeffs[17]);
      REQUIRE(p.coeff(A(2, 0, 1)) == coeffs[17]);
      REQUIRE(p.coeff(2, 1, 0) == coeffs[18]);
      REQUIRE(p.coeff(A(2, 1, 0)) == coeffs[18]);
      REQUIRE(p.coeff(3, 0, 0) == coeffs[19]);
      REQUIRE(p.coeff(A(3, 0, 0)) == coeffs[19]);
    }
  }
}
