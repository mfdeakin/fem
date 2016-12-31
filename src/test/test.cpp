
#include <iomanip>
#include <iostream>

#include <random>

#include <array.hpp>
#include <ctmath.hpp>
#include <polynomial.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace Numerical;

TEST_CASE("Change Polynomial Degree", "[Polynomial]") {
  constexpr const int dim = 3;
  constexpr const int degree_1 = 2;
  constexpr const int degree_2 = 3;
  using CoeffT = double;
  using P1 = Polynomial<CoeffT, degree_1, dim>;
  using P2 = Polynomial<CoeffT, degree_2, dim>;
  SECTION("Degree 3 to 2") {
    P2 p_start((Tags::Zero_Tag()));
    p_start.coeff(0, 0, 0) = 1.0;

    p_start.coeff(0, 0, 1) = 2.0;
    p_start.coeff(0, 1, 0) = 3.0;
    p_start.coeff(1, 0, 0) = 4.0;

    p_start.coeff(0, 0, 2) = 5.0;
    p_start.coeff(0, 1, 1) = 6.0;
    p_start.coeff(0, 2, 0) = 7.0;
    p_start.coeff(1, 0, 1) = 8.0;
    p_start.coeff(1, 1, 0) = 9.0;
    p_start.coeff(2, 0, 0) = 10.0;

    P1 reduced = p_start.change_degree(P1());
    reduced.coeff_iterator(
        [&](const Array<int, dim> &exponents) {
          REQUIRE(reduced.coeff(exponents) ==
                  p_start.coeff(exponents));
        });
  }
  SECTION("Degree 2 to 3") {
    P1 p_start((Tags::Zero_Tag()));
    p_start.coeff(0, 0, 0) = 1.0;

    p_start.coeff(0, 0, 1) = 2.0;
    p_start.coeff(0, 1, 0) = 3.0;
    p_start.coeff(1, 0, 0) = 4.0;

    p_start.coeff(0, 0, 2) = 5.0;
    p_start.coeff(0, 1, 1) = 6.0;
    p_start.coeff(0, 2, 0) = 7.0;
    p_start.coeff(1, 0, 1) = 8.0;
    p_start.coeff(1, 1, 0) = 9.0;
    p_start.coeff(2, 0, 0) = 10.0;

    P2 increased =
        p_start.change_degree(P2((Tags::Zero_Tag())));
    increased.coeff_iterator(
        [&](const Array<int, dim> &exponents) {
          if(CTMath::sum(exponents) <= degree_1) {
            REQUIRE(increased.coeff(exponents) ==
                    p_start.coeff(exponents));
          } else {
            REQUIRE(increased.coeff(exponents) == 0.0);
          }
        });
  }
  SECTION("Degree 2 to 2") {
    P1 p_start((Tags::Zero_Tag()));
    p_start.coeff(0, 0, 0) = 1.0;

    p_start.coeff(0, 0, 1) = 2.0;
    p_start.coeff(0, 1, 0) = 3.0;
    p_start.coeff(1, 0, 0) = 4.0;

    p_start.coeff(0, 0, 2) = 5.0;
    p_start.coeff(0, 1, 1) = 6.0;
    p_start.coeff(0, 2, 0) = 7.0;
    p_start.coeff(1, 0, 1) = 8.0;
    p_start.coeff(1, 1, 0) = 9.0;
    p_start.coeff(2, 0, 0) = 10.0;

    P1 copy = p_start.change_degree(P1((Tags::Zero_Tag())));
    copy.coeff_iterator(
        [&](const Array<int, dim> &exponents) {
          REQUIRE(copy.coeff(exponents) ==
                  p_start.coeff(exponents));
        });
  }
}

TEST_CASE("Polynomial Slice", "[Polynomial]") {
  constexpr const int dim_1 = 3;
  constexpr const int dim_2 = dim_1 - 1;
  constexpr const int degree_range = 2;
  using CoeffT = double;
  using P = Polynomial<CoeffT, degree_range, dim_1>;
  using S = Polynomial<CoeffT, degree_range, dim_2>;
  P x;
  x.coeff(0, 0, 0) = 2.0;

  x.coeff(0, 0, 1) = 3.0;
  x.coeff(0, 1, 0) = 5.0;
  x.coeff(1, 0, 0) = 7.0;

  x.coeff(0, 0, 2) = 1.0;
  x.coeff(0, 1, 1) = 0.0;
  x.coeff(0, 2, 0) = 1.0;
  x.coeff(1, 0, 1) = 0.0;
  x.coeff(1, 1, 0) = 0.0;
  x.coeff(2, 0, 0) = 1.0;

  SECTION("Slice 0") {
    S s = x.slice(0, 2.0);
    REQUIRE(s.coeff(0, 0) == 20.0);

    REQUIRE(s.coeff(0, 1) == 3.0);
    REQUIRE(s.coeff(1, 0) == 5.0);

    REQUIRE(s.coeff(0, 2) == 1.0);
    REQUIRE(s.coeff(1, 1) == 0.0);
    REQUIRE(s.coeff(2, 0) == 1.0);
  }

  SECTION("Slice 1") {
    S s = x.slice(1, 3.0);
    REQUIRE(s.coeff(0, 0) == 26.0);

    REQUIRE(s.coeff(0, 1) == 3.0);
    REQUIRE(s.coeff(1, 0) == 7.0);

    REQUIRE(s.coeff(0, 2) == 1.0);
    REQUIRE(s.coeff(1, 1) == 0.0);
    REQUIRE(s.coeff(2, 0) == 1.0);
  }

  SECTION("Slice 2") {
    S s = x.slice(2, 5.0);
    REQUIRE(s.coeff(0, 0) == 42.0);

    REQUIRE(s.coeff(0, 1) == 5.0);
    REQUIRE(s.coeff(1, 0) == 7.0);

    REQUIRE(s.coeff(0, 2) == 1.0);
    REQUIRE(s.coeff(1, 1) == 0.0);
    REQUIRE(s.coeff(2, 0) == 1.0);
  }
}

TEST_CASE("Polynomial Integral", "[Polynomial]") {
  constexpr const int dim = 3;
  constexpr const int degree_range_1 = 2;
  constexpr const int degree_range_2 = degree_range_1 - 1;
  using CoeffT = double;
  using P = Polynomial<CoeffT, degree_range_1, dim>;
  using D = Polynomial<CoeffT, degree_range_2, dim>;

  P x;
  x.coeff(0, 0, 0) = 2.0;

  x.coeff(0, 0, 1) = 3.0;
  x.coeff(0, 1, 0) = 5.0;
  x.coeff(1, 0, 0) = 7.0;

  x.coeff(0, 0, 2) = 1.0;
  x.coeff(0, 1, 1) = 0.0;
  x.coeff(0, 2, 0) = 1.0;
  x.coeff(1, 0, 1) = 0.0;
  x.coeff(1, 1, 0) = 0.0;
  x.coeff(2, 0, 0) = 1.0;

  SECTION("x0") {
    D y = x.differentiate(0);

    REQUIRE(y.coeff(0, 0, 0) == x.coeff(1, 0, 0));

    REQUIRE(y.coeff(0, 0, 1) == x.coeff(1, 0, 1));
    REQUIRE(y.coeff(0, 1, 0) == x.coeff(1, 1, 0));
    REQUIRE(y.coeff(1, 0, 0) == 2.0 * x.coeff(2, 0, 0));
  }

  SECTION("x1") {
    D y = x.differentiate(1);

    REQUIRE(y.coeff(0, 0, 0) == x.coeff(0, 1, 0));

    REQUIRE(y.coeff(0, 0, 1) == x.coeff(1, 0, 1));
    REQUIRE(y.coeff(0, 1, 0) == 2.0 * x.coeff(0, 2, 0));
    REQUIRE(y.coeff(1, 0, 0) == x.coeff(1, 1, 0));
  }

  SECTION("x2") {
    D y = x.differentiate(2);

    REQUIRE(y.coeff(0, 0, 0) == x.coeff(0, 0, 1));

    REQUIRE(y.coeff(0, 0, 1) == 2.0 * x.coeff(0, 0, 2));
    REQUIRE(y.coeff(0, 1, 0) == x.coeff(0, 1, 1));
    REQUIRE(y.coeff(1, 0, 0) == x.coeff(1, 0, 1));
  }
}

TEST_CASE("Polynomial Sum", "[Polynomial]") {
  constexpr const int dim = 3;
  constexpr const int degree_range_1 = 2;
  constexpr const int degree_range_2 = degree_range_1 + 1;
  using CoeffT = double;
  using P1 = Polynomial<CoeffT, degree_range_1, dim>;
  using P2 = Polynomial<CoeffT, degree_range_2, dim>;
  P1 x;
  x.coeff(0, 0, 0) = 2.0;

  x.coeff(0, 0, 1) = 3.0;
  x.coeff(0, 1, 0) = 5.0;
  x.coeff(1, 0, 0) = 7.0;

  x.coeff(0, 0, 2) = 1.0;
  x.coeff(0, 1, 1) = 0.0;
  x.coeff(0, 2, 0) = 1.0;
  x.coeff(1, 0, 1) = 0.0;
  x.coeff(1, 1, 0) = 0.0;
  x.coeff(2, 0, 0) = 1.0;

  P2 y;
  y.coeff(0, 0, 0) = 3.0;

  y.coeff(0, 0, 1) = 5.0;
  y.coeff(0, 1, 0) = 7.0;
  y.coeff(1, 0, 0) = 11.0;

  y.coeff(0, 0, 2) = 13.0;
  y.coeff(0, 1, 1) = 17.0;
  y.coeff(0, 2, 0) = 19.0;
  y.coeff(1, 0, 1) = 23.0;
  y.coeff(1, 1, 0) = 29.0;
  y.coeff(2, 0, 0) = 31.0;

  y.coeff(0, 0, 3) = 37.0;
  y.coeff(0, 1, 2) = 41.0;
  y.coeff(0, 2, 1) = 43.0;
  y.coeff(0, 3, 0) = 47.0;
  y.coeff(1, 0, 2) = 53.0;
  y.coeff(1, 1, 1) = 59.0;
  y.coeff(1, 2, 0) = 61.0;
  y.coeff(2, 0, 1) = 67.0;
  y.coeff(2, 1, 0) = 71.0;
  y.coeff(3, 0, 0) = 73.0;

  using PSum = Polynomial<CoeffT, degree_range_2, dim>;
  SECTION("dir 1") {
    PSum s = x.sum(y);
    REQUIRE(s.coeff(0, 0, 0) == 5.0);

    REQUIRE(s.coeff(0, 0, 1) == 8.0);
    REQUIRE(s.coeff(0, 1, 0) == 12.0);
    REQUIRE(s.coeff(1, 0, 0) == 18.0);

    REQUIRE(s.coeff(0, 0, 2) == 14.0);
    REQUIRE(s.coeff(0, 1, 1) == 17.0);
    REQUIRE(s.coeff(0, 2, 0) == 20.0);
    REQUIRE(s.coeff(1, 0, 1) == 23.0);
    REQUIRE(s.coeff(1, 1, 0) == 29.0);
    REQUIRE(s.coeff(2, 0, 0) == 32.0);

    REQUIRE(s.coeff(0, 0, 3) == 37.0);
    REQUIRE(s.coeff(0, 1, 2) == 41.0);
    REQUIRE(s.coeff(0, 2, 1) == 43.0);
    REQUIRE(s.coeff(0, 3, 0) == 47.0);
    REQUIRE(s.coeff(1, 0, 2) == 53.0);
    REQUIRE(s.coeff(1, 1, 1) == 59.0);
    REQUIRE(s.coeff(1, 2, 0) == 61.0);
    REQUIRE(s.coeff(2, 0, 1) == 67.0);
    REQUIRE(s.coeff(2, 1, 0) == 71.0);
    REQUIRE(s.coeff(3, 0, 0) == 73.0);
  }
  SECTION("dir 2") {
    PSum s = y.sum(x);
    REQUIRE(s.coeff(0, 0, 0) == 5.0);

    REQUIRE(s.coeff(0, 0, 1) == 8.0);
    REQUIRE(s.coeff(0, 1, 0) == 12.0);
    REQUIRE(s.coeff(1, 0, 0) == 18.0);

    REQUIRE(s.coeff(0, 0, 2) == 14.0);
    REQUIRE(s.coeff(0, 1, 1) == 17.0);
    REQUIRE(s.coeff(0, 2, 0) == 20.0);
    REQUIRE(s.coeff(1, 0, 1) == 23.0);
    REQUIRE(s.coeff(1, 1, 0) == 29.0);
    REQUIRE(s.coeff(2, 0, 0) == 32.0);

    REQUIRE(s.coeff(0, 0, 3) == 37.0);
    REQUIRE(s.coeff(0, 1, 2) == 41.0);
    REQUIRE(s.coeff(0, 2, 1) == 43.0);
    REQUIRE(s.coeff(0, 3, 0) == 47.0);
    REQUIRE(s.coeff(1, 0, 2) == 53.0);
    REQUIRE(s.coeff(1, 1, 1) == 59.0);
    REQUIRE(s.coeff(1, 2, 0) == 61.0);
    REQUIRE(s.coeff(2, 0, 1) == 67.0);
    REQUIRE(s.coeff(2, 1, 0) == 71.0);
    REQUIRE(s.coeff(3, 0, 0) == 73.0);
  }
}

TEST_CASE("Polynomial Product", "[Polynomial]") {
  constexpr const int dim = 3;
  constexpr const int degree_range_1 = 2;
  constexpr const int degree_range_2 = 2;
  using CoeffT = double;
  using P1 = Polynomial<CoeffT, degree_range_1, dim>;
  using P2 = Polynomial<CoeffT, degree_range_2, dim>;
  P1 x;
  x.coeff(0, 0, 0) = 2.0;

  x.coeff(0, 0, 1) = 3.0;
  x.coeff(0, 1, 0) = 5.0;
  x.coeff(1, 0, 0) = 7.0;

  x.coeff(0, 0, 2) = 1.0;
  x.coeff(0, 1, 1) = 0.0;
  x.coeff(0, 2, 0) = 1.0;
  x.coeff(1, 0, 1) = 0.0;
  x.coeff(1, 1, 0) = 0.0;
  x.coeff(2, 0, 0) = 1.0;

  P2 y;
  y.coeff(0, 0, 0) = 2.0;

  y.coeff(0, 0, 1) = 3.0;
  y.coeff(0, 1, 0) = 5.0;
  y.coeff(1, 0, 0) = 7.0;

  y.coeff(0, 0, 2) = 1.0;
  y.coeff(0, 1, 1) = 0.0;
  y.coeff(0, 2, 0) = 1.0;
  y.coeff(1, 0, 1) = 0.0;
  y.coeff(1, 1, 0) = 0.0;
  y.coeff(2, 0, 0) = 1.0;

  using PProd =
      Polynomial<CoeffT, degree_range_1 + degree_range_2,
                 dim>;
  PProd p = x.product(y);
  REQUIRE(p.coeff(0, 0, 0) == 4.0);

  REQUIRE(p.coeff(0, 0, 1) == 12.0);
  REQUIRE(p.coeff(0, 1, 0) == 20.0);
  REQUIRE(p.coeff(1, 0, 0) == 28.0);

  REQUIRE(p.coeff(0, 0, 2) == 13.0);
  REQUIRE(p.coeff(0, 1, 1) == 30.0);
  REQUIRE(p.coeff(0, 2, 0) == 29.0);
  REQUIRE(p.coeff(1, 0, 1) == 42.0);
  REQUIRE(p.coeff(1, 1, 0) == 70.0);
  REQUIRE(p.coeff(2, 0, 0) == 53.0);

  REQUIRE(p.coeff(0, 0, 3) == 6.0);
  REQUIRE(p.coeff(0, 1, 2) == 10.0);
  REQUIRE(p.coeff(0, 2, 1) == 6.0);
  REQUIRE(p.coeff(0, 3, 0) == 10.0);
  REQUIRE(p.coeff(1, 0, 2) == 14.0);
  REQUIRE(p.coeff(1, 1, 1) == 0.0);
  REQUIRE(p.coeff(1, 2, 0) == 14.0);
  REQUIRE(p.coeff(2, 0, 1) == 6.0);
  REQUIRE(p.coeff(2, 1, 0) == 10.0);
  REQUIRE(p.coeff(3, 0, 0) == 14.0);

  REQUIRE(p.coeff(0, 0, 4) == 1.0);
  REQUIRE(p.coeff(0, 1, 3) == 0.0);
  REQUIRE(p.coeff(0, 2, 2) == 2.0);
  REQUIRE(p.coeff(0, 3, 1) == 0.0);
  REQUIRE(p.coeff(0, 4, 0) == 1.0);
  REQUIRE(p.coeff(1, 0, 3) == 0.0);
  REQUIRE(p.coeff(1, 1, 2) == 0.0);
  REQUIRE(p.coeff(1, 2, 1) == 0.0);
  REQUIRE(p.coeff(1, 3, 0) == 0.0);
  REQUIRE(p.coeff(2, 0, 2) == 2.0);
  REQUIRE(p.coeff(2, 1, 1) == 0.0);
  REQUIRE(p.coeff(2, 2, 0) == 2.0);
  REQUIRE(p.coeff(3, 0, 1) == 0.0);
  REQUIRE(p.coeff(3, 1, 0) == 00.0);
  REQUIRE(p.coeff(4, 0, 0) == 1.0);
}

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

  p.coeff(0, 0, 2) = 1.0;
  p.coeff(0, 1, 1) = 0.0;
  p.coeff(0, 2, 0) = 1.0;
  p.coeff(1, 0, 1) = 0.0;
  p.coeff(1, 1, 0) = 0.0;
  p.coeff(2, 0, 0) = 1.0;

  REQUIRE(p.eval(0.0, 0.0, 0.0) == 2.0);

  REQUIRE(p.eval(0.0, 0.0, 1.0) == 6.0);
  REQUIRE(p.eval(0.0, 1.0, 0.0) == 8.0);
  REQUIRE(p.eval(1.0, 0.0, 0.0) == 10.0);
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
