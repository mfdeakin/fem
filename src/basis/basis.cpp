
#include <iostream>
#include <cmath>

#include "polynomial.hpp"

constexpr const int dim = 3;
using CoeffT = double;

template <typename P1, typename P2>
CoeffT dot_product(const P1 &x, const P2 &y) {
  auto idp =
      x.product(y).integrate(0).integrate(1).integrate(2);

  return idp.eval(1.0, 1.0, 1.0) - idp.eval(0.0, 0.0, 0.0);
}

template <typename real, typename integer>
real subfactorial(const integer &n) noexcept {
  /* Computes !n */
  real sum = 0.0;
  for(int k = 0; k <= n; k++) {
    real v =
        real(1.0) / real(CTMath::partialFactorial(1, k));
    if(k % 2 == 0) {
      sum += v;
    } else {
      sum -= v;
    }
  }
  real n_fact = real(CTMath::partialFactorial(1, n));
  return n_fact * sum;
}

CoeffT x_exp_integral(int exp) {
  constexpr const CoeffT E =
      2.7182818284590452353602874713526624977572470937;
  CoeffT v = E * subfactorial<CoeffT, int>(exp) -
             CTMath::partialFactorial(1, exp);
  if(exp % 2 == 1) {
    v = -v;
  }
  return v;
}

template <typename Poly>
CoeffT exp_dot_product(const Poly &p) {
  /* Computes the dot product of p with
   * exp(x1 + x2 + ... + xn)
   * Each integral is separable, so compute each dimension
   * separately then multiply them together
   * \int_0^1 x^n exp(x) dx = (-1)^n (e(!n)-n!)
   */
  CoeffT dp = 0.0;
  p.coeff_iterator(
      [&](const Array<int, Poly::dim> &exponents) {
        CoeffT term = p.coeff(exponents);
        for(int d = 0; d < Poly::dim; d++) {
          term *= x_exp_integral(exponents[d]);
        }
        dp += term;
      });
  return dp;
}

template <typename P1, typename P2>
P2 project(const P1 &p, const P2 &dir) {
  return (dot_product(p, dir) / dot_product(dir, dir)) *
         dir;
}

template <typename P1, typename P2>
P1 orthogonal_helper(const P1 &p, const P2 &next) {
  P2 projected =
      (dot_product(p, next) / dot_product(next, next)) *
      next;
  P1 reduced =
      projected.change_degree(P1((Tags::Zero_Tag())));
  return reduced;
}

template <typename P1, typename P2, typename... P_Prior>
P1 orthogonal_helper(const P1 &p, const P2 &next,
                     P_Prior... prior_basis) {
  const P1 prior_projection =
      orthogonal_helper(p, prior_basis...);
  const P2 projected =
      (dot_product(p, next) / dot_product(next, next)) *
      next;
  P1 reduced =
      projected.change_degree(P1((Tags::Zero_Tag())));
  return reduced + prior_projection;
}

template <typename P1, typename... P_Prior>
P1 orthogonal(const P1 &p, P_Prior... prior_basis) {
  P1 projected = orthogonal_helper(p, prior_basis...);
  P1 unscaled = p + (-projected);
  return unscaled *
         std::sqrt(1.0 / dot_product(unscaled, unscaled));
}

int main(int argc, char **argv) {
  Numerical::Polynomial<CoeffT, 0, dim> constant;
  constant.coeff(0, 0, 0) = 1.0;
  using linear_p = Numerical::Polynomial<CoeffT, 1, dim>;
  linear_p linearx((Tags::Zero_Tag()));
  linear_p lineary((Tags::Zero_Tag()));
  linear_p linearz((Tags::Zero_Tag()));
  linearx.coeff(0, 0, 1) = 1.0;
  lineary.coeff(0, 1, 0) = 1.0;
  linearz.coeff(1, 0, 0) = 1.0;

  linear_p lx_o_c = orthogonal(linearx, constant);
  linear_p ly_o_c = orthogonal(lineary, constant, lx_o_c);
  linear_p lz_o_c =
      orthogonal(linearz, constant, lx_o_c, ly_o_c);

  using quadratic_p = Numerical::Polynomial<CoeffT, 2, dim>;
  quadratic_p quad_xx((Tags::Zero_Tag()));
  quad_xx.coeff(2, 0, 0) = 1.0;
  quadratic_p quad_o_xx =
      orthogonal(quad_xx, constant, lx_o_c, ly_o_c, lz_o_c);

  quadratic_p quad_xy((Tags::Zero_Tag()));
  quad_xy.coeff(1, 1, 0) = 1.0;
  quadratic_p quad_o_xy = orthogonal(
      quad_xy, constant, lx_o_c, ly_o_c, lz_o_c, quad_o_xx);

  quadratic_p quad_xz((Tags::Zero_Tag()));
  quad_xz.coeff(1, 0, 1) = 1.0;
  quadratic_p quad_o_xz =
      orthogonal(quad_xz, constant, lx_o_c, ly_o_c, lz_o_c,
                 quad_o_xx, quad_o_xy);

  quadratic_p quad_yy((Tags::Zero_Tag()));
  quad_yy.coeff(0, 2, 0) = 1.0;
  quadratic_p quad_o_yy =
      orthogonal(quad_yy, constant, lx_o_c, ly_o_c, lz_o_c,
                 quad_o_xx, quad_o_xy, quad_o_xz);

  quadratic_p quad_yz((Tags::Zero_Tag()));
  quad_yz.coeff(0, 1, 1) = 1.0;
  quadratic_p quad_o_yz = orthogonal(
      quad_yz, constant, lx_o_c, ly_o_c, lz_o_c, quad_o_xx,
      quad_o_xy, quad_o_xz, quad_o_yy);

  quadratic_p quad_zz((Tags::Zero_Tag()));
  quad_zz.coeff(0, 0, 2) = 1.0;
  quadratic_p quad_o_zz = orthogonal(
      quad_zz, constant, lx_o_c, ly_o_c, lz_o_c, quad_o_xx,
      quad_o_xy, quad_o_xz, quad_o_yy, quad_o_yz);

  using cubic_p = Numerical::Polynomial<CoeffT, 3, dim>;
  
  quadratic_p exp_proj =
      exp_dot_product(constant) * constant +
      exp_dot_product(lx_o_c) * lx_o_c +
      exp_dot_product(ly_o_c) * ly_o_c +
      exp_dot_product(lz_o_c) * lz_o_c +
      exp_dot_product(quad_o_xx) * quad_o_xx +
      exp_dot_product(quad_o_xy) * quad_o_xy +
      exp_dot_product(quad_o_xz) * quad_o_xz +
      exp_dot_product(quad_o_yy) * quad_o_yy +
      exp_dot_product(quad_o_yz) * quad_o_yz +
      exp_dot_product(quad_o_zz) * quad_o_zz;

  constexpr CoeffT optimal_dp =
      32.6001889612617945069684599145663334796335791191178;
  std::cout << "Exponent Projection: "
            << exp_dot_product(exp_proj) << " vs "
            << optimal_dp << "; "
            << exp_dot_product(exp_proj + 0.001 * lx_o_c)
            << "; "
            << exp_dot_product(exp_proj + -0.001 * lx_o_c)
            << "; " << std::endl
            << dot_product(exp_proj, exp_proj) << std::endl;
  exp_proj.coeff_iterator(
      [&](const Array<int, dim> &exponents) {
        std::cout << exponents << " : "
                  << exp_proj.coeff(exponents) << std::endl;
      });

  return 0;
}
