
#include <iostream>

#include "polynomial.hpp"

constexpr const int dim = 3;
using CoeffT = double;

template <typename P1, typename P2>
CoeffT dot_product(const P1 &x, const P2 &y) {
  auto idp =
      x.product(y).integrate(0).integrate(1).integrate(2);

  return idp.eval(1.0, 1.0, 1.0) - idp.eval(0.0, 0.0, 0.0);
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
  return p + (-projected);
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

  std::cout << "c . c : " << dot_product(constant, constant)
            << std::endl;

  std::cout << "lx_o . c : "
            << dot_product(lx_o_c, constant) << std::endl;
  std::cout << "lx_o . lx_o : "
            << dot_product(lx_o_c, lx_o_c) << std::endl;
  std::cout << "ly_o . c : "
            << dot_product(ly_o_c, constant) << std::endl;
  std::cout << "ly_o . lx_o : "
            << dot_product(ly_o_c, lx_o_c) << std::endl;
  std::cout << "ly_o . ly_o : "
            << dot_product(ly_o_c, ly_o_c) << std::endl;
  std::cout << "lz_o . c : "
            << dot_product(lz_o_c, constant) << std::endl;
  std::cout << "lz_o . lx_o : "
            << dot_product(lz_o_c, lx_o_c) << std::endl;
  std::cout << "lz_o . ly_o : "
            << dot_product(lz_o_c, ly_o_c) << std::endl;
  std::cout << "lz_o . lz_o : "
            << dot_product(lz_o_c, lz_o_c) << std::endl;

  std::cout << "qxx_o . c : "
            << dot_product(quad_o_xx, constant)
            << std::endl;
  std::cout << "qxx . lx_o : "
            << dot_product(quad_o_xx, lx_o_c) << std::endl;
  std::cout << "qxx_o . ly_o : "
            << dot_product(quad_o_xx, ly_o_c) << std::endl;
  std::cout << "qxx_o . lz_o : "
            << dot_product(quad_o_xx, lz_o_c) << std::endl;
  std::cout << "qxx_o . qxx_o : "
            << dot_product(quad_o_xx, quad_o_xx)
            << std::endl;

  std::cout << "qxy_o . c : "
            << dot_product(quad_o_xy, constant)
            << std::endl;
  std::cout << "qxy . lx_o : "
            << dot_product(quad_o_xy, lx_o_c) << std::endl;
  std::cout << "qxy_o . ly_o : "
            << dot_product(quad_o_xy, ly_o_c) << std::endl;
  std::cout << "qxy_o . lz_o : "
            << dot_product(quad_o_xy, lz_o_c) << std::endl;
  std::cout << "qxy_o . qxx_o : "
            << dot_product(quad_o_xy, quad_o_xx)
            << std::endl;
  std::cout << "qxy_o . qxy_o : "
            << dot_product(quad_o_xy, quad_o_xy)
            << std::endl;

  std::cout << "qxz_o . c : "
            << dot_product(quad_o_xz, constant)
            << std::endl;
  std::cout << "qxz . lx_o : "
            << dot_product(quad_o_xz, lx_o_c) << std::endl;
  std::cout << "qxz_o . ly_o : "
            << dot_product(quad_o_xz, ly_o_c) << std::endl;
  std::cout << "qxz_o . lz_o : "
            << dot_product(quad_o_xz, lz_o_c) << std::endl;
  std::cout << "qxz_o . qxx_o : "
            << dot_product(quad_o_xz, quad_o_xx)
            << std::endl;
  std::cout << "qxz_o . qxy_o : "
            << dot_product(quad_o_xz, quad_o_xy)
            << std::endl;
  std::cout << "qxz_o . qxz_o : "
            << dot_product(quad_o_xz, quad_o_xz)
            << std::endl;

  std::cout << "qyy_o . c : "
            << dot_product(quad_o_yy, constant)
            << std::endl;
  std::cout << "qyy . lx_o : "
            << dot_product(quad_o_yy, lx_o_c) << std::endl;
  std::cout << "qyy_o . ly_o : "
            << dot_product(quad_o_yy, ly_o_c) << std::endl;
  std::cout << "qyy_o . lz_o : "
            << dot_product(quad_o_yy, lz_o_c) << std::endl;
  std::cout << "qyy_o . qxx_o : "
            << dot_product(quad_o_yy, quad_o_xx)
            << std::endl;
  std::cout << "qyy_o . qxy_o : "
            << dot_product(quad_o_yy, quad_o_xy)
            << std::endl;
  std::cout << "qyy_o . qxz_o : "
            << dot_product(quad_o_yy, quad_o_xz)
            << std::endl;
  std::cout << "qyy_o . qyy_o : "
            << dot_product(quad_o_yy, quad_o_yy)
            << std::endl;

  std::cout << "qyz_o . c : "
            << dot_product(quad_o_yz, constant)
            << std::endl;
  std::cout << "qyz . lx_o : "
            << dot_product(quad_o_yz, lx_o_c) << std::endl;
  std::cout << "qyz_o . ly_o : "
            << dot_product(quad_o_yz, ly_o_c) << std::endl;
  std::cout << "qyz_o . lz_o : "
            << dot_product(quad_o_yz, lz_o_c) << std::endl;
  std::cout << "qyz_o . qxx_o : "
            << dot_product(quad_o_yz, quad_o_xx)
            << std::endl;
  std::cout << "qyz_o . qxy_o : "
            << dot_product(quad_o_yz, quad_o_xy)
            << std::endl;
  std::cout << "qyz_o . qxz_o : "
            << dot_product(quad_o_yz, quad_o_xz)
            << std::endl;
  std::cout << "qyz_o . qyy_o : "
            << dot_product(quad_o_yz, quad_o_yy)
            << std::endl;
  std::cout << "qyz_o . qyz_o : "
            << dot_product(quad_o_yz, quad_o_yz)
            << std::endl;

  std::cout << "qzz_o . c : "
            << dot_product(quad_o_zz, constant)
            << std::endl;
  std::cout << "qzz . lx_o : "
            << dot_product(quad_o_zz, lx_o_c) << std::endl;
  std::cout << "qzz_o . ly_o : "
            << dot_product(quad_o_zz, ly_o_c) << std::endl;
  std::cout << "qzz_o . lz_o : "
            << dot_product(quad_o_zz, lz_o_c) << std::endl;
  std::cout << "qzz_o . qxx_o : "
            << dot_product(quad_o_zz, quad_o_xx)
            << std::endl;
  std::cout << "qzz_o . qxy_o : "
            << dot_product(quad_o_zz, quad_o_xy)
            << std::endl;
  std::cout << "qzz_o . qxz_o : "
            << dot_product(quad_o_zz, quad_o_xz)
            << std::endl;
  std::cout << "qzz_o . qyy_o : "
            << dot_product(quad_o_zz, quad_o_yy)
            << std::endl;
  std::cout << "qzz_o . qyz_o : "
            << dot_product(quad_o_zz, quad_o_yz)
            << std::endl;
  std::cout << "qzz_o . qzz_o : "
            << dot_product(quad_o_zz, quad_o_zz)
            << std::endl;

  return 0;
}
