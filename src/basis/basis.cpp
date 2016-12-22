
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

int main(int argc, char **argv) {
  using linear_p = Numerical::Polynomial<CoeffT, 1, dim>;
  Numerical::Polynomial<CoeffT, 0, dim> constant;
  constant.coeff(0, 0, 0) = 1.0;
  linear_p linearx((Tags::Zero_Tag()));
  linear_p lineary((Tags::Zero_Tag()));
  linear_p linearz((Tags::Zero_Tag()));
  linearx.coeff(0, 0, 0) = 1.0;
  lineary.coeff(0, 0, 0) = 1.0;
  linearz.coeff(0, 0, 0) = 1.0;

  linearx.coeff(0, 0, 1) = 1.0;
  lineary.coeff(0, 1, 0) = 1.0;
  linearz.coeff(1, 0, 0) = 1.0;

  linear_p lx_o_c = linearx +
                    -((dot_product(linearx, constant) /
                       dot_product(constant, constant)) *
                      constant);

  linear_p ly_o_c = lineary +
                    -((dot_product(lineary, lx_o_c) /
                       dot_product(lx_o_c, lx_o_c)) *
                          lx_o_c +
                      (dot_product(lineary, constant) /
                       dot_product(constant, constant)) *
                          constant);

  linear_p lz_o_c = linearz +
                    -((dot_product(linearz, lx_o_c) /
                       dot_product(lx_o_c, lx_o_c)) *
                          lx_o_c +
                      (dot_product(linearz, ly_o_c) /
                       dot_product(ly_o_c, ly_o_c)) *
                          ly_o_c +
                      (dot_product(linearz, constant) /
                       dot_product(constant, constant)) *
                          constant);

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

  return 0;
}
