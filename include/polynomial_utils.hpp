
#ifndef _POLYNOMIAL_UTILS_HPP_
#define _POLYNOMIAL_UTILS_HPP_

#include <tuple>

#include <ctmath.hpp>

namespace Numerical {

template <typename CoeffT, int _degree, int _dim>
class Polynomial;

namespace Utilities {

// Returns the total number of coefficients required to
// represent polynomials of the specified degree
template <typename int_t>
constexpr int_t poly_num_coeffs(int_t degree_range,
                                int_t dim) noexcept {
  return CTMath::n_choose_k<int_t>(dim + degree_range, dim);
}

// Returns the total number of coefficients required to
// represent polynomials only containing terms of the
// specified degree
template <typename int_t>
constexpr int_t poly_degree_num_coeffs(int_t degree,
                                       int_t dim) noexcept {
  return CTMath::n_choose_k<int_t>(dim + degree - 1,
                                   degree);
}

// A TMP for deducing the tuple type required to represent a
// basis of the specified degree
// Starts with degree 0 and continues in increasing order
template <typename CoeffT, int _max_degree, int _dim,
          int _index = 1,
          int _remaining_degree = _max_degree,
          typename... poly_pack>
class basis_tuple {
 public:
  static constexpr int current_degree =
      _max_degree - _remaining_degree;
  using poly_type =
      Polynomial<CoeffT, current_degree, _dim>;

  using next =
      basis_tuple<CoeffT, _max_degree, _dim, _index + 1,
                  _index >= poly_num_coeffs(current_degree,
                                            _dim)
                      ? _remaining_degree - 1
                      : _remaining_degree,
                  poly_pack..., poly_type>;

  using tuple_type = typename next::tuple_type;
};

template <typename CoeffT, int _max_degree, int _dim,
          int _index, typename... poly_pack>
class basis_tuple<CoeffT, _max_degree, _dim, _index, -1,
                  poly_pack...> {
 public:
  using poly_type = void;
  using next = void;
  using tuple_type = std::tuple<poly_pack...>;
  using tuple_ptr_type = std::tuple<poly_pack...>;
};

// Maps an integer index to the exponents
// 0 -> (degree, 0, 0, ...)
// 1 -> (degree - 1, 1, 0, 0, ...)
// 2 -> (degree - 1, 0, 1, 0, 0, ...)
// n -> (degree - 1, 0, 0, ..., 1)
// n + 1 -> (degree - 2, 2, 0, 0, ...)
// n + 2 -> (degree - 2, 1, 1, 0, 0, ...)
// n + 3 -> (degree - 2, 1, 0, 1, 0, ...)
// n + m -> (degree - 2, 1, 0, 0, ..., 1)
template <int dim>
void index_to_exponents(const int deg_idx, const int degree,
                        const int coeff,
                        Array<int, dim> &exponents) {
  assert(deg_idx >= 0);
  assert(degree >= 0);
  assert(coeff >= 0);
  assert((degree == 0) ? (deg_idx == 0) : true);
  if(deg_idx > 0) {
    for(exponents[coeff] = degree;
        deg_idx >=
            poly_num_coeffs(degree - exponents[coeff],
                            dim - coeff - 1) &&
        exponents[coeff] > 0;
        exponents[coeff]--) {
    }
    index_to_exponents(
        std::max(
            0, deg_idx - poly_num_coeffs(
                             degree - exponents[coeff] - 1,
                             dim - coeff - 1)),
        degree - exponents[coeff], coeff + 1, exponents);
  } else {
    exponents[coeff] = degree;
    for(int c_idx = coeff + 1; c_idx < dim; c_idx++) {
      exponents[c_idx] = 0.0;
    }
  }
}

template <typename CoeffT, int _max_degree, int _dim,
          int _cur_degree = 0, int _index = 1>
class BasisGenerators {
 public:
  static void unit_basis(
      typename basis_tuple<CoeffT, _max_degree,
                           _dim>::tuple_type &basis,
      Array<int, _dim> &exponents) {
    Polynomial<CoeffT, _cur_degree, _dim> &current =
        std::get<_index - 1>(basis);

    current = Polynomial<CoeffT, _cur_degree, _dim>(
        Tags::Zero_Tag());
    int deg_index =
        _cur_degree > 0
            ? std::max(0, _index -
                              poly_num_coeffs(
                                  _cur_degree - 1, _dim) -
                              1)

            : 0;
    index_to_exponents(deg_index, _cur_degree, 0,
                       exponents);
    current.coeff(exponents) = 1.0;

    BasisGenerators<CoeffT, _max_degree, _dim,
                    _index >= poly_num_coeffs(_cur_degree,
                                              _dim)
                        ? _cur_degree + 1
                        : _cur_degree,
                    _index + 1>::unit_basis(basis,
                                            exponents);
  }
};

template <typename CoeffT, int _max_degree, int _dim,
          int _index>
class BasisGenerators<CoeffT, _max_degree, _dim,
                      _max_degree, _index> {
 public:
  static void unit_basis(
      typename basis_tuple<CoeffT, _max_degree,
                           _dim>::tuple_type &basis,
      Array<int, _dim> &exponents) {}
};
}  // namespace Utilities
}  // namespace Numerical

#endif
