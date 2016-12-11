
#ifndef _POLYNOMIAL_HPP_
#define _POLYNOMIAL_HPP_

#include "array.hpp"
#include "ctmath.hpp"

#include <iostream>

#include <type_traits>

namespace Numerical {

namespace {
struct Zero_Tag {};
}

template <typename CoeffT, int _degree_range, int _dim>
class Polynomial {
 public:
  Polynomial() {
    static_assert(_degree_range >= 0,
                  "A polynomial's _degree_range (max "
                  "exponent-min exponent) must be at least "
                  "zero, otherwise it's degenerate");
    static_assert(_dim >= 0,
                  "A polynomial's _dimension must be at "
                  "least zero, otherwise it's degenerate");
  }

  Polynomial(const Zero_Tag &) : lower_degree(Zero_Tag()) {
    static_assert(_degree_range >= 0,
                  "A polynomial's _degree_range (max "
                  "exponent-min exponent) must be at least "
                  "zero, otherwise it's degenerate");
    static_assert(_dim >= 0,
                  "A polynomial's _dimension must be at "
                  "least zero, otherwise it's degenerate");
    for(int i = 0; i < num_coeffs; ++i) {
      coeffs[i] = 0;
    }
  }

  template <typename... int_list,
            typename std::enable_if<
                sizeof...(int_list) == _dim, int>::type = 0>
  CoeffT coeff(int_list... args) const noexcept {
    if(CTMath::sum(args...) == _degree_range) {
      int idx =
          get_coeff_idx_helper(_degree_range, args...);
      assert(idx >= 0);
      assert(idx < num_coeffs);
      return coeffs[idx];
    } else {
      return lower_degree.coeff(args...);
    }
  }

  template <typename... int_list,
            typename std::enable_if<
                sizeof...(int_list) == _dim, int>::type = 0>
  CoeffT &coeff(int_list... args) noexcept {
    if(CTMath::sum(args...) == _degree_range) {
      int idx =
          get_coeff_idx_helper(_degree_range, args...);
      assert(idx >= 0);
      assert(idx < num_coeffs);
      return coeffs[idx];
    } else {
      return lower_degree.coeff(args...);
    }
  }

  CoeffT coeff(const Array<int, _dim> &exponents) const
      noexcept {
    if(CTMath::sum(exponents) == _degree_range) {
      int idx =
          get_coeff_idx_helper(_degree_range, exponents);
      assert(idx >= 0);
      assert(idx < num_coeffs);
      return coeffs[idx];
    } else {
      return lower_degree.coeff(exponents);
    }
  }

  CoeffT &coeff(
      const Array<int, _dim> &exponents) noexcept {
    if(CTMath::sum(exponents) == _degree_range) {
      int idx =
          get_coeff_idx_helper(_degree_range, exponents);
      assert(idx >= 0);
      assert(idx < num_coeffs);
      return coeffs[idx];
    } else {
      return lower_degree.coeff(exponents);
    }
  }

  template <int other_degree_range>
  Polynomial<CoeffT, _degree_range + other_degree_range,
             _dim>
  product(const Polynomial<CoeffT, other_degree_range, _dim>
              &m) const {
    using FP = Polynomial<
        CoeffT, _degree_range + other_degree_range, _dim>;
    FP prod((Zero_Tag()));
    Array<int, _dim> buf;
    coeff_iterator(
        _degree_range, 0, buf,
        [&](const Array<int, _dim> &exponents) {
          Array<int, _dim> buf2;
          m.coeff_iterator(
              other_degree_range, 0, buf2,
              [&](const Array<int, _dim> &other_exponents) {
                Array<int, _dim> final_exponents;
                for(int i = 0; i < _dim; i++) {
                  final_exponents[i] =
                      exponents[i] + other_exponents[i];
                }
                prod.coeff(final_exponents) +=
                    coeff(exponents) *
                    m.coeff(other_exponents);
              });
        });
    return prod;
  }

  template <
      typename... subs_list,
      typename std::enable_if<sizeof...(subs_list) == _dim,
                              int>::type = 0>
  CoeffT eval(subs_list... vars) const {
    Array<int, _dim> exponents;
    return eval_helper(_degree_range, exponents, vars...);
  }

  template <typename, int, int>
  friend class Polynomial;

 private:
  template <typename Lambda>
  void coeff_iterator(const int exp_left, const int cur_dim,
                      Array<int, _dim> &exponents,
                      Lambda function) const {
    for(exponents[cur_dim] = 0;
        exponents[cur_dim] <= exp_left;
        exponents[cur_dim]++) {
      if(cur_dim == _dim - 1) {
        function(exponents);
      } else {
        coeff_iterator(exp_left - exponents[cur_dim],
                       cur_dim + 1, exponents, function);
      }
    }
  }

  template <typename... subs_list>
  CoeffT eval_helper(int exp_left,
                     Array<int, _dim> &exponents,
                     CoeffT cur_var,
                     subs_list... vars) const {
    constexpr const auto cur_dim =
        _dim - sizeof...(subs_list)-1;
    CoeffT factor = 1.0;
    CoeffT term_sum = 0.0;
    for(exponents[cur_dim] = 0;
        exponents[cur_dim] <= exp_left;
        ++exponents[cur_dim]) {
      term_sum += factor *
                  eval_helper(exp_left - exponents[cur_dim],
                              exponents, vars...);
      factor *= cur_var;
    }
    return term_sum;
  }

  CoeffT eval_helper(int exp_left,
                     Array<int, _dim> &exponents,
                     CoeffT cur_var) const noexcept {
    constexpr const auto cur_dim = _dim - 1;
    CoeffT term_sum = 0.0;
    CoeffT factor = 1.0;
    for(exponents[cur_dim] = 0;
        exponents[cur_dim] <= exp_left;
        ++exponents[cur_dim]) {
      term_sum += factor * coeff(exponents);
      factor *= cur_var;
    }
    return term_sum;
  }

  template <typename... int_list,
            typename std::enable_if<
                sizeof...(int_list) == _dim, int>::type = 0>
  static int get_coeff_idx(int_list... args) noexcept {
    const int term_degree = CTMath::sum(args...);
    if(term_degree >= 0 && term_degree <= _degree_range) {
      return get_coeff_idx_helper(term_degree, args...);
    } else {
      return DEGREE_RANGE_OVERFLOW;
    }
  }

  template <typename int_t, typename... int_list>
  static int get_coeff_idx_helper(
      int_t exp_left, int_t head,
      int_list... args) noexcept {
    int nck = _dim + exp_left - 1;
    int term_1 =
        nck * CTMath::n_choose_k(nck - 1, exp_left);
    int term_2 =
        (nck - head) *
        CTMath::n_choose_k(nck - head - 1, exp_left - head);
    int first_index = (term_1 - term_2) / (_dim - 1);
    return first_index +
           Polynomial<CoeffT, _degree_range, _dim - 1>::
               get_coeff_idx_helper(exp_left - head,
                                    args...);
  }

  template <typename int_t>
  static int get_coeff_idx_helper(int_t exp_left,
                                  int_t arg) noexcept {
    return arg == exp_left
               ? 0
               : arg > exp_left ? DEGREE_RANGE_OVERFLOW
                                : DEGREE_RANGE_UNDERFLOW;
  }

  template <typename int_t>
  static int get_coeff_idx_helper(
      int_t exp_left,
      const Array<int_t, _dim> &exponents) noexcept {
    int cur_sum = 0;
    for(int i = 0; i < _dim - 1; ++i) {
      int nck = _dim - i + exp_left - 1;
      int term_1 =
          nck * CTMath::n_choose_k(nck - 1, exp_left);
      exp_left -= exponents[i];
      int term_2 = (nck - exponents[i]) *
                   CTMath::n_choose_k(
                       nck - exponents[i] - 1, exp_left);
      cur_sum += (term_1 - term_2) / (_dim - i - 1);
    }
    return cur_sum;
  }

  static constexpr const int dim = _dim;

  static constexpr const int num_coeffs =
      CTMath::poly_degree_num_coeffs<int>(_degree_range,
                                          _dim);
  Array<CoeffT, num_coeffs> coeffs;

  Polynomial<CoeffT, _degree_range - 1, _dim> lower_degree;

  enum {
    DEGREE_RANGE_UNDERFLOW = -1,
    DEGREE_RANGE_OVERFLOW = -2,
  };
};

template <typename CoeffT, int _dim>
class Polynomial<CoeffT, 0, _dim> {
 public:
  Polynomial() {}

  Polynomial(const Zero_Tag &) : value(0) {}

  template <typename... int_list,
            typename std::enable_if<
                sizeof...(int_list) == _dim, int>::type = 0>
  CoeffT coeff(int_list... args) const noexcept {
    assert(CTMath::sum(args...) == 0);
    return value;
  }

  template <typename... int_list,
            typename std::enable_if<
                sizeof...(int_list) == _dim, int>::type = 0>
  CoeffT &coeff(int_list... args) noexcept {
    assert(CTMath::sum(args...) == 0);
    return value;
  }

  CoeffT coeff(const Array<int, _dim> &exponents) const
      noexcept {
    assert(CTMath::sum(exponents) == 0);
    return value;
  }

  CoeffT &coeff(
      const Array<int, _dim> &exponents) noexcept {
    assert(CTMath::sum(exponents) == 0);
    return value;
  }

  template <typename... int_list,
            typename std::enable_if<
                sizeof...(int_list) == _dim, int>::type = 0>
  static int get_coeff_idx(int_list... args) noexcept {
    return CTMath::sum(args...) == 0
               ? 0
               : CTMath::sum(args...) > 0
                     ? DEGREE_RANGE_OVERFLOW
                     : DEGREE_RANGE_UNDERFLOW;
  }

  template <typename int_t, typename... int_list>
  static int get_coeff_idx_helper(
      int_t exp_left, int_t head,
      int_list... args) noexcept {
    return get_coeff_idx(head, args...);
  }

  template <typename, int, int>
  friend class Polynomial;

 private:
  static constexpr const int dim = _dim;

  static constexpr const int num_coeffs = 1;
  CoeffT value;

  enum {
    DEGREE_RANGE_UNDERFLOW = -1,
    DEGREE_RANGE_OVERFLOW = -2,
  };
};

template <typename CoeffT, int _degree_range>
class Polynomial<CoeffT, _degree_range, 0> {
 public:
  Polynomial() {
    static_assert(_degree_range == 0,
                  "The degree range of a zero-d polynomial "
                  "must be zero");
  }

  Polynomial(const Zero_Tag &) : value(0) {
    static_assert(_degree_range == 0,
                  "The degree range of a zero-d polynomial "
                  "must be zero");
  }

  CoeffT coeff() const noexcept { return value; }

  CoeffT &coeff() noexcept { return value; }

  CoeffT coeff(const Array<int, 0> &exponents) const
      noexcept {
    return value;
  }

  CoeffT &coeff(const Array<int, 0> &exponents) noexcept {
    return value;
  }

  template <typename, int, int>
  friend class Polynomial;

 private:
  static int get_coeff_idx() noexcept { return 0; }

  static int get_coeff_idx_helper(int exp_left) noexcept {
    return exp_left == 0 ? 0 : exp_left > 0
                                   ? DEGREE_RANGE_OVERFLOW
                                   : DEGREE_RANGE_UNDERFLOW;
  }

  static constexpr const int dim = 0;

  static constexpr const int num_coeffs = 1;
  CoeffT value;

  enum {
    DEGREE_RANGE_UNDERFLOW = -1,
    DEGREE_RANGE_OVERFLOW = -2,
  };
};
}

#endif  //_POLYNOMIAL_HPP_