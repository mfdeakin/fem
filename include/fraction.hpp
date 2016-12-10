
#ifndef _FRACTION_HPP_
#define _FRACTION_HPP_

namespace Numerical {

class Fraction {
 public:
  Fraction(int num, int den) {
    this->num = num;
    this->den = den;
  }

	bool isValid() const {
		return den != 0;
	}

  Fraction operator+(int val) const {
    return Fraction(val * this->den + this->num, this->den);
  }

  Fraction operator-(int val) const {
    return Fraction(-val * this->den + this->num,
                    this->den);
  }

  Fraction operator*(int val) const {
    return Fraction(val * this->num, this->den);
  }

  Fraction operator/(int val) const {
    return Fraction(this->num, val * this->den);
  }

  Fraction operator+=(int val) {
    this->num += val * this->den;
    return *this;
  }

  Fraction operator-=(int val) {
    this->num -= val * this->den;
    return *this;
  }

  Fraction operator*=(int val) {
    this->num *= val;
    return *this;
  }

  Fraction operator/=(int val) {
    this->den *= val;
    return *this;
  }

 private:
  int num, den;
};
}

#endif  //_FRACTION_HPP_
