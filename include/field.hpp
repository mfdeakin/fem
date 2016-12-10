
#ifndef _FIELD_HPP_
#define _FIELD_HPP_

namespace Numerical {

class FieldType {
 public:
  constexpr static const FieldType &addId();
  constexpr static const FieldType &multiplyId();
  FieldType operator+(const FieldType) const;
	FieldType operator-() const;
  FieldType operator*(const FieldType) const;
};
}

#endif  // _FIELD_HPP_
