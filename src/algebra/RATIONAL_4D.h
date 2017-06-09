#ifndef RATIONAL_4D_H
#define RATIONAL_4D_H

#include "POLYNOMIAL_4D.h"

using namespace std;

class RATIONAL_4D {
public:
  RATIONAL_4D(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom);
  QUATERNION evaluate(const QUATERNION& iterate);
  QUATERNION evaluateLagrange(const QUATERNION& iterate, bool debug = false);

private:
  const POLYNOMIAL_4D& _top;
  const POLYNOMIAL_4D& _bottom;
};
#endif
