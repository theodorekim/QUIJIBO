#include "RATIONAL_4D.h"
#include "TIMER.h"

RATIONAL_4D::RATIONAL_4D(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom) :
  _top(top),
  _bottom(bottom)
{
}

QUATERNION RATIONAL_4D::evaluate(const QUATERNION& iterate)
{
  QUATERNION topEval = _top.evaluateScaledPowerFactored(iterate);
  QUATERNION bottomEval;

  QUATERNION result;  
  if (_bottom.totalRoots() > 0)
  {
    bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    result = (topEval / bottomEval);
  }
  else
    result = topEval;

  return result;
}

QUATERNION RATIONAL_4D::evaluateLagrange(const QUATERNION& iterate, const bool debug)
{
#if 1
  QUATERNION result(1,0,0,0);
  const int topRoots = _top.totalRoots();
  const int bottomRoots = _bottom.totalRoots();

  if (bottomRoots == 0)
    return _top.evaluateScaledPowerFactored(iterate);

  const int sameRoots = (topRoots > bottomRoots) ? bottomRoots : topRoots;
  if (debug)
    cout << " result 0: " << result << endl;

  // compute the Lagrange form for top and bottom
  for (int x = 0; x < sameRoots; x++)
  {
    QUATERNION topEval = _top.evaluateRoot(iterate,x, debug);
    if (debug)
      result.debugMultiply(topEval);
    else
      result *= topEval;
    if (debug) cout << " result " << x << " top: " << result << endl;
    result = result / _bottom.evaluateRoot(iterate,x, debug);
    if (debug) cout << " result " << x << " bottom: " << result << endl;
  }

  // compute the polynomial for the leftovers
  if (topRoots > bottomRoots)
  {
    for (int x = sameRoots; x < topRoots; x++)
    {
      result *= _top.evaluateRoot(iterate, x, debug);
      if (debug) cout << " result " << x << " top: " << result << endl;
    }
  }
  else
  {
    for (int x = sameRoots; x < bottomRoots; x++)
    {
      result = result / _bottom.evaluateRoot(iterate, x, debug);
      if (debug) cout << " result " << x << " bottom: " << result << endl;
    }
  }

  return result;
#else
  QUATERNION result;  
  QUATERNION topEval = _top.evaluateScaledPowerFactored(iterate);
  QUATERNION bottomEval;

  if (_bottom.totalRoots() > 0)
  {
    bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    result = (topEval / bottomEval);
  }
  else
    result = topEval;

  return result;
#endif
}
