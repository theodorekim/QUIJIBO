#include "QUATERNION_T.h"

#if 0
//////////////////////////////////////////////////////////////////////
// function specialization for mpreal
//////////////////////////////////////////////////////////////////////
template<>
QUATERNION_T<mpreal> QUATERNION_T<mpreal>::pow(const mpreal& exponent) const
{
  const mpreal partial = _x * _x + _y * _y + _z * _z;
  const mpreal qMagnitude = sqrt(partial + _w * _w);
  const mpreal vMagnitude = sqrt(partial);
  const mpreal vMagnitudeInv = (vMagnitude > 0) ? 1.0 / vMagnitude : 0;

  const mpreal scale = exponent * acos(_w / qMagnitude) * vMagnitudeInv;

  const mpreal magnitude = scale * vMagnitude;
  const mpreal magnitudeInv = (magnitude > 0) ? 1.0 / magnitude : 0;

  // namespace problem (mpfr vs. std) is what creates the need for
  // a functon specializaton here
  const mpreal exps = mpfr::exp(exponent * mpfr::log(qMagnitude));

  const mpreal scale2 = scale * exps * magnitudeInv * sin(magnitude);
  return QUATERNION_T<mpreal>(exps * cos(magnitude),
                    scale2 * _x,
                    scale2 * _y,
                    scale2 * _z);
}
#endif
