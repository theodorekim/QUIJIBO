#ifndef TOUPEE_SETTINGS_H
#define TOUPEE_SETTINGS_H

////////////////////////////////////////////////////////////////////////////////
// What precision do you want to use for the non-linear solver?
// (uncomment one)
////////////////////////////////////////////////////////////////////////////////
//#define USING_SINGLE_NLS 1
//#define USING_DOUBLE_NLS 1
#define USING_EXTENDED_NLS 1

////////////////////////////////////////////////////////////////////////////////
// What precision do you want to use for the linear (Eigen) solves?
// (uncomment one)
////////////////////////////////////////////////////////////////////////////////
//#define USING_SINGLE_EIGEN 1
#define USING_DOUBLE_EIGEN 1
//#define USING_EXTENDED_EIGEN 1

////////////////////////////////////////////////////////////////////////////////
// You should not need to touch anything below here, at least to get started.
////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXf;
using Eigen::VectorXf;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXe;
typedef Eigen::Matrix<long double, Eigen::Dynamic, 1> VectorXe;

#include <cassert>
#include <float.h>

#if USING_SINGLE_EIGEN
typedef Eigen::MatrixXf MATRIX;
typedef Eigen::VectorXf VECTOR;
#endif
#if USING_DOUBLE_EIGEN
typedef Eigen::MatrixXd MATRIX;
typedef Eigen::VectorXd VECTOR;
#endif
#if USING_EXTENDED_EIGEN
typedef MatrixXe MATRIX;
typedef VectorXe VECTOR;
#endif

#ifndef Real
#if USING_SINGLE_NLS
#define Real float
#endif
#if USING_DOUBLE_NLS
#define Real double
#endif
#if USING_EXTENDED_NLS
#define Real long double
#endif

// the other thresholds are too small to detect things in time
#define REAL_MIN FLT_MIN

// Note that quad precision usually falls back to software, so might as well use MPFR
//#define Real __float128
//#define QUAD_PRECISION 1
#endif

// If you really want to use quad precision, some of the comment math functions
// do not have the common names.
/*
#ifdef QUAD_PRECISION
#include <quadmath.h>
#include <iostream>

inline std::ostream& operator<<(std::ostream &out, const __float128& q)
{
  out << (double)q; 
  return out;
}

inline __float128 fabs(const __float128& q)  { return fabsq(q); }
inline __float128 log(const __float128& q)   { return logq(q); }
inline __float128 sqrt(const __float128& q)  { return sqrtq(q); }
inline __float128 floor(const __float128& q) { return floorq(q); }
inline __float128 ceil(const __float128& q)  { return ceilq(q); }
inline __float128 acos(const __float128& q)  { return acosq(q); }
inline __float128 sin(const __float128& q)   { return sinq(q); }
inline __float128 cos(const __float128& q)   { return cosq(q); }
inline __float128 exp(const __float128& q)   { return expq(q); }
inline __float128 atan2(const __float128& q, const __float128& p) { return atan2q(q,p); }
inline __float128 pow(const __float128& q, const __float128& p)   { return powq(q,p); }
#endif
*/

#ifndef NDEBUG
#define NDEBUG
#endif

#include <GL/glut.h>

#endif
