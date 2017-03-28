/*
QUIJIBO: Source code for the paper Symposium on Geometry Processing
         2015 paper "Quaternion Julia Set Shape Optimization"
Copyright (C) 2015  Theodore Kim

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _QUATERNION_T_H
#define _QUATERNION_T_H

#include <SETTINGS.h>
#include <VEC3.h>
#include <MATRIX3.h>
#include <iostream>
#include <QUATERNION.h>
//#include "mpreal.h"

using namespace std;
//using mpfr::mpreal;

//////////////////////////////////////////////////////////////////////
// _w and _z are the real component
//////////////////////////////////////////////////////////////////////
template<class Float>
class QUATERNION_T {

public:
  QUATERNION_T();
  QUATERNION_T(Float w, Float x, Float y, Float z);
  QUATERNION_T(Float w, Float x);
  QUATERNION_T(const VEC3F& vector);
  QUATERNION_T(const QUATERNION_T& q);
  QUATERNION_T(const QUATERNION& q);

  // overload operators
  QUATERNION_T& operator=(const VEC3F& vec);
  inline QUATERNION_T& operator=(const QUATERNION_T q)   { _w = q.w(); _x = q.x(); _y = q.y(); _z = q.z(); return *this;};
  inline QUATERNION_T& operator*=(const Float r)        { _w *= r; _x *= r; _y *= r; _z *= r;  return *this; };
  inline QUATERNION_T& operator*=(const QUATERNION_T& q) { 
    const Float& x = this->_y * q._z - this->_z * q._y + q._w * this->_x + this->_w * q._x;
    const Float& y = this->_z * q._x - this->_x * q._z + q._w * this->_y + this->_w * q._y;
    const Float& z = this->_x * q._y - this->_y * q._x + q._w * this->_z + this->_w * q._z;
    const Float& w = this->_w * q._w - this->_x * q._x - q._y * this->_y - this->_z * q._z;
    _x = x; _y = y; _z = z; _w = w;
    return *this;
  };
  inline QUATERNION_T& operator-=(const QUATERNION_T& q) { _w -= q.w(); _x -= q.x(); _y -= q.y(); _z -= q.z(); return *this; };
  inline QUATERNION_T& operator+=(const QUATERNION_T& q) { _w += q.w(); _x += q.x(); _y += q.y(); _z += q.z(); return *this; };
  inline Float& operator[](const int x) { 
    if (x == 0) return _w;
    if (x == 1) return _x;
    if (x == 2) return _y;
    return _z;
  };
  inline Float operator[](const int x) const { 
    if (x == 0) return _w;
    if (x == 1) return _x;
    if (x == 2) return _y;
    return _z;
  };

  // compute the conjugate
  QUATERNION_T conjugate() const;

  // compute the inverse
  QUATERNION_T inverse() const;

  // populate a VECTOR of size 4
  VECTOR toVector();
  
  // normalize
  void normalize();

  // negative the imaginary component, effectively retrieving the transpose
  // of the desired rotation
  void negateIm();

  // accessors
  const Float& x() const { return _x; };
  const Float& y() const { return _y; };
  const Float& z() const { return _z; };
  const Float& w() const { return _w; };
  
  Float& x() { return _x; };
  Float& y() { return _y; };
  Float& z() { return _z; };
  Float& w() { return _w; };

  // 2 norm of components
  inline Float magnitude() const { return sqrt(_w * _w + _x * _x + _y * _y + _z * _z); };

  // in-place equality
  void equals(const QUATERNION_T& q);

  // check if any components are nan
  inline bool anyNans() { return _x != _x || _y != _y || _z != _z || _w != _w; };

  // do a julia set iteration
  void juliaIteration(const QUATERNION_T& c);

  void write(FILE* file) const;
  void read(FILE* file);

  inline void multiplyAdd(const QUATERNION_T& point, const QUATERNION_T& add)
  {
    const Float x = _x;
    const Float y = _y;
    const Float z = _z;
    const Float w = _w;
    _x = y * point._z - z * point._y + point._w * x + w * point._x + add._x;
    _y = z * point._x - x * point._z + point._w * y + w * point._y + add._y;
    _z = x * point._y - y * point._x + point._w * z + w * point._z + add._z;
    _w = w * point._w - x * point._x - point._y * y - z * point._z + add._w;
  };

  static bool wCompare(const QUATERNION_T& i, const QUATERNION_T& j) { return i._w < j._w; };
  static bool xCompare(const QUATERNION_T& i, const QUATERNION_T& j) { return i._x < j._x; };
  static bool yCompare(const QUATERNION_T& i, const QUATERNION_T& j) { return i._y < j._y; };
  static bool zCompare(const QUATERNION_T& i, const QUATERNION_T& j) { return i._z < j._z; };

  // take the exponential
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION_T exp() const;

  // take the log
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION_T log() const;

  // take the power
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION_T pow(const Float& exponent) const;

  // dot against another quaternion
  inline Float dot(const QUATERNION_T& rhs) const { return _w * rhs._w + _x * rhs._x + _y * rhs._y + _z * rhs._z; };

private:
  /*
  // _w and _entries[0] are the real component
  union {
     struct { Float _w,_x,_y,_z; };
     Float _entries[4];
     //__m128 _v;
  };
  */
  Float _w, _x, _y, _z;
};

template<class Float>
QUATERNION_T<Float> operator/(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right);

template<class Float>
QUATERNION_T<Float> operator/(const QUATERNION_T<Float>& left, const Float& right);

template<class Float>
Float operator^(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right);

template<class Float>
ostream &operator<<(ostream &out, const QUATERNION_T<Float>& q);

template<class Float>
inline QUATERNION_T<Float> operator*(const QUATERNION_T<Float>& left, const Float& right) {
  return QUATERNION_T<Float>(left.w() * right, left.x() * right,
                             left.y() * right, left.z() * right);
};

template<class Float>
inline QUATERNION_T<Float> operator*(const Float& left, const QUATERNION_T<Float>& right) {
  return QUATERNION_T<Float>(right.w() * left, right.x() * left,
                             right.y() * left, right.z() * left);
};

template<class Float>
inline QUATERNION_T<Float> operator-(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right) {
  return QUATERNION_T<Float>(left.w() - right.w(), left.x() - right.x(),
                             left.y() - right.y(), left.z() - right.z());
};

template<class Float>
inline QUATERNION_T<Float> operator+(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right) {
  return QUATERNION_T<Float>(left.w() + right.w(), left.x() + right.x(),
                             left.y() + right.y(), left.z() + right.z());
};

template<class Float>
QUATERNION_T<Float> operator*(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right);


template<class Float>
QUATERNION_T<Float>::QUATERNION_T() :
  _w(0), _x(0), _y(0), _z(0)
{
}

template<class Float>
QUATERNION_T<Float>::QUATERNION_T(Float w, Float x, Float y, Float z) :
  _w(w), _x(x), _y(y), _z(z)
{
}

template<class Float>
QUATERNION_T<Float>::QUATERNION_T(Float w, Float x) :
  _w(w), _x(x), _y(0), _z(0)
{
}

template<class Float>
QUATERNION_T<Float>::QUATERNION_T(const VEC3F& vector) :
  _w(0), _x(vector[0]), _y(vector[1]), _z(vector[2])
{
}

template<class Float>
QUATERNION_T<Float>::QUATERNION_T(const QUATERNION_T& q) :
  _w(q._w), _x(q._x), _y(q._y), _z(q._z)
{
}

template<class Float>
QUATERNION_T<Float>::QUATERNION_T(const QUATERNION& q) :
  _w(q.w()), _x(q.x()), _y(q.y()), _z(q.z())
{
}

template<class Float>
QUATERNION_T<Float> QUATERNION_T<Float>::conjugate() const
{
  QUATERNION_T final(_w, -_x, -_y, -_z);
  return final;
}

template<class Float>
QUATERNION_T<Float>& QUATERNION_T<Float>::operator=(const VEC3F& vec)
{
  this->_x = vec[0];
  this->_y = vec[1];
  this->_z = vec[2];
  this->_w = 0.0;
  return *this;
}

template<class Float>
QUATERNION_T<Float> operator*(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right)
{
  QUATERNION_T<Float> final;
  final.x() = left.y() * right.z() - left.z() * right.y() + right.w() * left.x() + left.w() * right.x();
  final.y() = left.z() * right.x() - left.x() * right.z() + right.w() * left.y() + left.w() * right.y();
  final.z() = left.x() * right.y() - left.y() * right.x() + right.w() * left.z() + left.w() * right.z();
  final.w() = left.w() * right.w() - left.x() * right.x() - right.y() * left.y() - left.z() * right.z();
  return final;
}

//////////////////////////////////////////////////////////////////////
// compute the inverse
//////////////////////////////////////////////////////////////////////
template<class Float>
QUATERNION_T<Float> QUATERNION_T<Float>::inverse() const
{
  const Float magnitudeSq = (*this)[0] * (*this)[0] + (*this)[1] * (*this)[1] +
                           (*this)[2] * (*this)[2] + (*this)[3] * (*this)[3];
  return QUATERNION_T<Float>(conjugate() / magnitudeSq);
}

template<class Float>
QUATERNION_T<Float> operator/(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right)
{
  // build the reciprocal of the right hand side
  Float magnitudeSq = right[0] * right[0] + right[1] * right[1] +
                      right[2] * right[2] + right[3] * right[3];
  QUATERNION_T<Float> rightInv = right.conjugate() / magnitudeSq;

  return left * rightInv;
}

template<class Float>
Float operator^(const QUATERNION_T<Float>& left, const QUATERNION_T<Float>& right)
{
  return left._x * right._x + left._y * right._y +
         left._z * right._z + left._w * right._w;
}

template<class Float>
ostream& operator<<(ostream &out, const QUATERNION_T<Float>& q)
{
  out << "(" << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << ")";
  return out;
}

template<class Float>
QUATERNION_T<Float> operator/(const QUATERNION_T<Float>& left, const Float& right)
{
  const Float rightInv = 1.0 / right;
  return left * rightInv;
}

template<class Float>
void QUATERNION_T<Float>::normalize()
{
  Float invMagnitude = 1.0 / sqrtf(_x * _x + _y * _y + _z * _z + _w * _w);
  _x *= invMagnitude;
  _y *= invMagnitude;
  _z *= invMagnitude;
  _w *= invMagnitude;
}

//////////////////////////////////////////////////////////////////////
// negate the imaginary component -- effetively transpose the rotation
//////////////////////////////////////////////////////////////////////
template<class Float>
void QUATERNION_T<Float>::negateIm()
{
  _x *= -1.0;
  _y *= -1.0;
  _z *= -1.0;
}

//////////////////////////////////////////////////////////////////////
// populate a VECTOR of size 4
//////////////////////////////////////////////////////////////////////
template<class Float>
VECTOR QUATERNION_T<Float>::toVector()
{
  VECTOR final(4);
  final[0] = _w;
  final[1] = _x;
  final[2] = _y;
  final[3] = _z;
  return final;
}

//////////////////////////////////////////////////////////////////////
// in-place equality
//////////////////////////////////////////////////////////////////////
template<class Float>
void QUATERNION_T<Float>::equals(const QUATERNION_T& q)
{
  _w = q._w;
  _x = q._x;
  _y = q._y;
  _z = q._z;
}

//////////////////////////////////////////////////////////////////////
// do a julia set iteration
//////////////////////////////////////////////////////////////////////
template<class Float>
void QUATERNION_T<Float>::juliaIteration(const QUATERNION_T& c)
{
  const Float copy[] = {_x, _y, _z, _w};
  
  //const Float copy01 = copy[0] * copy[1];
  //const Float copy12 = copy[1] * copy[2];
  //const Float copy02 = copy[0] * copy[2];
  const Float copy03 = 2.0f * copy[0] * copy[3];
  const Float copy31 = 2.0f * copy[3] * copy[1];
  const Float copy32 = 2.0f * copy[3] * copy[2];

  //_x = copy12 - copy12 + copy03 + copy03 + c[0];
  //_y = copy02 - copy02 + copy31 + copy31 + c[1];
  //_z = copy01 - copy01 + copy32 + copy32 + c[2];
  _x = copy03 + c[0];
  _y = copy31 + c[1];
  _z = copy32 + c[2];
  _w = copy[3] * copy[3] - copy[0] * copy[0] - copy[1] * copy[1] - copy[2] * copy[2] + c[3];
}

template<class Float>
void QUATERNION_T<Float>::write(FILE* file) const
{
  if (sizeof(Float) == sizeof(double))
  {
    fwrite((void*)&_w, sizeof(Float), 1, file);
    fwrite((void*)&_x, sizeof(Float), 1, file);
    fwrite((void*)&_y, sizeof(Float), 1, file);
    fwrite((void*)&_z, sizeof(Float), 1, file);
  }
  else
  {
    double entries[] = {_w, _x, _y, _z};
    fwrite((void*)&entries[0], sizeof(double), 1, file);
    fwrite((void*)&entries[1], sizeof(double), 1, file);
    fwrite((void*)&entries[2], sizeof(double), 1, file);
    fwrite((void*)&entries[3], sizeof(double), 1, file);
  }
}

template<class Float>
void QUATERNION_T<Float>::read(FILE* file)
{
  if (sizeof(Float) == sizeof(double))
  {
    fread((void*)&_w, sizeof(Float), 1, file);
    fread((void*)&_x, sizeof(Float), 1, file);
    fread((void*)&_y, sizeof(Float), 1, file);
    fread((void*)&_z, sizeof(Float), 1, file);
  }
  else
  {
    double entries[4];
    fread((void*)&entries[0], sizeof(double), 1, file);
    fread((void*)&entries[1], sizeof(double), 1, file);
    fread((void*)&entries[2], sizeof(double), 1, file);
    fread((void*)&entries[3], sizeof(double), 1, file);

    _w = entries[0];
    _x = entries[1];
    _y = entries[2];
    _z = entries[3];
  }
}

//////////////////////////////////////////////////////////////////////
// take the exponential
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//////////////////////////////////////////////////////////////////////
template<class Float>
QUATERNION_T<Float> QUATERNION_T<Float>::exp() const
{
  const Float& s =  _w;
  const VEC3F v(_x, _y, _z);

  const Float magnitude = v.magnitude();
  const Float exps = std::exp(s);
  const VEC3F vFinal = (v / magnitude) * sin(magnitude);

  return exps * QUATERNION_T(cos(magnitude), vFinal[0], vFinal[1], vFinal[2]);
}

//////////////////////////////////////////////////////////////////////
// take the log 
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//////////////////////////////////////////////////////////////////////
template<class Float>
QUATERNION_T<Float> QUATERNION_T<Float>::log() const
{
  const Float& s =  _w;
  const VEC3F v(_x, _y, _z);

  const Float qMagnitude = magnitude();
  const Float vMagnitude = v.magnitude();

  const VEC3F vNormalized = (vMagnitude > 0) ? v / (Real)vMagnitude : VEC3F(0,0,0);

  const VEC3F vFinal = vNormalized * (Real)acos(s / (Real)qMagnitude);

  return QUATERNION_T(std::log(qMagnitude), vFinal[0], vFinal[1], vFinal[2]);
}

//////////////////////////////////////////////////////////////////////
// take the power
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//
// some care has been taken to optimize this ... 
//////////////////////////////////////////////////////////////////////
template<class Float>
QUATERNION_T<Float> QUATERNION_T<Float>::pow(const Float& exponent) const
{
  const Float partial = _x * _x + _y * _y + _z * _z;
  const Float qMagnitude = sqrt(partial + _w * _w);
  const Float vMagnitude = sqrt(partial);
  const Float vMagnitudeInv = (vMagnitude > 0) ? 1.0 / vMagnitude : 0;

  const Float scale = exponent * acos(_w / qMagnitude) * vMagnitudeInv;

  const Float magnitude = scale * vMagnitude;
  const Float magnitudeInv = (magnitude > 0) ? 1.0 / magnitude : 0;
  const Float exps = std::exp(exponent * std::log(qMagnitude));

  const Float scale2 = scale * exps * magnitudeInv * sin(magnitude);
  return QUATERNION_T(exps * cos(magnitude),
                    scale2 * _x,
                    scale2 * _y,
                    scale2 * _z);
}

#if 0
//////////////////////////////////////////////////////////////////////
// function specialization for mpreal
//////////////////////////////////////////////////////////////////////
template<>
QUATERNION_T<mpreal> QUATERNION_T<mpreal>::pow(const mpreal& exponent) const;
#endif

typedef QUATERNION_T<float> QUATERNIONF;
typedef QUATERNION_T<double> QUATERNIOND;
typedef QUATERNION_T<long double> QUATERNIONE;
//typedef QUATERNION_T<mpreal> QUATERNIONMP;

#endif
