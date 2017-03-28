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
#ifndef POLYNOMIAL_4D_T_H
#define POLYNOMIAL_4D_T_H

#include <cmath>
#include <vector>
#include "VEC3.h"
#include "QUATERNION_T.h"
#include "TIMER.h"

#include <POLYNOMIAL_4D.h>

using namespace std;

template<class Float>
class POLYNOMIAL_4D_T {
public:
  POLYNOMIAL_4D_T();
  POLYNOMIAL_4D_T(const POLYNOMIAL_4D& polynomial);
  ~POLYNOMIAL_4D_T();

  // compute coefficients from the roots
  POLYNOMIAL_4D_T(const vector<QUATERNION_T<Float> >& roots);
  POLYNOMIAL_4D_T(const vector<QUATERNION_T<Float> >& roots, const vector<Float> powers);

  // set the coefficients directly - powers increase with index, i.e.
  // coeffs[0] -> x^0
  // coeffs[1] -> x^1
  POLYNOMIAL_4D_T(const vector<float>& coeffs);

  QUATERNION_T<Float> evaluate(const QUATERNION_T<Float>& point) const;
  QUATERNION_T<Float> evaluateDerivative(const QUATERNION_T<Float>& point) const;
  QUATERNION_T<Float> evaluateFactoredDerivative(const QUATERNION_T<Float>& point) const;
  QUATERNION_T<Float> evaluateSecondDerivative(const QUATERNION_T<Float>& point) const;
  void evaluateMultiple(const QUATERNION_T<Float>& point, QUATERNION_T<Float>& poly, QUATERNION_T<Float>& deriv) const;
  void evaluateMultiple(const QUATERNION_T<Float>& point, QUATERNION_T<Float>& poly, QUATERNION_T<Float>& deriv, QUATERNION_T<Float>& secondDeriv) const;

  // use the brute force nested formulation
  QUATERNION_T<Float> evaluateFactored(const QUATERNION_T<Float>& point) const;
  QUATERNION_T<Float> evaluatePowerFactored(const QUATERNION_T<Float>& point) const;
  QUATERNION_T<Float> evaluateScaledPowerFactored(const QUATERNION_T<Float>& point) const;
  QUATERNION_T<Float> evaluateFactoredDouble(const QUATERNION_T<Float>& point) const;
  QUATERNION_T<Float> evaluateFactoredPositive(const QUATERNION_T<Float>& point) const;
  void computeNestedCoeffs();

  // use the brute force nested formulation, but cache the multiplies
  QUATERNION_T<Float> evaluatePowerFactored(const QUATERNION_T<Float>& point, vector<QUATERNION_T<Float> >& forward, vector<QUATERNION_T<Float> >& backward) const;
  QUATERNION_T<Float> evaluateScaledPowerFactored(const QUATERNION_T<Float>& point, vector<QUATERNION_T<Float> >& forward, vector<QUATERNION_T<Float> >& backward) const;

  // add a new root
  void addRoot(const QUATERNION_T<Float>& newRoot);
  void addRoot(const QUATERNION_T<Float>& newRoot, const Float& power);
  void addFrontRoot(const QUATERNION_T<Float>& newRoot);
  
  // modify an existing root
  void modifyRoot(const int whichRoot, const QUATERNION_T<Float>& newRoot);

  // accessors
  const int totalRoots() const { return _totalRoots; };
  const vector<QUATERNION_T<Float> >& roots() const { return _roots; };
  const vector<QUATERNION_T<Float> >& coeffs() const { return _coeffs; };
  const vector<QUATERNION_T<Float> >& derivs() const { return _derivs; };
  const vector<QUATERNION_T<Float> >& secondDerivs() const { return _secondDerivs; };
  const vector<Float>& powers() const { return _rootPowers; };
  vector<QUATERNION_T<Float> >& rootsMutable() { return _roots; };
  vector<Float>& powersMutable() { return _rootPowers; };
  Float& powerScalar() { return _powerScalar; };
  const Float powerScalar() const { return _powerScalar; };

  // run a unit test for a known rational function
  static void rationalTest();

  // evaluate rational function
  static void evaluateRational(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime, QUATERNION_T<Float>& pPrimePrime);
  static void evaluateRational(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime);
  
  static void evaluateFactoredRational(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime);

  // for debugging purposes
  static void evaluateFactoredQuadratic(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime);

  // resize the polynomial to a given number of roots
  void resizeAndWipe(int totalRoots);

  // file IO
  void write(FILE* file) const;
  void read(FILE* file);

  // get the condition number of the polynomial
  Float conditionNumber();

  // sum of all the root coefficients
  Float rootSum() const;

  // overloaded operators
  POLYNOMIAL_4D_T& operator*=(const Float& alpha);
  POLYNOMIAL_4D_T& operator-=(const VEC3F& v);
  POLYNOMIAL_4D_T& operator+=(const VEC3F& v);

  void translateExceptFirst(const VEC3F& v);

  // change the power of a root
  void changePower(const int& whichRoot, const Float& newPower) {
    assert(whichRoot < _totalRoots);
    _rootPowers[whichRoot] = newPower;
  };

  // take the derivative with respect to a root
  QUATERNION_T<Float> powerDerivative(const QUATERNION_T<Float>& point, const int whichRoot) const;
  QUATERNION_T<Float> inversePowerDerivative(const QUATERNION_T<Float>& point, const int whichRoot) const;

  // compute the gradient with respect to each root
  VECTOR powerGradient(const QUATERNION_T<Float>& point);

  // test out taking the derivative of a power
  static void testSingleDerivative();
  static void testPolynomialDerivative();
  static void testBulkDerivative();

  // stomp everything
  void clear();

private:
  int _totalRoots;
  vector<QUATERNION_T<Float> > _coeffs;
  vector<QUATERNION_T<Float> > _derivs;
  vector<QUATERNION_T<Float> > _secondDerivs;
  vector<QUATERNION_T<Float> > _roots;

  vector<Float> _rootPowers;

  // the nested linear factorization
  vector<QUATERNION_T<Float> > _ws;

  // a uniform scalar to multiply all the powers by
  Float _powerScalar;

  // compute the polynomial coefficients
  void computeCoeffs();
 
  // compute the derivative coefficients
  void computeDerivativeCoeffs();

  // compute the polynomial coefficients
  void computeCoeffsFast();

  // support function for coefficients
  vector<QUATERNION_T<Float> > choose(const vector<QUATERNION_T<Float> > m, const int n);
};

template <class Float>
ostream& operator<<(ostream &out, const POLYNOMIAL_4D_T<Float>& poly);

template <class Float>
POLYNOMIAL_4D_T<Float>::POLYNOMIAL_4D_T(const POLYNOMIAL_4D& polynomial)
{
  _totalRoots = polynomial.totalRoots();

  for (unsigned int x = 0; x < polynomial.coeffs().size(); x++)
    _coeffs.push_back(QUATERNION_T<Float>(polynomial.coeffs()[x]));
  
  for (unsigned int x = 0; x < polynomial.derivs().size(); x++)
    _derivs.push_back(QUATERNION_T<Float>(polynomial.derivs()[x]));
  
  for (unsigned int x = 0; x < polynomial.secondDerivs().size(); x++)
    _secondDerivs.push_back(QUATERNION_T<Float>(polynomial.secondDerivs()[x]));
  
  for (unsigned int x = 0; x < polynomial.roots().size(); x++)
    _roots.push_back(QUATERNION_T<Float>(polynomial.roots()[x]));
  
  for (unsigned int x = 0; x < polynomial.powers().size(); x++)
    _rootPowers.push_back(polynomial.powers()[x]);

  _powerScalar = polynomial.powerScalar();  
}

template <class Float>
POLYNOMIAL_4D_T<Float>::POLYNOMIAL_4D_T(const vector<QUATERNION_T<Float> >& roots)
{
  _powerScalar = 1.0;

  _roots = roots;
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  if (roots.size() == 0)
    return;

  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < _totalRoots; x++)
    _rootPowers[x] = 1.0;

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

template <class Float>
POLYNOMIAL_4D_T<Float>::POLYNOMIAL_4D_T(const vector<QUATERNION_T<Float> >& roots, const vector<Float> powers)
{
  assert(roots.size() == powers.size());
  _powerScalar = 1.0;

  _roots = roots;
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  if (roots.size() == 0)
    return;

  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < _totalRoots; x++)
    _rootPowers[x] = powers[x];

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

///////////////////////////////////////////////////////////////////////
// set the coefficients directly - powers increase with index, i.e.
// coeffs[0] -> x^0
// coeffs[1] -> x^1
///////////////////////////////////////////////////////////////////////
template <class Float>
POLYNOMIAL_4D_T<Float>::POLYNOMIAL_4D_T(const vector<float>& coeffs)
{
  _coeffs.resize(coeffs.size());
  _powerScalar = 1.0;

  for (unsigned int x = 0; x < coeffs.size(); x++)
    _coeffs[x] = QUATERNION_T<Float>(coeffs[x], 0);
  _totalRoots = _coeffs.size() - 1;

  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < _totalRoots; x++)
    _rootPowers[x] = 1.0;

  computeDerivativeCoeffs();
}

template <class Float>
POLYNOMIAL_4D_T<Float>::POLYNOMIAL_4D_T()
{
  _totalRoots = -1;
  _powerScalar = 1.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
template <class Float>
POLYNOMIAL_4D_T<Float>::~POLYNOMIAL_4D_T()
{
}

///////////////////////////////////////////////////////////////////////
// m choose n
///////////////////////////////////////////////////////////////////////
template <class Float>
vector<QUATERNION_T<Float> > POLYNOMIAL_4D_T<Float>::choose(const vector<QUATERNION_T<Float> > m, const int n)
{
  if (n == 1)
    return m;

  vector<QUATERNION_T<Float> > final;
  for (unsigned int x = 0; x < m.size(); x++)
  {
    vector<QUATERNION_T<Float> > subset;
    for (unsigned int y = x + 1; y < m.size(); y++)
      subset.push_back(m[y]);

    vector<QUATERNION_T<Float> > chooseOneLess = choose(subset, n - 1);

    for (unsigned int y = 0; y < chooseOneLess.size(); y++)
      final.push_back(chooseOneLess[y] * m[x]);
  }
  return final;
}

//////////////////////////////////////////////////////////////////////
// resize the polynomial to a given number of roots
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::resizeAndWipe(int totalRoots)
{
  _roots.clear();
  for (int x = 0; x < totalRoots; x++)
    _roots.push_back(QUATERNION_T<Float>(0,0,0,0));

  _totalRoots = totalRoots;
  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < totalRoots; x++)
    _rootPowers[x] = 1.0;

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// generically compute new root coefficients
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::computeCoeffs()
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " SHOULDN'T BE CALLING THIS " << endl;
  assert(_totalRoots > 0);

  vector<QUATERNION_T<Float> > chosen;

  for (int x = 0; x < _totalRoots; x++)
    _coeffs[x] = 0;

  for (int x = 0; x < _totalRoots; x++)
  {
    chosen = choose(_roots, _totalRoots - x);
    for (unsigned int y = 0; y < chosen.size(); y++)
      _coeffs[x] += chosen[y];

    if (x % 2 != _totalRoots % 2)
      _coeffs[x] *= -1;
  }
  //_coeffs[_totalRoots] = QUATERNION_T<Float>(1.0,0);
  _coeffs[_totalRoots] = QUATERNION_T<Float>(0,0,0,1.0);
}

//////////////////////////////////////////////////////////////////////
// compute the polynomial coefficients
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::computeCoeffsFast()
{
  //TIMER functionTimer(__FUNCTION__);
  vector<QUATERNION_T<Float> > coeffs;
  vector<QUATERNION_T<Float> > coeffsOld;

  coeffs.push_back(QUATERNION_T<Float>(1.0, 0.0));
  coeffs.push_back(-1.0 * _roots[0]);

  for (int x = 1; x < _totalRoots; x++)
  {
    coeffsOld = coeffs;
    coeffs.clear();

    QUATERNION_T<Float> alpha = _roots[x];
    coeffs.push_back(coeffsOld[0]);

    for (unsigned int y = 1; y < coeffsOld.size(); y++)
      coeffs.push_back(coeffsOld[y] - alpha * coeffsOld[y - 1]);

    coeffs.push_back(-1.0 * alpha * coeffsOld.back());
  }

  _coeffs.clear();
  for (unsigned int x = 0; x < coeffs.size(); x++)
    _coeffs.push_back(coeffs[coeffs.size() - 1 - x]);
}

//////////////////////////////////////////////////////////////////////
// compute the derivative coefficients
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::computeDerivativeCoeffs()
{
  _derivs.resize(_totalRoots);
  _secondDerivs.resize(_totalRoots - 1);

  for (int x = 0; x < _totalRoots; x++)
    _derivs[x] = (x + 1.0) * _coeffs[x + 1];
  for (int x = 0; x < _totalRoots - 1; x++)
    _secondDerivs[x] = (x + 1.0) * _derivs[x + 1];
}

//////////////////////////////////////////////////////////////////////
// evaluate the polynomial at a point
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluate(const QUATERNION_T<Float>& point) const
{
  assert(_totalRoots > 0);

  int roots = totalRoots();
  QUATERNION_T<Float> final(_coeffs[roots]);
  for (int x = roots - 1; x >= 0; x--)
    //g = g * point + topCoeffs[x];
    final.multiplyAdd(point, _coeffs[x]);

  return final;
}

//////////////////////////////////////////////////////////////////////
// evaluate the derivative at a point
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateDerivative(const QUATERNION_T<Float>& point) const
{
  assert(_totalRoots > 0);

  int roots = totalRoots();
  QUATERNION_T<Float> final(_derivs[roots - 1]);
  for (int x = roots - 2; x >= 0; x--)
    //gPrime = gPrime * point + topDerivs[x];
    final.multiplyAdd(point, _derivs[x]);

  return final;
}

//////////////////////////////////////////////////////////////////////
// evaluate the second derivative at a point
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateSecondDerivative(const QUATERNION_T<Float>& point) const
{
  assert(_totalRoots > 0);

  // compute all the needed powers
  QUATERNION_T<Float> powers[_totalRoots];
  powers[0] = point;
  for (int x = 1; x < _totalRoots; x++)
    powers[x] = powers[x - 1] * point;

  // compute the second derivative
  QUATERNION_T<Float> final = _secondDerivs[0];
  for (int x = 1; x < _totalRoots - 1; x++)
    final += _secondDerivs[x] * powers[x-1];

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::evaluateMultiple(const QUATERNION_T<Float>& point, QUATERNION_T<Float>& poly, QUATERNION_T<Float>& deriv) const
{
  assert(_totalRoots > 0);

  // compute all the needed powers
  QUATERNION_T<Float> powers[_totalRoots];
  powers[0] = point;
  for (int x = 1; x < _totalRoots; x++)
    powers[x] = powers[x - 1] * point;

  // compute the polynomial
  poly = _coeffs[0];
  for (int x = 1; x < _totalRoots + 1; x++)
    poly += _coeffs[x] * powers[x - 1];
  
  // compute the derivative
  deriv = _derivs[0];
  for (int x = 1; x < _totalRoots; x++)
    deriv += _derivs[x] * powers[x-1];
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::evaluateMultiple(const QUATERNION_T<Float>& point, QUATERNION_T<Float>& poly, QUATERNION_T<Float>& deriv, QUATERNION_T<Float>& secondDeriv) const
{
  assert(_totalRoots > 0);

  // compute all the needed powers
  QUATERNION_T<Float> powers[_totalRoots];
  powers[0] = point;
  for (int x = 1; x < _totalRoots; x++)
    powers[x] = powers[x - 1] * point;

  // compute the polynomial
  poly = _coeffs[0];
  for (int x = 1; x < _totalRoots + 1; x++)
    poly += _coeffs[x] * powers[x - 1];
  
  // compute the derivative
  deriv = _derivs[0];
  for (int x = 1; x < _totalRoots; x++)
    deriv += _derivs[x] * powers[x-1];
  
  // compute the second derivative
  secondDeriv = _secondDerivs[0];
  for (int x = 1; x < _totalRoots - 1; x++)
    secondDeriv += _secondDerivs[x] * powers[x-1];
}

//////////////////////////////////////////////////////////////////////
// evaluate rational function
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::evaluateRational(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime, QUATERNION_T<Float>& pPrimePrime)
{
  QUATERNION_T<Float> g, gPrime, gPrimePrime;
  QUATERNION_T<Float> h, hPrime, hPrimePrime;

  const int topRoots = top.totalRoots();
  const int bottomRoots = bottom.totalRoots();
  const int maxRoots = (topRoots > bottomRoots) ? topRoots : bottomRoots;

  // compute all the needed powers
  QUATERNION_T<Float> powers[maxRoots];
  powers[0] = point;
  for (int x = 1; x < maxRoots; x++)
    powers[x] = powers[x - 1] * point;
  
  // compute the g polynomial
  const vector<QUATERNION_T<Float> >& topCoeffs = top.coeffs();
  const vector<QUATERNION_T<Float> >& topDerivs = top.derivs();
  const vector<QUATERNION_T<Float> >& topSecond = top.secondDerivs();
  g = topCoeffs[0];
  for (int x = 1; x < topRoots + 1; x++)
    g += topCoeffs[x] * powers[x - 1];
  gPrime = topDerivs[0];
  for (int x = 1; x < topRoots; x++)
    gPrime += topDerivs[x] * powers[x-1];
  gPrimePrime = topSecond[0];
  for (int x = 1; x < topRoots- 1; x++)
    gPrimePrime += topSecond[x] * powers[x-1];

  // compute the h polynomial
  const vector<QUATERNION_T<Float> >& bottomCoeffs = bottom.coeffs();
  const vector<QUATERNION_T<Float> >& bottomDerivs = bottom.derivs();
  const vector<QUATERNION_T<Float> >& bottomSecond = bottom.secondDerivs();
  h = bottomCoeffs[0];
  for (int x = 1; x < bottomRoots + 1; x++)
    h += bottomCoeffs[x] * powers[x - 1];
  hPrime = bottomDerivs[0];
  for (int x = 1; x < bottomRoots; x++)
    hPrime += bottomDerivs[x] * powers[x-1];
  hPrimePrime = bottomSecond[0];
  for (int x = 1; x < bottomRoots- 1; x++)
    hPrimePrime += bottomSecond[x] * powers[x-1];
 
  // compute the rational derivatives 
  const QUATERNION_T<Float> numerator = gPrime * h - g * hPrime;
  const QUATERNION_T<Float> hSq = h * h;

  const QUATERNION_T<Float> inverse = h.inverse();
  const QUATERNION_T<Float> inverseSq = inverse * inverse;
  const QUATERNION_T<Float> inverseCubed = inverseSq * inverse;

  p = g * inverse;
  pPrime = numerator * inverseSq;

  //const QUATERNION_T<Float> hCubed = hSq * h;
  const QUATERNION_T<Float> first = gPrimePrime * hSq;
  const QUATERNION_T<Float> second = -2.0 * gPrime * hPrime * h;
  const QUATERNION_T<Float> third = g * (2.0 * hPrime * hPrime - h * hPrimePrime);

  pPrimePrime = (first + second + third) * inverseCubed;
}

//////////////////////////////////////////////////////////////////////
// evaluate rational function
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::evaluateRational(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime)
{
  //const int topRoots = top.totalRoots();
  //const int bottomRoots = bottom.totalRoots();
  //const int maxRoots = (topRoots > bottomRoots) ? topRoots : bottomRoots;

  /*
  // compute all the needed powers
  QUATERNION_T<Float> powers[maxRoots];
  powers[0] = point;
  for (int x = 1; x < maxRoots; x++)
    powers[x] = powers[x - 1] * point;
    */
  
  // compute the g polynomial
  //const vector<QUATERNION_T<Float> >& topCoeffs = top.coeffs();
  //const vector<QUATERNION_T<Float> >& topDerivs = top.derivs();
  //const vector<QUATERNION_T<Float> >& topSecond = top.secondDerivs();

  /*
  //QUATERNION_T<Float> g(topCoeffs[0]);
  //for (int x = 1; x < topRoots + 1; x++)
  //  g += topCoeffs[x] * powers[x - 1];
  QUATERNION_T<Float> g(topCoeffs[topRoots]);
  for (int x = topRoots - 1; x >= 0; x--)
    //g = g * point + topCoeffs[x];
    g.multiplyAdd(point, topCoeffs[x]);
  */
  QUATERNION_T<Float> g = top.evaluate(point);

  /*
  //QUATERNION_T<Float> gPrime(topDerivs[0]);
  //for (int x = 1; x < topRoots; x++)
  //  gPrime += topDerivs[x] * powers[x-1];
  QUATERNION_T<Float> gPrime(topDerivs[topRoots - 1]);
  for (int x = topRoots - 2; x >= 0; x--)
    //gPrime = gPrime * point + topDerivs[x];
    gPrime.multiplyAdd(point, topDerivs[x]);
  */
  QUATERNION_T<Float> gPrime = top.evaluateDerivative(point);

  // compute the h polynomial
  //const vector<QUATERNION_T<Float> >& bottomCoeffs = bottom.coeffs();
  //const vector<QUATERNION_T<Float> >& bottomDerivs = bottom.derivs();
  //const vector<QUATERNION_T<Float> >& bottomSecond = bottom.secondDerivs();

  /*
  //QUATERNION_T<Float> h(bottomCoeffs[0]);
  //for (int x = 1; x < bottomRoots + 1; x++)
  //  h += bottomCoeffs[x] * powers[x - 1];
  QUATERNION_T<Float> h(bottomCoeffs[bottomRoots]);
  for (int x = bottomRoots - 1; x >=0; x--)
    //h = h * point + bottomCoeffs[x];
    h.multiplyAdd(point, bottomCoeffs[x]);
  */
  QUATERNION_T<Float> h = bottom.evaluate(point);

  /*
  //QUATERNION_T<Float> hPrime(bottomDerivs[0]);
  //for (int x = 1; x < bottomRoots; x++)
  //  hPrime += bottomDerivs[x] * powers[x-1];
  QUATERNION_T<Float> hPrime(bottomDerivs[bottomRoots - 1]);
  for (int x = bottomRoots - 2; x >= 0; x--)
    //hPrime = hPrime * point + bottomDerivs[x];
    hPrime.multiplyAdd(point, bottomDerivs[x]);
  */
  QUATERNION_T<Float> hPrime = bottom.evaluateDerivative(point);
 
  // compute the rational derivatives 
  const QUATERNION_T<Float> numerator = gPrime * h - g * hPrime;
  const QUATERNION_T<Float> inverse = h.inverse();
  const QUATERNION_T<Float> inverseSq = inverse * inverse;
  p = g * inverse;
  pPrime = numerator * inverseSq;
  //const QUATERNION_T<Float> hSq = h * h;

  //p = g / h;
  //pPrime = numerator / hSq;
}

//////////////////////////////////////////////////////////////////////
// run a unit test for a known rational function
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::rationalTest()
{
  vector<float> coeffs;
  coeffs.push_back(-1.0);
  coeffs.push_back(0.0);
  coeffs.push_back(1.0);
  POLYNOMIAL_4D_T top = POLYNOMIAL_4D_T(coeffs);

  coeffs.clear();
  coeffs.push_back(0);
  coeffs.push_back(0);
  coeffs.push_back(1);
  POLYNOMIAL_4D_T bottom = POLYNOMIAL_4D_T(coeffs);

  QUATERNION_T<Float> point(1.23,4.56);

  QUATERNION_T<Float> g, gPrime, gPrimePrime;
  top.evaluateMultiple(point, g, gPrime, gPrimePrime);
  QUATERNION_T<Float> h, hPrime, hPrimePrime;
  bottom.evaluateMultiple(point, h, hPrime, hPrimePrime);

  QUATERNION_T<Float> rational = g / h;

  // compute (1 - 1 / z^2) directly
  QUATERNION_T<Float> ground = QUATERNION_T<Float>(1,0) - QUATERNION_T<Float>(1,0) / (point * point);

  QUATERNION_T<Float> diff = rational - ground;

  cout << " ground: " << ground << endl;
  cout << " computed: " << rational << endl;
  cout << " diff: " << diff.magnitude() << endl;

  QUATERNION_T<Float> numerator = gPrime * h - g * hPrime;
  QUATERNION_T<Float> hSq = h * h;
  QUATERNION_T<Float> pPrime = numerator / hSq;

  // compute (2 / z^3) directly
  ground = QUATERNION_T<Float>(2,0) / (point * point * point);
  diff = pPrime - ground;
  
  cout << " ground derivative: " << ground << endl;
  cout << " computed derivative: " << pPrime << endl;
  cout << " diff: " << diff.magnitude() << endl;

  QUATERNION_T<Float> hCubed = hSq * h;
  QUATERNION_T<Float> first = gPrimePrime * hSq;
  QUATERNION_T<Float> second = 2.0 * gPrime * hPrime * h;
  QUATERNION_T<Float> third = g * (2.0 * hPrime * hPrime - h * hPrimePrime);
  QUATERNION_T<Float> pPrimePrime = (first - second + third) / hCubed;

  // compute (-6 / z^4) directly
  ground = QUATERNION_T<Float>(-6,0) / (point * point * point * point);
  diff = pPrimePrime - ground;
  cout << " ground second derivative: " << ground << endl;
  cout << " computed second derivative: " << pPrimePrime << endl;
  cout << " diff: " << diff.magnitude() << endl;

}

//////////////////////////////////////////////////////////////////////
// add a new root
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::addRoot(const QUATERNION_T<Float>& newRoot)
{
  _roots.push_back(newRoot);
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  _rootPowers.push_back(1.0);

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// add a new root
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::addFrontRoot(const QUATERNION_T<Float>& newRoot)
{
  vector<QUATERNION_T<Float> > roots;
  roots.push_back(newRoot);
  for (unsigned int x = 0; x < _roots.size(); x++)
    roots.push_back(_roots[x]);
  _roots = roots;

  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  vector<Float> rootPowers;
  rootPowers.push_back(1.0);
  for (unsigned int x = 0; x < _rootPowers.size(); x++)
    rootPowers.push_back(_rootPowers[x]);
  _rootPowers = rootPowers;

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// add a new root
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::addRoot(const QUATERNION_T<Float>& newRoot, const Float& power)
{
  _roots.push_back(newRoot);
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  _rootPowers.push_back(power);

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// modify an existing root
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::modifyRoot(const int whichRoot, const QUATERNION_T<Float>& newRoot)
{
  assert(whichRoot < _totalRoots);

  if (whichRoot >= _totalRoots) return;

  _roots[whichRoot] = newRoot;
  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::write(FILE* file) const
{
  fwrite((void*)&_totalRoots, sizeof(int), 1, file);

  int size = _coeffs.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _coeffs[x].write(file);

  size = _derivs.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _derivs[x].write(file);

  size = _secondDerivs.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _secondDerivs[x].write(file);
  
  size = _roots.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _roots[x].write(file);

  size = _rootPowers.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
  {
    double power = _rootPowers[x];
    fwrite((void*)&power, sizeof(double), 1, file);
    //_rootPowers[x].write(file);
  }

  double powerScalar = _powerScalar;
  fwrite((void*)&powerScalar, sizeof(double), 1, file);
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::read(FILE* file)
{
  fread((void*)&_totalRoots, sizeof(int), 1, file);
  cout << " Reading in 4D polynomial with " << _totalRoots << " roots " << endl;

  int size;
  fread((void*)&size, sizeof(int), 1, file);
  _coeffs.resize(size);
  for (int x = 0; x < size; x++)
    _coeffs[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _derivs.resize(size);
  for (int x = 0; x < size; x++)
    _derivs[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _secondDerivs.resize(size);
  for (int x = 0; x < size; x++)
    _secondDerivs[x].read(file);
  
  fread((void*)&size, sizeof(int), 1, file);
  _roots.resize(size);
  for (int x = 0; x < size; x++)
    _roots[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _rootPowers.resize(size);
  for (int x = 0; x < _totalRoots; x++)
  {
    double power = 0;
    fread((void*)&power, sizeof(double), 1, file);
    _rootPowers[x] =power;
  }

#if 1
  double powerScalar;
  fread((void*)&powerScalar, sizeof(double), 1, file);
  _powerScalar = powerScalar;
#endif
}

//////////////////////////////////////////////////////////////////////
// Print polynomial to stream
//////////////////////////////////////////////////////////////////////
template <class Float>
ostream& operator<<(ostream &out, const POLYNOMIAL_4D_T<Float>& poly)
{
  for (unsigned int x = 0; x < poly.roots().size(); x++)
    out << " Root " << x << ": " << poly.roots()[x] << endl;
    /*
  for (unsigned int x = 0; x < poly.coeffs().size(); x++)
    out << " Coeff " << x << ": " << poly.coeffs()[x] << endl;
  for (unsigned int x = 0; x < poly.derivs().size(); x++)
    out << " Deriv " << x << ": " << poly.derivs()[x] << endl;
    */

  return out;
}

//////////////////////////////////////////////////////////////////////
// get the condition number of the polynomial
//////////////////////////////////////////////////////////////////////
template <class Float>
Float POLYNOMIAL_4D_T<Float>::conditionNumber()
{
  vector<Float> derivs;

  for (int x = 0; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> derivative = evaluateDerivative(_roots[x]);
    derivs.push_back(derivative.magnitude());
  }

  vector<QUATERNION_T<Float> > rootPowers = _roots;
  Float maxFound = 0;
  for (int i = 0; i < _totalRoots; i++)
  {
    for (int j = 0; j < _totalRoots; j++)
    {
      Float condition = (_coeffs[i] * rootPowers[j]).magnitude() / derivs[j];
      maxFound = (condition < maxFound) ? maxFound : condition;
    }

    // add another power to the roots
    for (int j = 0; j < _totalRoots; j++)
      rootPowers[j] = rootPowers[j] * _roots[j];
  }

  return maxFound;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::computeNestedCoeffs()
{
  _ws.clear();

  // the first w is the same as the root
  _ws.push_back(_roots[0]);

  /*
  for (unsigned int x = 1; x < _roots.size(); x++)
  {
    // evaluate the polynomial up until now (P_{k-1})
    QUATERNION_T<Float> Pkm1 = (_roots[x] - _ws[0]);
    assert(_ws.size() == x);
    for (unsigned int y = 1; y < _ws.size(); y++)
    {
      QUATERNION_T<Float> linear = (_roots[x] - _ws[y]);
      Pkm1 = linear * Pkm1;
    }

    QUATERNION_T<Float> w = Pkm1 * _roots[x] * Pkm1.inverse();
    _ws.push_back(w);
    cout << " inverse product: " << Pkm1 * Pkm1.inverse() << endl;
  }

  cout << " Computed ws: " << endl;
  for (unsigned int x = 0; x < _ws.size(); x++)
    cout << _ws[x] << endl;
    */

  //QUATERNION_T<Float> w = (_roots[1] * (_roots[1] - _roots[0])) * (_roots[1] - _roots[0]).inverse();
  //QUATERNION_T<Float> w = _roots[1] * ((_roots[1] - _roots[0]) * (_roots[1] - _roots[0]).inverse());
  QUATERNION_T<Float> w = _roots[1];
  _ws.push_back(w);
  w = _roots[2];
  _ws.push_back(w);
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateFactoredPositive(const QUATERNION_T<Float>& point) const
{
  QUATERNION_T<Float> result = point + _roots[0];

  for (int x = 1; x < _totalRoots; x++)
  {
    //cout  << " temporary: " << result << endl;
    //result = (point - _roots[x]) * result;
    result = result * (point + _roots[x]);
  }

  return result;

  //assert(_totalRoots == 2);
  //return (point - _ws[1]) * (point - _ws[0]);
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateFactoredDouble(const QUATERNION_T<Float>& point) const
{
  QUATERNION_T<Float> result = point - _roots[0];

  int power = 2;
  for (int x = 0; x < power; x++)
    result = result * result;

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> currentRoot = (point - _roots[x]);
    currentRoot = currentRoot * currentRoot;
    //result = result * currentRoot;
    result = currentRoot * result;
  }

  return result;

  //assert(_totalRoots == 2);
  //return (point - _ws[1]) * (point - _ws[0]);
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateScaledPowerFactored(const QUATERNION_T<Float>& point) const
{
  assert(_roots.size() == _rootPowers.size());
  QUATERNION_T<Float> result = point - _roots[0];
  result = result.pow(_powerScalar * _rootPowers[0]);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> term = (point - _roots[x]);
    result *= term.pow(_powerScalar * _rootPowers[x]);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluatePowerFactored(const QUATERNION_T<Float>& point) const
{
  assert(_roots.size() == _rootPowers.size());
  QUATERNION_T<Float> result = point - _roots[0];
  result = result.pow(_rootPowers[0]);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> term = (point - _roots[x]);
    result *= term.pow(_rootPowers[x]);
  }

  return result;
  /*
  assert(_roots.size() == _rootPowers.size());
  QUATERNION_T<Float> result = point - _roots[0];

  int power = _rootPowers[0];
  QUATERNION_T<Float> original = result;
  for (int x = 1; x < power; x++)
    result = result * original;

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> term = (point - _roots[x]);
    original = term;

    for (int y = 1; y < (int)_rootPowers[x]; y++)
      term = term * original;

    result = result * term;
  }

  return result;
  */
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation, but cache the
// multiplies
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluatePowerFactored(const QUATERNION_T<Float>& point, vector<QUATERNION_T<Float> >& forward, vector<QUATERNION_T<Float> >& backward) const
{
  assert(_roots.size() == _rootPowers.size());
  vector<QUATERNION_T<Float> > powers;
  forward.clear();
  backward.clear();

  QUATERNION_T<Float> result = point - _roots[0];
  result = result.pow(_rootPowers[0]);
  powers.push_back(result);
  forward.push_back(result);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> term = (point - _roots[x]).pow(_rootPowers[x]);
    powers.push_back(term);
    result *= term;
    forward.push_back(result);
  }

  // cache out the multiplies
  QUATERNION_T<Float> current = powers[_totalRoots - 1];
  backward.resize(_totalRoots);
  backward[_totalRoots - 1] = current;
  for (int x = _totalRoots - 2; x >= 0; x--)
  {
    current = powers[x] * current;
    backward[x] = current;
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation, but cache the
// multiplies
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateScaledPowerFactored(const QUATERNION_T<Float>& point, vector<QUATERNION_T<Float> >& forward, vector<QUATERNION_T<Float> >& backward) const
{
  assert(_roots.size() == _rootPowers.size());
  vector<QUATERNION_T<Float> > powers;
  forward.clear();
  backward.clear();

  QUATERNION_T<Float> result = point - _roots[0];
  result = result.pow(_powerScalar * _rootPowers[0]);
  powers.push_back(result);
  forward.push_back(result);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> term = (point - _roots[x]).pow(_powerScalar * _rootPowers[x]);
    powers.push_back(term);
    result *= term;
    forward.push_back(result);
  }

  // cache out the multiplies
  QUATERNION_T<Float> current = powers[_totalRoots - 1];
  backward.resize(_totalRoots);
  backward[_totalRoots - 1] = current;
  for (int x = _totalRoots - 2; x >= 0; x--)
  {
    current = powers[x] * current;
    backward[x] = current;
  }
  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateFactored(const QUATERNION_T<Float>& point) const
{
  QUATERNION_T<Float> result = point - _roots[0];

  for (int x = 1; x < _totalRoots; x++)
    result = result * (point - _roots[x]);

  return result;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::evaluateFactoredRational(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime)
{
  const vector<QUATERNION_T<Float> > topRoots = top.roots();
  assert(topRoots.size() > 0);

  // assume ordering is (x - _roots[0]) * (x - _roots[1]) ...
  QUATERNION_T<Float> g = (point - topRoots[0]);
  for (unsigned int x = 1; x < topRoots.size(); x++)
    g = g * (point - topRoots[x]);

  const vector<QUATERNION_T<Float> > bottomRoots = bottom.roots();
  assert(bottomRoots.size() > 0);

  // assume ordering is (x - _roots[0]) * (x - _roots[1]) ...
  QUATERNION_T<Float> h = (point - bottomRoots[0]);
  for (unsigned int x = 1; x < bottomRoots.size(); x++)
    h = h * (point - bottomRoots[x]);
  
  QUATERNION_T<Float> gPrime = top.evaluateFactoredDerivative(point);
  QUATERNION_T<Float> hPrime = bottom.evaluateFactoredDerivative(point);
  
  const QUATERNION_T<Float> numerator = gPrime * h - g * hPrime;
  const QUATERNION_T<Float> inverse = h.inverse();
  const QUATERNION_T<Float> inverseSq = inverse * inverse;
  p = g * inverse;
  pPrime = numerator * inverseSq;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::evaluateFactoredQuadratic(const POLYNOMIAL_4D_T& top, const POLYNOMIAL_4D_T& bottom, const QUATERNION_T<Float>& point, QUATERNION_T<Float>& p, QUATERNION_T<Float>& pPrime)
{
  // assume ordering is (x - _roots[0]) * (x - _roots[1]) ...
  const vector<QUATERNION_T<Float> > topRoots = top.roots();
  assert(topRoots.size() > 0);
  QUATERNION_T<Float> g = (point - topRoots[0]) * (point - topRoots[1]);

  const vector<QUATERNION_T<Float> > bottomRoots = bottom.roots();
  assert(bottomRoots.size() > 0);
  QUATERNION_T<Float> h = (point - bottomRoots[0]) * (point - bottomRoots[1]);
  
  QUATERNION_T<Float> gPrime = (point - topRoots[0]) + (point - topRoots[1]);
  QUATERNION_T<Float> hPrime = (point - bottomRoots[0]) + (point -  bottomRoots[1]);
  
  const QUATERNION_T<Float> numerator = gPrime * h - g * hPrime;
  const QUATERNION_T<Float> inverse = h.inverse();
  const QUATERNION_T<Float> inverseSq = inverse * inverse;
  p = g * inverse;
  pPrime = numerator * inverseSq;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::evaluateFactoredDerivative(const QUATERNION_T<Float>& point) const
{
  QUATERNION_T<Float> final;

  //  which one to knock out?
  for (unsigned int x = 0; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> knockedOut(1,0,0,0);
    for (unsigned int y = 0; y < _totalRoots; y++)
    {
      if (y == x) continue;
      knockedOut = knockedOut * (point - _roots[y]);
    }

    final += knockedOut;
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// sum of all the root coefficients
//////////////////////////////////////////////////////////////////////
template <class Float>
Float POLYNOMIAL_4D_T<Float>::rootSum() const
{
  Float final = 0;
  for (unsigned int x = 0; x < _totalRoots; x++)
  {
    final += _roots[x].x();
    final += _roots[x].y();
    final += _roots[x].z();
    final += _roots[x].w();
  }
  return final; 
}

//////////////////////////////////////////////////////////////////////
// overloaded operators
//////////////////////////////////////////////////////////////////////
template <class Float>
POLYNOMIAL_4D_T<Float>& POLYNOMIAL_4D_T<Float>::operator*=(const Float& alpha)
{
  for (unsigned int x = 0; x < _totalRoots; x++)
    _roots[x] *= alpha;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// overloaded operators
//////////////////////////////////////////////////////////////////////
template <class Float>
POLYNOMIAL_4D_T<Float>& POLYNOMIAL_4D_T<Float>::operator-=(const VEC3F& v)
{
  QUATERNION_T<Float> q(v[0], v[1], v[2], 0);

  for (unsigned int x = 0; x < _totalRoots; x++)
    _roots[x] -= q;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// overloaded operators
//////////////////////////////////////////////////////////////////////
template <class Float>
POLYNOMIAL_4D_T<Float>& POLYNOMIAL_4D_T<Float>::operator+=(const VEC3F& v)
{
  QUATERNION_T<Float> q(v[0], v[1], v[2], 0);

  for (unsigned int x = 0; x < _totalRoots; x++)
    _roots[x] += q;

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::translateExceptFirst(const VEC3F& v)
{
  QUATERNION_T<Float> q(v[0], v[1], v[2], 0);

  for (unsigned int x = 1; x < _totalRoots; x++)
    _roots[x] += q;
}

//////////////////////////////////////////////////////////////////////
// test out taking the derivative of a power
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::testSingleDerivative()
{
  int res = 100;
  vector<QUATERNION_T<Float> > roots;
  roots.push_back(QUATERNION_T<Float>(1.0, 2.0, 3.0, 0));
  POLYNOMIAL_4D_T polynomial(roots);

  vector<Float> xs;
  vector<Float> ys;
  vector<Float> squares;
  vector<Float> gradients;
  vector<Float> temps;

  //QUATERNION_T<Float> root(1.0, 2.0, 3.0, 4.0);
  QUATERNION_T<Float> root(2, 2,3,4);
  //QUATERNION_T<Float> root(2, 10,2,10);
  //root.normalize();
#if 1
  //for (Float x = -10; x < 10; x += 0.25)
  //for (Float x = -10; x < 10; x += 0.1)
  for (Float x = -10; x < 10; x += 0.1)
  {
    QUATERNION_T<Float> value = root.pow(x);
    QUATERNION_T<Float> derivative = value * root.log();

    Float gradient = (1.0 / sqrt(value.dot(value))) * value.dot(derivative);
    assert(!(value.anyNans()) && !(derivative.anyNans()));

    xs.push_back(x);
    ys.push_back(sqrt(value.dot(value)));
    gradients.push_back(gradient);
  }
#else
  for (Float x = -10; x < 10; x += 0.25)
  {
    QUATERNION_T<Float> value = root.pow(x);

    // separate into two components
    Float magnitude = root.magnitude();
    QUATERNION_T<Float> qhat = root;
    qhat.normalize();

    QUATERNION_T<Float> qhatPow = qhat.pow(x);

    QUATERNION_T<Float> qhatDerivative = qhatPow * qhat.log();
    QUATERNION_T<Float> derivative = pow(magnitude, x) * qhatDerivative + x * pow(magnitude, x - 1) * qhatPow;
    //QUATERNION_T<Float> derivative = pow(magnitude, x) * qhatDerivative;

    Float temp = value.dot(derivative);

    Float square = value.dot(derivative) + derivative.dot(value);
    //Float gradient = sqrt(value.dot(derivative) + derivative.dot(value));
    Float gradient = value.dot(derivative) + derivative.dot(value);

    xs.push_back(x);
    //ys.push_back(value.magnitude());
    ys.push_back(value.dot(value));
    gradients.push_back(gradient);
    squares.push_back(square);
    temps.push_back(temp);
  }
#endif

  // try to do some convergence
  Float dx = 0.1;
  Float fixedPoint = 2.0;

  QUATERNION_T<Float> value = root.pow(fixedPoint);
  QUATERNION_T<Float> derivative = value * root.log();

  //Float gradient = value.dot(derivative) + derivative.dot(value);
  //Float gradient = 0.5 * pow(value.dot(value), (Float)-0.5) * (value.dot(derivative) + derivative.dot(value));
  Float gradient = pow(value.dot(value), (Float)-0.5) * value.dot(derivative);
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Computed analytic: " << gradient << endl;

  vector<Float> dxs;
  vector<Float> diffs;

  for (int x = 0; x < 10; x++)
  {
    //Float R = value.dot(value);
    Float R = sqrt(value.dot(value));

    QUATERNION_T<Float> displaced = root.pow(fixedPoint + dx);
    //Float Rdx = displaced.dot(displaced);
    Float Rdx = sqrt(displaced.dot(displaced));
    
    Float numerical = (Rdx - R) / dx;
    //cout << " Numerical for dx " << dx << ": " << numerical << endl;
    dx *= 0.1;

    dxs.push_back(dx);
    diffs.push_back((numerical - gradient) / gradient);

    cout << " diff: " << (numerical - gradient) / gradient << "\t dx: " << dx << "\t numerical: " << numerical << endl;
  }
  cout << " dx = " << VECTOR(dxs) << endl;
  cout << " diff = " << VECTOR(diffs) << endl;

  /*
  QUATERNION_T<Float> root(20.0, 20.0, 20.0, 20.0);
  root.normalize();

  //for (Float x = -10; x < 10; x += 0.1)
  //for (Float x = -10; x < 10; x += 0.1)
  //for (Float x = -10; x < 10; x +=1)
  for (Float x = -10; x < 0; x += 0.25)
  //for (Float x = 0; x < 10; x += 0.25)
  //for (Float x = -10; x < 0; x += 0.5)
  {
    //QUATERNION_T<Float> root(1.0, 2.0, 3.0, 4.0);
    QUATERNION_T<Float> value = root.pow(x);

    // separate into two components
    Float magnitude = root.magnitude();
    QUATERNION_T<Float> qhat = root;
    qhat.normalize();

    QUATERNION_T<Float> qhatPow = qhat.pow(x);

    QUATERNION_T<Float> qhatDerivative = qhatPow * qhat.log();
    QUATERNION_T<Float> derivative = pow(magnitude, x) * qhatDerivative + x * pow(magnitude, x - 1) * qhatPow;
    //QUATERNION_T<Float> derivative = pow(magnitude, x) * qhatDerivative;

    Float temp = value.dot(derivative);

    Float square = value.dot(derivative) + derivative.dot(value);
    //Float gradient = sqrt(value.dot(derivative) + derivative.dot(value));
    Float gradient = value.dot(derivative) + derivative.dot(value);

    xs.push_back(x);
    //ys.push_back(value.magnitude());
    ys.push_back(value.dot(value));
    gradients.push_back(gradient);
    squares.push_back(square);
    temps.push_back(temp);
  }
  */
  /*
  VECTOR xsMatlab(xs);
  VECTOR ysMatlab(ys);
  VECTOR gradientsMatlab(gradients);
  VECTOR squaresMatlab(squares);
  VECTOR tempsMatlab(temps);
  xsMatlab.writeMatlab("xPower1D.m", "x");
  ysMatlab.writeMatlab("yPower1D.m", "y");
  gradientsMatlab.writeMatlab("gradientsPower1D.m", "gradient");
  squaresMatlab.writeMatlab("squaresPower1D.m", "squares");
  tempsMatlab.writeMatlab("tempsPower1D.m", "temps");
  */

  return;
}

//////////////////////////////////////////////////////////////////////
// take the derivative with respect to a root
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::powerDerivative(const QUATERNION_T<Float>& point, const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  assert(whichRoot >= 0);
  assert(whichRoot < _totalRoots);

  QUATERNION_T<Float> left(1,0,0,0);
  QUATERNION_T<Float> right(1,0,0,0);
  
  for (int x = 0; x < whichRoot; x++)
    left *= (point - _roots[x]).pow(_rootPowers[x]);
  
  for (int x = whichRoot + 1; x < _totalRoots; x++)
    right *= (point - _roots[x]).pow(_rootPowers[x]);

  // get the derivative of the root in question
  QUATERNION_T<Float> base = (point - _roots[whichRoot]);
  QUATERNION_T<Float> power = base.pow(_rootPowers[whichRoot]);
  QUATERNION_T<Float> middle = power * base.log();

  return left * middle * right;
}

//////////////////////////////////////////////////////////////////////
// take the derivative with respect to a root
//////////////////////////////////////////////////////////////////////
template <class Float>
QUATERNION_T<Float> POLYNOMIAL_4D_T<Float>::inversePowerDerivative(const QUATERNION_T<Float>& point, const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  assert(whichRoot >= 0);
  assert(whichRoot < _totalRoots);

  QUATERNION_T<Float> left(1,0,0,0);
  QUATERNION_T<Float> right(1,0,0,0);
  
  for (int x = 0; x < whichRoot; x++)
    left *= (point - _roots[x]).pow(_rootPowers[x]);
  
  for (int x = whichRoot + 1; x < _totalRoots; x++)
    right *= (point - _roots[x]).pow(_rootPowers[x]);

  // get the derivative of the root in question
  QUATERNION_T<Float> base = (point - _roots[whichRoot]);
  QUATERNION_T<Float> power = base.pow(_rootPowers[whichRoot]);
  QUATERNION_T<Float> middle = power * base.log().inverse();

  return QUATERNION_T<Float>(-1,0,0,0) * left * middle * right;
}

//////////////////////////////////////////////////////////////////////
// compute the gradient with respect to each root
//////////////////////////////////////////////////////////////////////
template <class Float>
VECTOR POLYNOMIAL_4D_T<Float>::powerGradient(const QUATERNION_T<Float>& point)
{
  VECTOR gradient(_totalRoots);

  // TODO: tons of tree-caching optimizations possible here
  QUATERNION_T<Float> value = evaluatePowerFactored(point);
  for (int x = 0; x < _totalRoots; x++)
  {
    QUATERNION_T<Float> derivative = powerDerivative(point, x);
    gradient[x] = pow(value.dot(value), (Float)-0.5) * value.dot(derivative);
  }

  return gradient;
}

//////////////////////////////////////////////////////////////////////
// test out taking the derivative of a power
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::testPolynomialDerivative()
{
  vector<QUATERNION_T<Float> > roots;
  roots.push_back(QUATERNION_T<Float>(1.0, 2.0, 3.0, 4.0));
  roots.push_back(QUATERNION_T<Float>(4.0, 5.0, 6.0, 7.0));
  roots.push_back(QUATERNION_T<Float>(0.8, 0.9, 0.10, 0.11));

  Float fixedPower = 3;
  vector<Float> powers;
  powers.push_back(2);
  powers.push_back(fixedPower);
  powers.push_back(4);

  POLYNOMIAL_4D_T polynomial(roots, powers);

  //QUATERNION_T<Float> fixedPoint(1,1,1,1);
  QUATERNION_T<Float> fixedPoint(-1,0.3,-0.1234,2);

  // analytical derivative
  int whichRoot = 1;
  //QUATERNION_T<Float> value = polynomial.evaluatePowerFactored(fixedPoint);
  //QUATERNION_T<Float> derivative = polynomial.powerDerivative(fixedPoint, whichRoot);
  QUATERNION_T<Float> value = fixedPoint * polynomial.evaluatePowerFactored(fixedPoint);
  QUATERNION_T<Float> derivative = fixedPoint * polynomial.powerDerivative(fixedPoint, whichRoot);

  Float gradient = pow(value.dot(value), (Float)-0.5) * value.dot(derivative);
  //Float gradient = 2.0 * value.dot(derivative);

  cout << " Analytical derivative: " << derivative << endl;
  cout << " Analytical gradient:   " << gradient << endl;

  // does the numerical converge?
  vector<Float> dxs;
  vector<Float> diffs;
  Float dx = 0.1;

  vector<Float>& rootPowers = polynomial.powersMutable();

  Float R = sqrt(value.dot(value));
  for (int x = 0; x < 10; x++)
  {
    const Float originalPower = rootPowers[whichRoot];
    rootPowers[whichRoot] += dx;
    //QUATERNION_T<Float> displaced = polynomial.evaluatePowerFactored(fixedPoint);
    QUATERNION_T<Float> displaced = fixedPoint * polynomial.evaluatePowerFactored(fixedPoint);
    rootPowers[whichRoot] = originalPower;

    Float Rdx = sqrt(displaced.dot(displaced));
    
    Float numerical = (Rdx - R) / dx;
    dx *= 0.1;

    dxs.push_back(dx);
    diffs.push_back((numerical - gradient) / gradient);

    cout << " diff: " << (numerical - gradient) / gradient << "\t dx: " << dx << "\t numerical: " << numerical << endl;
  }

  VECTOR fullGradient = polynomial.powerGradient(fixedPoint);
  cout << " Full analytical gradient: " << fullGradient << endl;

  return;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::testBulkDerivative()
{
  cout << " ========================================================================= " << endl;
  cout << "  TESTING DERIVATIVE OF ALL POWERS AT ONCE " << endl;
  cout << " ========================================================================= " << endl;
  vector<QUATERNION_T<Float> > roots;
  roots.push_back(QUATERNION_T<Float>(1.0, 2.0, 3.0, 4.0));
  roots.push_back(QUATERNION_T<Float>(4.0, 5.0, 6.0, 7.0));
  roots.push_back(QUATERNION_T<Float>(0.8, 0.9, 0.10, 0.11));

  Float fixedPower = 3;
  vector<Float> powers;
  powers.push_back(1);
  powers.push_back(1);
  powers.push_back(1);

  POLYNOMIAL_4D_T polynomial(roots, powers);

  //QUATERNION_T<Float> fixedPoint(1,1,1,1);
  QUATERNION_T<Float> fixedPoint(-1,0.3,-0.1234,2);

  QUATERNION_T<Float> value = fixedPoint * polynomial.evaluatePowerFactored(fixedPoint);

  // analytical derivative
  Float totalGradient = 0;
  for (int x = 0; x < powers.size(); x++)
  {
    QUATERNION_T<Float> derivative = fixedPoint * polynomial.powerDerivative(fixedPoint, x);
    Float gradient = pow(value.dot(value), (Float)-0.5) * value.dot(derivative);

    totalGradient += gradient;
  }

  //cout << " Analytical derivative: " << derivative << endl;
  cout << " Analytical gradient:   " << totalGradient << endl;

  // does the numerical converge?
  //vector<Float> dxs;
  //vector<Float> diffs;
  Float dx = 0.1;

  vector<Float>& rootPowers = polynomial.powersMutable();

  Float R = sqrt(value.dot(value));
  for (int x = 0; x < 10; x++)
  {
    // do the dx displacement
    vector<Float> originalPowers;
    for (int y = 0; y < powers.size(); y++)
    {
      originalPowers.push_back(rootPowers[y]);
      rootPowers[y] += dx;
    }
    QUATERNION_T<Float> displaced = fixedPoint * polynomial.evaluatePowerFactored(fixedPoint);

    // wind things back
    for (int y = 0; y < powers.size(); y++)
      rootPowers[y] = originalPowers[y];

    Float Rdx = sqrt(displaced.dot(displaced));
    Float numerical = (Rdx - R) / dx;
    dx *= 0.1;

    //dxs.push_back(dx);
    //diffs.push_back((numerical - totalGradient) / totalGradient);

    cout << " diff: " << (numerical - totalGradient) / totalGradient << "\t dx: " << dx << "\t numerical: " << numerical << endl;
  }

  //VECTOR fullGradient = polynomial.powerGradient(fixedPoint);
  //cout << " Full analytical gradient: " << fullGradient << endl;

  return;
}

//////////////////////////////////////////////////////////////////////
// stomp everything
//////////////////////////////////////////////////////////////////////
template <class Float>
void POLYNOMIAL_4D_T<Float>::clear()
{
  _totalRoots = 0;
  _coeffs.clear();
  _derivs.clear();
  _secondDerivs.clear();
  _roots.clear();
  _rootPowers.clear();
  _ws.clear();
}

typedef POLYNOMIAL_4D_T<float> POLYNOMIAL_4DF;
typedef POLYNOMIAL_4D_T<double> POLYNOMIAL_4DD;
typedef POLYNOMIAL_4D_T<long double> POLYNOMIAL_4DE;
//typedef POLYNOMIAL_4D_T<mpreal> POLYNOMIAL_4DMP;

#endif
