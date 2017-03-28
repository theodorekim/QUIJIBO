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
#ifndef OPTIMIZE_3D_H
#define OPTIMIZE_3D_H

#include <cmath>
#include <string>
#include <iostream>

#include <FIELD_3D.h>
#include <POLYNOMIAL_4D.h>
#include "POLYNOMIAL_CACHE.h"
#include "QUATERNION_CACHE.h"

#include <POLYNOMIAL_CACHE_T.h>
#include <POLYNOMIAL_4D_T.h>

using namespace std;

#ifndef POLYNOMIALVIEW3D
#define POLYNOMIALVIEW3D(x) OPTIMIZE_3D::polynomialViewer(x); sleep(1);
#endif
#ifndef LOGPOLYNOMIALVIEW3D
#define LOGPOLYNOMIALVIEW3D(x,y) OPTIMIZE_3D::logPolynomialViewer(x,y); sleep(1);
#endif

class OPTIMIZE_3D {

public:
  OPTIMIZE_3D();
  OPTIMIZE_3D(int xRes, int yRes, int zRes, const VEC3F& center, const VEC3F& lengths);
  virtual ~OPTIMIZE_3D();

  // accessors
  FIELD_3D& distanceField()  { return _distanceField; };
  FIELD_3D& curvatureField() { return _curvatureField; };
  FIELD_3D& fractal()   { return _fractal; };
  FIELD_3D& auxiliary() { return _auxiliary; };
  POLYNOMIAL_4D& top()      { return _top; };
  POLYNOMIAL_4D& bottom()   { return _bottom; };
  POLYNOMIAL_4D top() const    { return _top; };
  POLYNOMIAL_4D bottom() const { return _bottom; };
  //virtual const int totalDOFs()     { return 4 * (_top.totalRoots() + _bottom.totalRoots()); };
  virtual const int totalDOFs()     { return 4 * _top.totalRoots(); };
  VEC3F& center() { return _center; };
  VEC3F& lengths() { return _lengths; };
  //VEC3F& flowCenter() { return _flowCenter; };
  //Real& flowRadius() { return _flowRadius; };
  vector<VECTOR>& history() { return _historyTop; };
  Real dx() { return _lengths[0] / _xRes; };
  int& maxIterations() { return _maxIterations; };
  int& xRes() { return _xRes; };
  int& yRes() { return _yRes; };
  int& zRes() { return _zRes; };
  int& binaryBandwidth() { return _binaryBandwidth; };
  int& curvatureBandwidth() { return _curvatureBandwidth; };
  Real& expScaling() { return _expScaling; };
  const Real expScaling() const { return _expScaling; };
  string& lastLoadedFilename() { return _lastLoadedFilename; };
  Real& powerDx() { return _powerDx; };
  FIELD_3D& inverseDistanceField()  { return _inverseDistanceField; };
  FIELD_3D& weightField()  { return _weightField; };
  FIELD_3D& overflowed() { return _overflowed; };
  bool& adaptiveFired() { return _adaptiveFired; };
  Real& slice() { return _quaternionSlice; };

  void setScoreFields(const FIELD_3D& distanceField, const FIELD_3D& curvatureField);
  void setScoreFields(const FIELD_3D& distanceField, const FIELD_3D& curvatureField, const FIELD_3D& inverseField);
  void setLengthCenter(const Real& newLength, const VEC3F& newCenter);

  // get all the roots, all at once
  virtual VECTOR getRoots() const;
  virtual VECTOR getRationalRoots() const;
  virtual VECTOR getBottomRoots() const;

  // set all the roots, all at once
  virtual void setRoots(const VECTOR& newRoots);
  virtual void setRationalRoots(const VECTOR& newRoots);
  virtual void setBottomRoots(const VECTOR& newRoots);

  // drawing functions
  void drawRoots();
  void drawPowerRoots();
  void drawRoots2D();
  void drawRootsSliceZ(const int zSlice);
  void drawRootsSliceY(const int ySlice);
  /*
  void drawGrid();
  void drawDirections();
  void drawCenterOfMass();
  */

  // energy derivatives
  virtual VECTOR computeCenteredGradient();
  virtual VECTOR computeCenteredGradientCUDA();
  virtual VECTOR computeBandedCenteredGradient();
  virtual VECTOR computeBandedCenteredGradientCUDA();

  VECTOR computeLogCenteredGradient();
  VECTOR computeLogBandedCenteredGradient();
  VECTOR computeLogRationalCenteredGradient();
  VECTOR computeLogRationalBottomCenteredGradient();
  VECTOR computePowerCenteredGradient();

  // self derivatives
  VECTOR computeCenteredSelfGradient();
  VECTOR computeCenteredSelfGradientCUDA();

  // compute the derivative with respect to scaling and translation
  VECTOR computeScaleTranslateGradient();
  VECTOR computeLogScaleTranslateGradient();
  VECTOR computeScaleTranslateGradientCUDA();

  // recompute _fractal
  virtual void computePolynomial() { computePolynomial(_top, _bottom, _fractal); };
  void computePolynomial(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal);
  /*
  void computePolynomial(const int xRes, const int yRes, const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& insideOutside, FIELD_3D& iterationCount);
  void computePolynomial(FIELD_3D& insideOutside, FIELD_3D& iterationCount);
  void computePolynomialBanded(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, const int bandwidth);
  */
  void computeFactoredPolynomial();

  void computeMap() { computeMap(_top, _fractal, _auxiliary, false); };
  void computeMap(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, FIELD_3D& auxiliary, bool verbose = false);
  void computeMapCUDA() { computeMapCUDA(_top, _fractal, _auxiliary); };
  void computeMapCUDA(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, FIELD_3D& auxiliary);
  void computeKittyMap();
  void computeLogMap(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, bool verbose = false);
  void computeLogMap(bool verbose = false) { computeLogMap(_top, _fractal, verbose); };
  void computeLogRationalMap(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false);
  void computeLogRationalMap(bool verbose = false) { computeLogRationalMap(_top, _bottom, _fractal, verbose); };
  void computeLogRationalBottomMap(const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false);
  void computeLogRationalBottomMap(bool verbose = false) { computeLogRationalBottomMap(_bottom, _fractal, verbose); };
  void computeLogBandedMap(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, bool verbose = false);
  void computeLogBandedMap(bool verbose = false) { computeLogBandedMap(_top, _fractal, verbose); };
  FIELD_2D computeKittyMap2D(int xRes, int yRes);
  void computeLogPolynomial();

  void computeLogRationalMapDistanceEstimated(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false);
  void computeLogRationalMapDistanceEstimated(bool verbose = false) { computeLogRationalMapDistanceEstimated(_top, _bottom, _fractal, verbose); };

  void computeLogPowerRationalMapCached(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false);
  void computeLogScaledPowerRationalMapCached(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false);
  void computeLogScaledPowerRationalMap(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false);
  void computeLogPowerRationalMap(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false);
  //void computeLogPowerRationalMap(bool verbose = false) { computeLogPowerRationalMap(_top, _bottom, _fractal, verbose); };
  void computeLogPowerRationalMap(bool verbose = false) { computeLogPowerRationalMapCached(_top, _bottom, _fractal, verbose); };

  // THIS ONE
  void computeLogScaledPowerRationalMap(bool verbose = false) { computeLogScaledPowerRationalMap(_top, _bottom, _fractal, verbose); };
  void computeAdaptivePowerRationalMap(bool verbose = false);
  void computeAdaptivePowerRationalMap(int powerBits, int mantissaBits, bool verbose);

  virtual Real computeScore(const VECTOR& roots);
  virtual Real computeScoreCUDA(const VECTOR& roots);
  virtual Real computeBandedScore(const VECTOR& roots);
  virtual Real computeBandedScoreCUDA(const VECTOR& roots);
  virtual Real computeLogRationalScore(const VECTOR& roots);
  virtual Real computeLogRationalBottomScore(const VECTOR& roots);
  Real computeLogRationalBottomScoreVarying(const VECTOR& roots);
  Real computeLogRationalBottomScoreVarying();
  Real computeLogPowerRationalScore();
  Real computeLogScaledPowerRationalScore();
  Real computeLogPowerRationalScoreExpGradient();
  Real computeLogPowerScaledRationalScoreExpGradient();
  virtual Real computePowerScore(const VECTOR& roots);
 
  // adaptive precision (AP) versions of various functions
  Real computeScoreAP();                                   // computeLogScaledPowerRationalScore() 
  void computeRationalMapAP(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose = false); // computeLogScaledPowerRationalMapCached
  Real computeExpGradientAP();                       // computeLogPowerScaledRationalScoreExpGradient
  Real computeBottomBulkScalarGradientAP();          // computeBottomScaledPowerGradient
  Real computeTopBulkScalarGradientAP();             // computeTopScaledPowerGradient
  Real computeBottomGradientAP(const int whichRoot); // computeBottomScaledPowerGradient
  Real computeTopGradientAP(const int whichRoot);    // computeTopScaledPowerGradient

  Real computeMaximumScaledPowerRationalScore();

  //void initializeCircleRoots(const int totalTopRoots, const int totalBottomRoots, const bool probe = false);
  void initializeDebugCircleRoots(const int totalTopRoots, const int totalBottomRoots);

  void initRandomSphere();
  void initPolynomial();
  virtual void addToHistory();
  virtual void clearHistory() { _historyTop.clear(); };
  virtual void write(const string& filename) const;
  void read(const string& filename);
  void readDistanceField(const string& filename);
  void initCube();
  void initMap();

  Real computeScore() { return computeScore(_top, _fractal); };
  Real computeScoreCUDA() { return computeScoreCUDA(_top, _fractal); };
  Real computeBandedScore() { return computeBandedScore(_top, _fractal); };
  Real computeBandedScoreCUDA() { return computeBandedScoreCUDA(_top, _fractal); };
  Real computeScore(const POLYNOMIAL_4D& top, FIELD_3D& fractal);
  Real computeScoreCUDA(const POLYNOMIAL_4D& top, FIELD_3D& fractal);
  Real computeBandedScore(const POLYNOMIAL_4D& top, FIELD_3D& fractal);
  Real computeBandedScoreCUDA(const POLYNOMIAL_4D& top, FIELD_3D& fractal);

  // compute log versions of the score
  Real computeLogScore(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal);
  Real computeLogScore() { return computeLogScore(_top, _fractal); };
  Real computeLogScore(const VECTOR& roots) { setRoots(roots); return computeLogScore(); };
  Real computeLogBandedScore(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal);
  Real computeLogBandedScore() { return computeLogBandedScore(_top, _fractal); };
  Real computeLogBandedScore(const VECTOR& roots) { setRoots(roots); return computeLogBandedScore(); };
  Real computeLogRationalScore(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal);
  Real computeLogRationalScore() { return computeLogRationalScore(_top, _bottom, _fractal); };
  Real computeLogRationalBottomScore(const POLYNOMIAL_4D& bottom, FIELD_3D& fractal);
  Real computeLogRationalBottomScore() { return computeLogRationalBottomScore(_bottom, _fractal); };

  // compute a power score
  Real computePowerScore(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal);
  Real computePowerScore() { return computePowerScore(_top, _bottom, _fractal); };

  Real computeBandedBinaryScore() { return computeBinaryBandedScore(_fractal); };
  Real computeBandedCurvatureScore() { return computeCurvatureBandedScore(_fractal); };

  // compute the center of a signed distance field
  VEC3F computeCenter(const FIELD_3D& field);
  void computeCenterAndRadius(const FIELD_3D& field, VEC3F& center, Real& radius);
  void computeDistanceCenterAndRadius(VEC3F& center, Real& radius) { computeCenterAndRadius(_distanceField, center, radius); };
  void computeFractalCenterAndRadius(VEC3F& center, Real& radius) { computeMap(); computeCenterAndRadius(_fractal, center, radius); };
 
  // compute the error on the outside (assuming the convex hull has been built already)
  FIELD_3D computeExteriorErrorField();

  // do some CUDA comparisons
  void debugCUDA();

  // list of roots within this z slice
  vector<int> visibleRootsZ(const int zSlice);
  vector<int> visibleRootsY(const int ySlice);

  // verify that CUDA is still consistent
  bool verifyCUDA();

  static void polynomialViewer(const OPTIMIZE_3D& optimize3D);
  static void logPolynomialViewer(const OPTIMIZE_3D& optimize3D, const int expPow);

  // IO for the scale translate optimization
  void writeScaleTranslate(const string& path, int subdivisions, bool includeOrigin) const;
  bool readScaleTranslate(const string& path, int subdivisions, bool includeOrigin);

  // compute an error field along the fractal borer
  FIELD_3D errorField();

  // scale the entire optimization problem
  void scaleEverything(const Real& scale);
  
  // translate the entire optimization problem
  void translateEverything(const VEC3F& translate);
  void translateEverythingExceptFirst(const VEC3F& translate);

  // load up the kitty as the polynomial
  void loadKitty();

  // get a vector of all the root exponents
  VECTOR getPowers() const;
  void setPowers(const VECTOR& powers);

  // take a first shot at estimating the conformal radius
  Real estimateConformalRadius();

  // get the inverse distance score
  Real computeInverseDistanceScore() { return computeInverseDistanceScore(_fractal); };

  // plot one of the scoring functions with respect to one variable
  void plotExpEnergy1D();
  static void plotEnergy2D(OPTIMIZE_3D& input);

  // draw the output of a score computation
  void drawInverseDistanceScore();

  // debug the derivative of a rational function
  void debugRationalTopDerivative(int whichRoot);
  void debugRationalLogTopDerivative(int whichRoot);
  void debugRationalBottomDerivative(int whichRoot);
  void debugRationalTopScoreDerivative(int whichRoot);
  void debugRationalBottomScoreDerivative(int whichRoot);
  void debugRationalTopSummedScoreDerivative(int whichRoot);
  void debugRationalBottomSummedScoreDerivative(int whichRoot);
  void debugConformalRadiusDerivative();
  void debugRationalBulkTopDerivative();
  void debugRationalBulkBottomDerivative();

  // try with a global power scaling
  void debugRationalTopScaledPowerDerivative(int whichRoot);
  void debugRationalBottomScaledPowerDerivative(int whichRoot);
  
  // try with a global power scaling, get the derivative of
  // _powerScalar itself
  void debugRationalTopScaledPowerDerivative();
  void debugRationalBottomScaledPowerDerivative();

  // debug one of the Hessian entries
  void debugHessianConformalConformal();
  void debugHessianTopConformal(const int whichRoot);
  void debugHessianBottomConformal(const int whichRoot);
  void debugHessianTopDiagonal(const int whichRoot);
  void debugHessianBottomDiagonal(const int whichRoot);
  void debugHessianMixedTop(const int firstRoot, const int secondRoot);
  void debugHessianMixedBottom(const int firstRoot, const int secondRoot);
  void debugHessianTopBottom(const int topRoot, const int bottomRoot);

  // compute the derivative of a root power
  Real computeTopPowerGradient(const int whichRoot) const;
  Real computeBottomPowerGradient(const int whichRoot) const;
 
  // derivative of particular root power, incorporating _powerScalar 
  Real computeTopScaledPowerGradient(const int whichRoot) const;
  Real computeBottomScaledPowerGradient(const int whichRoot) const;
  
  // derivative of _powerScalar itself
  Real computeTopScaledPowerGradient() const;
  Real computeBottomScaledPowerGradient() const;
  
  // version with OMP turned off so we can step through each cell
  Real computeTopScaledPowerGradientDebug() const;
  Real computeBottomScaledPowerGradientDebug(const int whichRoot) const;
  Real computeTopScaledPowerGradientDebug(const int whichRoot) const;
 
  Real computeTopBulkPowerGradient() const;
  Real computeBottomBulkPowerGradient() const;

  // compute entries in the Hessian
  Real computeConformalConformalHessian() const;

  Real computeTopConformalHessian(const int whichRoot) const;
  Real computeTopConformalHessianFast(const int whichRoot) const;

  Real computeBottomConformalHessian(const int whichRoot) const;
  Real computeBottomConformalHessianFast(const int whichRoot) const;

  Real computeTopDiagonalHessian(const int whichRoot) const;
  Real computeTopDiagonalHessianFast(const int whichRoot) const;

  Real computeMixedTopHessian(const int firstRoot, const int secondRoot) const;
  Real computeMixedTopHessianFast(const int firstRoot, const int secondRoot) const;

  Real computeBottomDiagonalHessian(const int whichRoot) const;
  Real computeBottomDiagonalHessianFast(const int whichRoot) const;

  Real computeMixedBottomHessian(const int firstRoot, const int secondRoot) const;
  Real computeMixedBottomHessianFast(const int firstRoot, const int secondRoot) const;

  Real computeTopBottomHessian(const int topRoot, const int bottomRoot) const;
  Real computeTopBottomHessianFast(const int topRoot, const int bottomRoot) const;

  // print out the current powers of the roots
  void printRootPowers() const;

  // print out all the root information
  void printRoots() const;

  // reweight the score weight field with curvature
  void reweightWithCurvature(const FIELD_3D& curvatureField);

  // build the whole Hessian
  MATRIX buildHessian();
  MATRIX buildHessianFast();

protected:
  // resolution of the field
  int _xRes;
  int _yRes;
  int _zRes;

  VEC3F _center;
  VEC3F _lengths;

  // the field being drawn and manipulated
  FIELD_3D _fractal;

  // an aux field that stores other iteration information
  FIELD_3D _auxiliary;

  // the polynomials in the rational
  POLYNOMIAL_4D _top;
  POLYNOMIAL_4D _bottom;
  POLYNOMIAL_4DE _topExtended;
  POLYNOMIAL_4DE _bottomExtended;

  // the distance field's center of mass
  VEC3F _centerOfMass;

  // how many iterations to use when computing the fractal?
  int _maxIterations;

  // which quaternion slice to start from?
  Real _quaternionSlice;

  // distance field to compare fractal against
  FIELD_3D _distanceField;
  
  // inverse distance field to compare fractal against
  FIELD_3D _inverseDistanceField;
  Real _inverseFieldSum;
  
  FIELD_3D _weightField;

  // curvature field to compare fractal against
  FIELD_3D _curvatureField;

  // original center and radius of the conformal flow
  //VEC3F _flowCenter;
  //Real _flowRadius;

  // history of polynomial coefficients
  vector<VECTOR> _historyTop;

  int _binaryBandwidth;
  int _curvatureBandwidth;

  // exponent scaling for the log version
  Real _expScaling;

  // the most recently loaded filename
  string _lastLoadedFilename;

  // finite difference dx to use when computing power score
  Real _powerDx;

  // parameters of the tanh function
  Real _tanhAlpha;
  Real _tanhThreshold;

  // cache out the polynomial evaluations
  POLYNOMIAL_CACHE _topCache;
  POLYNOMIAL_CACHE _bottomCache;
  POLYNOMIAL_CACHEE _topCacheExtended;
  POLYNOMIAL_CACHEE _bottomCacheExtended;

  // cache out the log evaluations
  QUATERNION_CACHE _topLogCache;
  QUATERNION_CACHE _bottomLogCache;

  // cache out other quantities shared across Hessian evaluations
  FIELD_3D _RdotR;
  FIELD_3D _conformalConformal;
  FIELD_3D _dSdF;
  QUATERNION_CACHE _topDerivativeCache;
  QUATERNION_CACHE _bottomDerivativeCache;

  // track whether or not a cell overflowed during polynomial computation
  FIELD_3D _overflowed;
  FIELD_3D _nans;
  FIELD_3D _infs;

  // did the adaptive precision machinery fire?
  bool _adaptiveFired;

  // remove the root from the list that is nearest to this position
  void removeNearest(vector<QUATERNION>& roots, const VEC3F& position);

  // try inserting the given root into the list
  void insertNewRoot(vector<QUATERNION>& roots, const QUATERNION& newRoot, int guess = -1);
  
  // try inserting the given root into the list
  void insertBestScoreRoot(vector<QUATERNION>& roots, const QUATERNION& newRoot, const Real radius);

  // compute the score for a sphere with the given radius
  Real computeSphereBinaryScore(const FIELD_3D& fractal, const Real radius);

  // compute the score with respect to a distance field
  Real computeBinaryScore(const FIELD_3D& fractal) const;
  Real computeBinaryBandedScore(const FIELD_3D& fractal) const;
  Real computeInverseDistanceScore(const FIELD_3D& fractal);
  
  // compute the score with respect to a curvature field
  Real computeCurvatureScore(const FIELD_3D& fractal) const;
  Real computeCurvatureBandedScore(const FIELD_3D& fractal) const;

  // compute a score based on the distance field values
  Real computeRootDistanceScore() const;

  // add to a real using Kahan summation
  void kahanAdd(Real& original, const Real& add, Real& error) const;

  // compute the fitting score for a given threshold
  Real inverseDistanceScore(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField, const Real& inverseSum, const Real& threshold);
  Real inverseDistanceScoreDebug(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField, const Real& inverseSum, const Real& threshold);
  Real inverseDistanceScoreConformalGradient(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField, const Real& inverseSum, const Real& threshold);

  Real weightScore(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField);
  Real weightScoreConformalGradient(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField);

  // compute some arbitrary weight field
  void computeWeightField();

  // refresh polynomial evaluation caches
  void refreshPolynomialCache(POLYNOMIAL_CACHE& cache, const POLYNOMIAL_4D& polynomial, const FIELD_3D& fractal);
  void computeLogCache(QUATERNION_CACHE& cache, const POLYNOMIAL_4D& polynomial, const FIELD_3D& fractal, const bool inverse = false);

  // compute some cached quantities for use later in the Hessian
  // assumes that the _topCache and _bottomCache are up to date!!!!
  void computeHessianCaches();
};

#endif
