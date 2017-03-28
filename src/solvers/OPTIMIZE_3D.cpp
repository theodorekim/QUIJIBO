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

#include "OPTIMIZE_3D.h"
//#include <omp.h>
#include <MERSENNETWISTER.h>
#include <algorithm>
#include <float.h>

#include <fenv.h>

#define EXPLICIT_LEADING_ROOT 1
#define OMP_ENABLED 0

#define RECENTER_ITERATE 1

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
OPTIMIZE_3D::OPTIMIZE_3D()
{
  _center = VEC3F(0,0,0);

  //_lengths = VEC3F(1,1,1);
  _lengths = VEC3F(3,3,3);

  //_maxIterations = 100;
  _maxIterations = 1;
  _quaternionSlice = 0;

  //_fractal = FIELD_3D(300,300,300);
  //_fractal = FIELD_3D(200,200,200);
  //_fractal = FIELD_3D(100,100,100);
  //_fractal = FIELD_3D(10,10,10);
  _fractal = FIELD_3D(50,50,50);
  _auxiliary = _fractal;

  _xRes = _fractal.xRes();
  _yRes = _fractal.yRes();
  _zRes = _fractal.zRes();

  //_flowRadius = 1.0;
  //_flowCenter = VEC3F(0,0,0);
  _binaryBandwidth = 16;
  _curvatureBandwidth = 1;
  //_curvatureBandwidth = 4;
  //_binaryBandwidth = 16;
  //_curvatureBandwidth = 4;
 
  //_expScaling = exp(-108);
  _expScaling = 1;

  _powerDx = 0.5;

  _tanhAlpha = 10;
  _tanhThreshold = 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
OPTIMIZE_3D::OPTIMIZE_3D(int xRes, int yRes, int zRes, const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths)
{
  //_maxIterations = 100;
  _maxIterations = 1;
  _quaternionSlice = 0;

  _fractal = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _auxiliary = _fractal;

  _binaryBandwidth = 16;
  _curvatureBandwidth = 1;
 
  //_curvaturePower = 1.0;
  _expScaling = exp(-108);
  _powerDx = 0.5;
  
  _tanhAlpha = 10;
  _tanhThreshold = 0;
}

OPTIMIZE_3D::~OPTIMIZE_3D()
{
}

///////////////////////////////////////////////////////////////////////
// get all the roots, all at once
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::getRoots() const
{
  int totalRoots = _top.totalRoots();
  int totalDOFs = 4 * totalRoots;

  VECTOR final(totalDOFs);

  for (int x = 0; x < _top.totalRoots(); x++)
  {
    final[4 * x] = _top.roots()[x][0];
    final[4 * x + 1] = _top.roots()[x][1];
    final[4 * x + 2] = _top.roots()[x][2];
    final[4 * x + 3] = _top.roots()[x][3];
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::getRationalRoots() const
{
  int totalRoots = _top.totalRoots() + _bottom.totalRoots();
  int totalDOFs = 4 * totalRoots;

  VECTOR final(totalDOFs);

  for (int x = 0; x < _top.totalRoots(); x++)
  {
    final[4 * x] = _top.roots()[x][0];
    final[4 * x + 1] = _top.roots()[x][1];
    final[4 * x + 2] = _top.roots()[x][2];
    final[4 * x + 3] = _top.roots()[x][3];
  }
  int index = 4 * _top.totalRoots();
  for (int x = 0; x < _bottom.totalRoots(); x++, index += 4)
  {
    final[index] = _bottom.roots()[x][0];
    final[index + 1] = _bottom.roots()[x][1];
    final[index + 2] = _bottom.roots()[x][2];
    final[index + 3] = _bottom.roots()[x][3];
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::getBottomRoots() const
{
  int totalRoots = _bottom.totalRoots();
  int totalDOFs = 4 * totalRoots;

  VECTOR final(totalDOFs);

  int index = 0; 
  for (int x = 0; x < _bottom.totalRoots(); x++, index += 4)
  {
    final[index] = _bottom.roots()[x][0];
    final[index + 1] = _bottom.roots()[x][1];
    final[index + 2] = _bottom.roots()[x][2];
    final[index + 3] = _bottom.roots()[x][3];
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// set all the roots, all at once
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::setRoots(const VECTOR& newRoots)
{
  int totalRoots = _top.totalRoots();
  int totalDOFs = 4 * totalRoots;

  assert(newRoots.size() == totalDOFs);

  for (int x = 0; x < _top.totalRoots(); x++)
  {
    QUATERNION root;
    root[0] = newRoots[4 * x];
    root[1] = newRoots[4 * x + 1];
    root[2] = newRoots[4 * x + 2];
    root[3] = newRoots[4 * x + 3];
    _top.modifyRoot(x, root);
  }
}

///////////////////////////////////////////////////////////////////////
// set all the rational roots, all at once
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::setRationalRoots(const VECTOR& newRoots)
{
  int totalRoots = _top.totalRoots() + _bottom.totalRoots();
  int totalDOFs = 4 * totalRoots;

  assert(newRoots.size() == totalDOFs);

  for (int x = 0; x < _top.totalRoots(); x++)
  {
    QUATERNION root;
    root[0] = newRoots[4 * x];
    root[1] = newRoots[4 * x + 1];
    root[2] = newRoots[4 * x + 2];
    root[3] = newRoots[4 * x + 3];
    _top.modifyRoot(x, root);
  }
  int index = 4 * _top.totalRoots();
  for (int x = 0; x < _bottom.totalRoots(); x++, index += 4)
  {
    QUATERNION root;
    root[0] = newRoots[index];
    root[1] = newRoots[index + 1];
    root[2] = newRoots[index + 2];
    root[3] = newRoots[index + 3];
    _bottom.modifyRoot(x, root);
  }
}

///////////////////////////////////////////////////////////////////////
// set all the rational roots, all at once
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::setBottomRoots(const VECTOR& newRoots)
{
  int totalRoots = _bottom.totalRoots();
  int totalDOFs = 4 * totalRoots;

  assert(newRoots.size() == totalDOFs);

  int index = 0;
  for (int x = 0; x < _bottom.totalRoots(); x++, index += 4)
  {
    QUATERNION root;
    root[0] = newRoots[index];
    root[1] = newRoots[index + 1];
    root[2] = newRoots[index + 2];
    root[3] = newRoots[index + 3];
    _bottom.modifyRoot(x, root);
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::initCube()
{
  vector<QUATERNION> roots;

  Real last = 0;
  roots.push_back(QUATERNION(-0.5,-0.5,-0.5, last));
  //roots.push_back(QUATERNION(-0.3,-0.3,-0.4, last));
  roots.push_back(QUATERNION(-0.5,-0.5,0.5, last));
  roots.push_back(QUATERNION(-0.5,0.5,-0.5, last));
  roots.push_back(QUATERNION(-0.5,0.5,0.5, last));
  roots.push_back(QUATERNION(0.5,-0.5,-0.5, last));
  roots.push_back(QUATERNION(0.5,-0.5,0.5, last));
  roots.push_back(QUATERNION(0.5,0.5,-0.5, last));
  roots.push_back(QUATERNION(0.5,0.5,0.5, last));

  _top = POLYNOMIAL_4D(roots);

  for (unsigned int x = 0; x < roots.size(); x++)
    roots[x] *= 0.25;

  roots.resize(roots.size() - 3);

  //_bottom = POLYNOMIAL_4D(roots);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::initRandomSphere()
{
  MERSENNETWISTER twister(123456);

  int totalTopRoots = 50;
  //int totalBottomRoots = 0;

  vector<QUATERNION> topRoots;
  for (int x = 0; x < totalTopRoots; x++)
  {
    QUATERNION newRoot(twister.rand() - 0.5, twister.rand() - 0.5, twister.rand() - 0.5, 0);
    newRoot.normalize();
    //newRoot[3] = twister.rand() - 0.5;

    //newRoot += QUATERNION(0.5, 0.5, 0.5, 0.5);
    //newRoot += QUATERNION(1.5, 1.5, 1.5, 1.5);
    topRoots.push_back(newRoot);
  }

  for (int x = 0; x < totalTopRoots; x++)
  {
    topRoots[x] *= 2;
    topRoots[x] *= 0.2;
  }
  _top = POLYNOMIAL_4D(topRoots);

/*
  vector<QUATERNION> bottomRoots;
  for (int x = 0; x < totalBottomRoots; x++)
  {
    QUATERNION newRoot(2 * twister.rand() - 0.5, twister.rand() - 0.5, twister.rand() - 0.5, 0);
    //newRoot.normalize();
    newRoot[3] = twister.rand() - 0.5;
    //newRoot += QUATERNION(1.5, 1.5, 1.5, 1.5);
    bottomRoots.push_back(newRoot);
  }
  for (unsigned int x = 0; x < totalBottomRoots; x++)
    bottomRoots[x] *= 0.25;

  _bottom = POLYNOMIAL_4D(bottomRoots);
  */
}

//////////////////////////////////////////////////////////////////////
// try inserting the given root into the list
//////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::insertNewRoot(vector<QUATERNION>& roots, const QUATERNION& newRoot, int guess)
{
  vector<QUATERNION> bestSeen;
  int bestIndex = -1;
  Real bestSeenConditionNumber = 0;
  cout << " original roots size: " << roots.size() << endl;

  // try inserting the new root in every possible position
  //for (int x = 0; x < roots.size() + 1; x++)
  unsigned int x = 11;
  {
    if (guess > 0)
      x = guess;

    vector<QUATERNION> testRoots;

    // insert all the roots
    for (unsigned int y = 0; y < roots.size() + 1; y++)
    {
      if (x == y)
        testRoots.push_back(newRoot);

      if (y < roots.size())
        testRoots.push_back(roots[y]);
    }

    cout << " Testing " << testRoots.size() << " roots " << endl;
    POLYNOMIAL_4D testPolynomial(testRoots);
    Real conditionNumber = testPolynomial.conditionNumber();

    if (conditionNumber > bestSeenConditionNumber || x == 0)
    {
      bestSeen = testRoots;
      bestIndex = x;
      bestSeenConditionNumber = conditionNumber;
    }

    cout << " Condition number " << x << ": " << conditionNumber << endl;

    if (guess > 0)
      x = roots.size() + 1;
  }
  roots = bestSeen;

  cout << " Best seen condition number " << bestIndex << ": " << bestSeenConditionNumber << endl; 
  //exit(0);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::initPolynomial()
{
  /*
  _top.addRoot(QUATERNION(0,0,0,-1));
  _top.addRoot(QUATERNION(0,0,0,1));

  _bottom.addRoot(QUATERNION(0,0,0.1,0.1));
  _bottom.addRoot(QUATERNION(0,0,0,-0.1));
  */

  // try one where everybody is in the complex plane
  _top.addRoot(QUATERNION(-1, 0,0,0));
  _top.addRoot(QUATERNION(1, 0,0,0));
  
  _bottom.addRoot(QUATERNION(0.1,0.1,0,0));
  _bottom.addRoot(QUATERNION(-2,0,0,0));
}

//////////////////////////////////////////////////////////////////////
// Print state to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, const OPTIMIZE_3D& optimize2D)
{
  out << " Top roots: " << endl;
  out << optimize2D.top() << endl;
  out << " Bottom roots: " << endl;
  out << optimize2D.bottom() << endl;

  return out;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computePolynomial(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing rational 3D Newton ..." << flush;

  //const int maxIterations = 100;
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();
  //cout << " res: " << xRes << " " << yRes << " " << zRes << " " << flush;
  //cout << " lengths: " << _lengths << endl;

  const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif 

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
  {
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif
        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        
        Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        Real pPrimeMagnitude = 1;
        int totalIterations = 0;
        
        QUATERNION p, pPrime;
        //bool hasNan = isnan(magnitude);

        /*
        bool debug = (x == 36 && y == 25 && z == 25);
        if (debug)
          cout << " Original iterate: " << iterate << endl;
          */

        //while (magnitude > 1e-5 && !hasNan && totalIterations < maxIterations)
        while (magnitude > eps && totalIterations < _maxIterations && pPrimeMagnitude > eps)
        {
          POLYNOMIAL_4D::evaluateRational(top, bottom, iterate, p, pPrime);

          QUATERNION div(p / pPrime);
          iterate -= div;

          pPrimeMagnitude = pPrime.magnitude();
          distanceTravelled += div.magnitude();
          //hasNan = isnan(magnitude);
          //hasNan = magnitude != magnitude;
          magnitude = p.magnitude();
          totalIterations++;

          /*
          if (debug)
          {
            cout << "================================" << endl;
            cout << " Iterate " << totalIterations << endl;
            cout << "================================" << endl;
            cout << " p:       " << p << endl;
            cout << " pPrime:  " << pPrime << endl;
            cout << " iterate: " << iterate << endl;
          }
          */
        }
        /*
        if (debug)  {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " DEBUGGING " << endl;
          exit(0);
        }
        */

        fractal(x,y,z) = (totalIterations == _maxIterations || 
                          magnitude > eps || 
                          pPrimeMagnitude < eps) ? -1.0 : 1.0;

        //_auxiliary(x,y,z) = distanceTravelled;
        _auxiliary(x,y,z) = totalIterations;
        //fractal(x,y,z) = magnitude;
        //fractal(x,y,z) = totalIterations;
        //fractal(x,y,z) = pPrime.anyNans();
        
        //fractal(x,y,z) = iterate.anyNans() ? 0 : totalIterations;
        //fractal(x,y,z) = iterate.anyNans() ? 0 : distanceTravelled;
      }
    if (z % (int)(zRes / 10) == 0)
      cout << 100 * ((Real)z / zRes) << "% " << flush;
  }
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// list of roots within this z slice
///////////////////////////////////////////////////////////////////////
vector<int> OPTIMIZE_3D::visibleRootsZ(const int zSlice)
{
  vector<int> final;

  VEC3F corner = _center - (Real)0.5 * _lengths;
  Real dz = _lengths[2] / _zRes;
  VEC3F dxs(_lengths[0] / _xRes, _lengths[1] / _yRes, _lengths[2] / _zRes);
  corner += (Real)0.5 * dxs;

  for (int x = 0; x < _top.totalRoots(); x++)
  {
    QUATERNION center = _top.roots()[x];
    VEC3F center3(center.w(), center.x(), center.y());
    
    
    center3 -= corner;
    center3[2] *= 1.0 / dz;
    int z0 = (int)center3[2] + 1;
    z0 = (z0 < 0) ? 0 : z0;
    z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

    if (z0 != zSlice) continue;

    final.push_back(x);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// list of roots within this z slice
///////////////////////////////////////////////////////////////////////
vector<int> OPTIMIZE_3D::visibleRootsY(const int ySlice)
{
  vector<int> final;

  VEC3F corner = _center - (Real)0.5 * _lengths;
  Real dy = _lengths[1] / _yRes;
  VEC3F dxs(_lengths[0] / _xRes, _lengths[1] / _yRes, _lengths[2] / _zRes);
  corner += (Real)0.5 * dxs;

  for (int x = 0; x < _top.totalRoots(); x++)
  {
    QUATERNION center = _top.roots()[x];
    VEC3F center3(center.w(), center.x(), center.y());
    
    center3 -= corner;
    center3[1] *= 1.0 / dy;
    int y0 = (int)center3[1] + 1;
    y0 = (y0 < 0) ? 0 : y0;
    y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;

    if (y0 != ySlice) continue;

    final.push_back(x);
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// init the roots of the rational
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::initializeDebugCircleRoots(const int totalTopRoots, const int totalBottomRoots)
{
  assert(totalTopRoots <= totalBottomRoots);

  Real radius = 0.5;

  vector<QUATERNION> topRoots;
  vector<QUATERNION> bottomRoots;

  Real topRadius = 0.9;
  QUATERNION center(0,0);

  for (int i = 0; i < totalTopRoots; i++)
  {
    Real fraction = 2.0 * M_PI * (Real)i / totalTopRoots;
    Real xReal = sin(fraction);
    Real yReal = cos(fraction);

    //QUATERNION point(xReal, yReal, 0,0);
    //QUATERNION point(0, xReal, yReal, 0);
    //QUATERNION point(xReal, 0, yReal, 0);
    QUATERNION point(xReal, 0,0,  yReal);
    point = topRadius * radius * point + center;
    topRoots.push_back(point);
  }
  for (int i = 0; i < totalBottomRoots; i++)
  {
    Real fraction = 2.0 * M_PI * (Real)i / totalBottomRoots;
    Real xReal = sin(fraction);
    Real yReal = cos(fraction);

    //QUATERNION point(xReal, yReal,0,0);
    //QUATERNION point(0, xReal, yReal,0);
    //QUATERNION point(xReal, 0, yReal,0);
    QUATERNION point(xReal, 0, 0, yReal);
    point = radius * point + center;
    bottomRoots.push_back(point);
  }
  _top = POLYNOMIAL_4D(topRoots);
  _bottom = POLYNOMIAL_4D(bottomRoots);
}

///////////////////////////////////////////////////////////////////////
// remove the root from the list that is nearest to this position
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::removeNearest(vector<QUATERNION>& roots, const VEC3F& position)
{
  Real nearestMagnitude = 0;
  vector<QUATERNION>::iterator iter;
  vector<QUATERNION>::iterator nearest;
  //int nearestIndex = 0;

  iter = roots.begin();
  for (unsigned int x = 0; x < roots.size(); x++, iter++)
  {
    VEC3F root(roots[x][0], roots[x][1], roots[x][2]);
    VEC3F diff = position - root;

    Real magnitude = diff.magnitude();
    if (magnitude < nearestMagnitude || x == 0)
    {
      nearest = iter;
      nearestMagnitude = magnitude;
      //nearestIndex = x;
    }
  }

  cout << " Removing nearest to " << position << ", " << *nearest << endl;

  roots.erase(nearest);
}

///////////////////////////////////////////////////////////////////////
// try inserting the given root into the list
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::insertBestScoreRoot(vector<QUATERNION>& roots, const QUATERNION& newRoot, const Real radius)
{
  const vector<QUATERNION> originalRoots = roots;

  Real bestScore = 1;
  int bestPosition = 0;
  
  // try inserting the new root in every possible position
  for (unsigned int x = 0; x < originalRoots.size() + 1; x++)
  {
    vector<QUATERNION> testRoots;

    // insert all the roots
    for (unsigned int y = 0; y < originalRoots.size() + 1; y++)
    {
      if (x == y)
        testRoots.push_back(newRoot);

      if (y < originalRoots.size())
        testRoots.push_back(originalRoots[y]);
    }

    POLYNOMIAL_4D bottomTest(testRoots);
  
    for (unsigned int y = 0; y < testRoots.size(); y++)
      testRoots[y] *= 0.5;
    testRoots.resize(testRoots.size() - 4);

    POLYNOMIAL_4D topTest(testRoots);

    FIELD_3D fractal(50, 50, 50, VEC3F(0,0,0), VEC3F(1,1,1));
    computePolynomial(topTest, bottomTest, fractal);

    Real score = computeSphereBinaryScore(fractal, radius);

    cout << " Score: " << score << endl;

    if (score < bestScore || x == 0)
    {
      bestScore = score;
      bestPosition = x;
      roots = bottomTest.roots();
    }
  }
  cout << " Best score seen: " << bestScore << " in position: " << bestPosition << endl;
}

///////////////////////////////////////////////////////////////////////
// score based on a binary fractal and weight function
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeSphereBinaryScore(const FIELD_3D& fractal, const Real radius)
{
  TIMER functionTimer(__FUNCTION__);
  // compute the score
  Real score = 0;

  int xRes = fractal.xRes();
  int yRes = fractal.yRes();
  int zRes = fractal.zRes();

  FIELD_3D product(xRes, yRes, zRes, VEC3F(0,0,0), VEC3F(1,1,1));
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        VEC3F center = fractal.cellCenter(x,y,z);

        Real magnitude = center.magnitude();

        double thresholded = (magnitude < radius) ? -1 : 1; 
        product(x,y,z) = fractal(x,y,z) * thresholded;
        //product(x,y,z) = thresholded;
        score += fractal(x,y,z) * thresholded;
      }

  /*
  FIELDVIEW3D(product);
  FIELDVIEW3D(fractal);

  exit(0);
  */

  return score / fractal.totalCells();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeFactoredPolynomial()
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing factored rational 3D Newton ..." << flush;

  POLYNOMIAL_4D& top = _top;
  POLYNOMIAL_4D& bottom = _bottom;
  FIELD_3D& fractal = _fractal;
  
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5; 
#endif

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
  {
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = ((Real)x / xRes * _lengths[0] + xMin);
        Real yReal = ((Real)y / yRes * _lengths[1] + yMin);
        Real zReal = ((Real)z / zRes * _lengths[2] + zMin);
#endif
        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        
        Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        Real pPrimeMagnitude = 1;
        int totalIterations = 0;
        
        QUATERNION p, pPrime;
        //bool hasNan = isnan(magnitude);

        /*
        bool debug = (x == 36 && y == 25 && z == 25);
        if (debug)
        {
          cout << " Original iterate: " << iterate << endl;
          cout << " Top:    " << top << endl;
          cout << " Bottom: " << bottom << endl;
        }
        */

        //while (magnitude > 1e-5 && !hasNan && totalIterations < maxIterations)
        while (magnitude > eps && totalIterations < _maxIterations && pPrimeMagnitude > eps)
        {
          //POLYNOMIAL_4D::evaluateRational(top, bottom, iterate, p, pPrime);
          POLYNOMIAL_4D::evaluateFactoredRational(top, bottom, iterate, p, pPrime);
          //POLYNOMIAL_4D::evaluateFactoredQuadratic(top, bottom, iterate, p, pPrime);

          QUATERNION div(p / pPrime);
          iterate -= div;

          pPrimeMagnitude = pPrime.magnitude();
          distanceTravelled += div.magnitude();
          magnitude = p.magnitude();
          totalIterations++;
          /*
          if (debug)
          {
            cout << "================================" << endl;
            cout << " Iterate " << totalIterations << endl;
            cout << "================================" << endl;
            cout << " p:       " << p << endl;
            cout << " pPrime:  " << pPrime << endl;
            cout << " iterate: " << iterate << endl;
          }
          */
        }
        /*
        if (debug)  
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " DEBUGGING " << endl;
          exit(0);
        }
        */
        fractal(x,y,z) = (totalIterations == _maxIterations || 
                          magnitude > eps || 
                          pPrimeMagnitude < eps) ? -1.0 : 1.0;

        //_auxiliary(x,y,z) = distanceTravelled;
        _auxiliary(x,y,z) = totalIterations;
        //fractal(x,y,z) = magnitude;
        //fractal(x,y,z) = totalIterations;
        //fractal(x,y,z) = pPrime.anyNans();
        
        //fractal(x,y,z) = iterate.anyNans() ? 0 : totalIterations;
        //fractal(x,y,z) = iterate.anyNans() ? 0 : distanceTravelled;
      }
    if (z % (int)(zRes / 10) == 0)
      cout << 100 * ((Real)z / zRes) << "% " << flush;
  }
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeMap(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, FIELD_3D& auxiliary, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  //_lengths = VEC3F(6,6,6);
  //_lengths = VEC3F(3,3,3);

  //POLYNOMIAL_4D& top = _top;
  //POLYNOMIAL_4D& bottom = _bottom;
  //FIELD_3D& fractal = _fractal;
  //FIELD_3D& auxiliary = _auxiliary;
  
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  //Real xMin = _center[0] - _lengths[0] * 0.5;
  //Real yMin = _center[1] - _lengths[1] * 0.5; 
  //Real zMin = _center[2] - _lengths[2] * 0.5;

  //Real escape = 2.0;
  //Real escape = 4.0;
  Real escape = 10.0;
  //Real escape = 100.0;

  if (verbose)
    cout << " Computing 3D Newton map ..." << flush;
#define FACTORED 1

  if (verbose)
#if FACTORED
  cout << " Using factored form " << flush;
#else
  cout << " Using simple form " << flush;
#endif

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
  {
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        //QUATERNION iterate(0, xReal, yReal, zReal);
        
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;
 
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
        /*
          if (_bottom.totalRoots() > 0 && _top.totalRoots() > 0)
          {
#if FACTORED
            QUATERNION top = _top.evaluateFactored(iterate);
            QUATERNION bottom = _bottom.evaluateFactored(iterate);
#else
            QUATERNION top = _top.evaluate(iterate);
            QUATERNION bottom = _bottom.evaluate(iterate);
#endif
            iterate = top / bottom;
          }
          else if (_top.totalRoots() > 0)
          {
#if FACTORED
            QUATERNION top = _top.evaluateFactored(iterate);
#else
            QUATERNION top = _top.evaluate(iterate);
#endif
            iterate = top;
          }
          else
          {
#if FACTORED
            QUATERNION bottom = _bottom.evaluateFactored(iterate);
#else
            QUATERNION bottom = _bottom.evaluate(iterate);
#endif
            iterate = bottom;
          }
          */
          iterate = polynomial.evaluateFactored(iterate);
          magnitude = iterate.magnitude();
          totalIterations++;
        }
        fractal(x,y,z) = (totalIterations == _maxIterations || magnitude < escape) ? -1.0 : 1.0;
        //fractal(x,y,z) = iterate.x();
        //fractal(x,y,z) = magnitude;
        //fractal(x,y,z) = escape;

        //_auxiliary(x,y,z) = distanceTravelled;
        auxiliary(x,y,z) = totalIterations;
        //fractal(x,y,z) = magnitude;
        //fractal(x,y,z) = totalIterations;
        //fractal(x,y,z) = pPrime.anyNans();
        
        //fractal(x,y,z) = iterate.anyNans() ? 0 : totalIterations;
        //fractal(x,y,z) = iterate.anyNans() ? 0 : distanceTravelled;
      }

    if (verbose)
      if (z % (int)(zRes / 10) == 0)
        cout << 100 * ((Real)z / zRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogBandedMap(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();
  const int slabSize = xRes * yRes;

  int maxRes = (xRes > yRes) ? xRes : yRes;
  maxRes = (maxRes > zRes) ? maxRes : zRes;

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  //int maxIterations = 20;
  //int maxIterations = 1;
  //int maxIterations = 2;
  if (verbose)
    cout << " Computing Julia set ... " << flush;

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;

        double sign = (_distanceField[index] < 0) ? 1 : -1; 
        int cellsInside = fabs(_distanceField[index] * maxRes);

        if (sign < 0 && cellsInside > _binaryBandwidth)
          continue;
        if (sign > 0 && cellsInside > _binaryBandwidth)
          continue;

#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;
 
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          // iterate = polynomial.evaluateFactored(iterate);
          //iterate = scaling * iterate * _top.evaluateFactored(iterate);
          iterate = _expScaling * iterate * polynomial.evaluateFactored(iterate);
          magnitude = iterate.magnitude();
          totalIterations++;
        }
        //fractal(x,y,z) = magnitude;

        fractal(x,y,z) = log(magnitude);

/*
        if (fabs(magnitude) > 0)
          fractal(x,y,z) = log(magnitude);
        else
          fractal(x,y,z) = 0;
          */

        //if (isinf(fractal(x,y,z)))  
        //  fractal(x,y,z) = 0;

        //fractal(x,y,z) = (totalIterations == _maxIterations || magnitude < escape) ? -1.0 : 1.0;
        //fractal(x,y,z) = totalIterations;
      }

    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;

  //FIELDVIEW3D(fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogMap(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  //int maxIterations = 20;
  //int maxIterations = 1;
  //int maxIterations = 2;
  if (verbose)
    cout << " Computing Julia set ... " << flush;

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;
 
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          iterate = _expScaling * iterate * polynomial.evaluateFactored(iterate);
          magnitude = iterate.magnitude();
          totalIterations++;
        }
        //fractal(x,y,z) = magnitude;

        fractal(x,y,z) = log(magnitude);

/*
        if (fabs(magnitude) > 0)
          fractal(x,y,z) = log(magnitude);
        else
          fractal(x,y,z) = 0;
          */

        //if (isinf(fractal(x,y,z)))  
        //  fractal(x,y,z) = 0;

        //fractal(x,y,z) = (totalIterations == _maxIterations || magnitude < escape) ? -1.0 : 1.0;
        //fractal(x,y,z) = totalIterations;
      }

    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;

  //FIELDVIEW3D(fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogPowerRationalMapCached(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();
  //FIELDVIEW3D(_weightField);
  //exit(0);

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  if (verbose)
    cout << " Computing powerized rational Julia set ... " << flush;

  _topCache = POLYNOMIAL_CACHE(xRes, yRes, zRes);
  _bottomCache = POLYNOMIAL_CACHE(xRes, yRes, zRes);

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        //Real magnitudeTrial = 0;
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          vector<QUATERNION>& forwardCacheTop = _topCache.forward(x,y,z);
          vector<QUATERNION>& backwardCacheTop = _topCache.backward(x,y,z);
#if EXPLICIT_LEADING_ROOT
          QUATERNION topEval = top.evaluatePowerFactored(iterate, forwardCacheTop, backwardCacheTop);
#else
          QUATERNION topEval = iterate * top.evaluatePowerFactored(iterate, forwardCacheTop, backwardCacheTop);
#endif

/*
          if (x == 80 && y == 54 && z == 48)
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << " cache size: " << forwardCacheTop.size() << endl;
          }
          */

          QUATERNION bottomEval;
          
          if (bottom.totalRoots() > 0)
          {
            vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(x,y,z);
            vector<QUATERNION>& backwardCacheBottom = _bottomCache.backward(x,y,z);
#if EXPLICIT_LEADING_ROOT
            bottomEval = bottom.evaluatePowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
#else
            bottomEval = iterate * bottom.evaluatePowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
#endif
            iterate = (topEval / bottomEval);
          }
          else
            iterate = topEval;

          //magnitudeTrial = iterate.magnitude();

          iterate *= _expScaling;
          magnitude = iterate.magnitude();

          totalIterations++;
        }

        //Real logKnown = log(magnitude);
        //Real logTrial = log(magnitudeTrial) + log(_expScaling);

        fractal(x,y,z) = log(magnitude);
      }
    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogScaledPowerRationalMap(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //FIELDVIEW3D(_weightField);
  //exit(0);

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  if (verbose)
    cout << " Computing powerized rational Julia set ... " << flush;

  _overflowed = FIELD_3D(xRes, yRes, zRes);
  _overflowed = (Real)0.0;

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        /*
        vector<QUATERNION> iterates;
        vector<QUATERNION> tops;
        vector<QUATERNION> bottoms;
        vector<Real> magnitudes;
        */
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          // DEBUG
          //feclearexcept(FE_OVERFLOW);

          const QUATERNION input = iterate;

          QUATERNION topEval = top.evaluateScaledPowerFactored(iterate);
          QUATERNION bottomEval;
          
          if (bottom.totalRoots() > 0)
          {
            bottomEval = bottom.evaluateScaledPowerFactored(iterate);
            iterate = (topEval / bottomEval);
          }
          else
            iterate = topEval;

          // DEBUG
          //if (fetestexcept(FE_OVERFLOW))
          //  _overflowed(x,y,z) = (Real)1.0;

          //magnitudeTrial = iterate.magnitude();

          iterate *= _expScaling;
          magnitude = iterate.magnitude();

          /*
          iterates.push_back(iterate);
          magnitudes.push_back(magnitude);
          tops.push_back(topEval);
          bottoms.push_back(bottomEval);

          // DEBUG
          if (isinf(log(magnitude)))
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << " Found inf: " << log(magnitude) << std::endl;
            cout << " Threshold: " << 10.0 * REAL_MIN << endl;
            cout << " single: " << FLT_MIN << endl;
            cout << " double: " << DBL_MIN << endl;
            cout << " extended: " << LDBL_MIN << endl;
            cout << " total iterations: " << totalIterations << endl;
            cout << " (x,y,z): " << x << "," << y << "," << z << endl;
            const QUATERNION original(xReal, yReal, zReal, _quaternionSlice);
            cout << " original: " << original << endl;
            cout << " iterates: " << endl;
            for (int i = 0; i < iterates.size(); i++)
              cout << iterates[i] << endl;
            cout << " magnitudes: " << endl;
            for (int i = 0; i < magnitudes.size(); i++)
              cout << magnitudes[i] << endl;
            cout << " tops: " << endl;
            for (int i = 0; i < tops.size(); i++)
              cout << tops[i] << endl;
            cout << " bottoms: " << endl;
            for (int i = 0; i < bottoms.size(); i++)
              cout << bottoms[i] << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            //exit(0);
          }
          if (isnan(magnitude))
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << " Found nan: " << magnitude << std::endl;
            cout << " Threshold: " << 10.0 * REAL_MIN << endl;
            cout << " single: " << FLT_MIN << endl;
            cout << " double: " << DBL_MIN << endl;
            cout << " extended: " << LDBL_MIN << endl;
            cout << " total iterations: " << totalIterations << endl;
            cout << " (x,y,z): " << x << "," << y << "," << z << endl;
            const QUATERNION original(xReal, yReal, zReal, _quaternionSlice);
            cout << " input: " << input << endl;

            cout << " original: " << original << endl;
            cout << endl;
            cout << " iterates: " << endl;
            for (int i = 0; i < iterates.size(); i++)
              cout << iterates[i] << endl;
            cout << endl;
            cout << " magnitudes: " << endl;
            for (int i = 0; i < magnitudes.size(); i++)
              cout << magnitudes[i] << endl;
            cout << endl;
            cout << " tops: " << endl;
            for (int i = 0; i < tops.size(); i++)
              cout << tops[i] << endl;
            cout << endl;
            cout << " bottoms: " << endl;
            for (int i = 0; i < bottoms.size(); i++)
              cout << bottoms[i] << endl;
            cout << endl;
         
            cout << " TOP: " << endl; 
            QUATERNION redo = top.evaluateScaledPowerFactoredVerbose(input);
            cout << endl;
            cout << " BOTTOM: " << endl; 
            QUATERNION redoBottom = bottom.evaluateScaledPowerFactoredVerbose(input);
            cout << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            exit(0);
          }
          */

          totalIterations++;

          // see if it fell into the black hole at the origin
          if (magnitude < 10.0 * REAL_MIN)
            totalIterations = _maxIterations;
        }

        //Real logKnown = log(magnitude);
        //Real logTrial = log(magnitudeTrial) + log(_expScaling);

        fractal(x,y,z) = log(magnitude);
      }
    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeRationalMapAP(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  if (verbose)
    cout << " Computing powerized rational Julia set ... " << flush;

  _topCache = POLYNOMIAL_CACHE(xRes, yRes, zRes);
  _bottomCache = POLYNOMIAL_CACHE(xRes, yRes, zRes);

  _overflowed = FIELD_3D(xRes, yRes, zRes);
  _overflowed = 0;
  _nans = _overflowed;
  _infs = _overflowed;

  _topCacheExtended = POLYNOMIAL_CACHEE(xRes, yRes, zRes);
  _bottomCacheExtended = POLYNOMIAL_CACHEE(xRes, yRes, zRes);

  _topExtended = POLYNOMIAL_4DE(top);
  _bottomExtended = POLYNOMIAL_4DE(bottom);

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
      
        Real magnitude = iterate.magnitude();
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        int totalIterations = 0;

        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          vector<QUATERNION>& forwardCacheTop = _topCache.forward(x,y,z);
          vector<QUATERNION>& backwardCacheTop = _topCache.backward(x,y,z);

          // DEBUG
          feclearexcept(FE_OVERFLOW);

#if EXPLICIT_LEADING_ROOT
          QUATERNION topEval = top.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
#else
          QUATERNION topEval = iterate * top.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
#endif
          QUATERNION bottomEval;
          
          if (bottom.totalRoots() > 0)
          {
            vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(x,y,z);
            vector<QUATERNION>& backwardCacheBottom = _bottomCache.backward(x,y,z);
#if EXPLICIT_LEADING_ROOT
            bottomEval = bottom.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
#else
            bottomEval = iterate * bottom.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
#endif
            iterate = (topEval / bottomEval);

          }
          else
            iterate = topEval;

          if (fetestexcept(FE_OVERFLOW))
            _overflowed(x,y,z) = 1;

          iterate *= _expScaling;
          magnitude = iterate.magnitude();

          if (isinf(magnitude) || isnan(magnitude))
            _overflowed(x,y,z) = 1;

          totalIterations++;
        }

        fractal(x,y,z) = log(magnitude);
      }
    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;

  // if there were no overflows, we're all done.
  int overflows = _overflowed.sum();
  if (overflows == 0) return;
  _adaptiveFired = true;

  // second adaptive precision pass
  cout << " Recomputing " << overflows << " overflows " << flush;

  FIELD_3D stillNans(_overflowed);
  stillNans = 0;

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        // only recompute overflowed cells
        if (_overflowed(x,y,z) < 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        long double xReal = point[0];
        long double yReal = point[1];
        long double zReal = point[2];
#else
        long double xReal = (long double)x / xRes * _lengths[0] + xMin + _center[0];
        long double yReal = (long double)y / yRes * _lengths[1] + yMin + _center[1];
        long double zReal = (long double)z / zRes * _lengths[2] + zMin + _center[2];
#endif
        QUATERNIONE iterate(xReal, yReal, zReal, _quaternionSlice);
      
        long double magnitude = iterate.magnitude();
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        int totalIterations = 0;

        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          vector<QUATERNIONE>& forwardCacheTop = _topCacheExtended.forward(x,y,z);
          vector<QUATERNIONE>& backwardCacheTop = _topCacheExtended.backward(x,y,z);

          //feclearexcept(FE_OVERFLOW);

          QUATERNIONE topEval = _topExtended.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
          QUATERNIONE bottomEval;
          
          if (bottom.totalRoots() > 0)
          {
            vector<QUATERNIONE>& forwardCacheBottom = _bottomCacheExtended.forward(x,y,z);
            vector<QUATERNIONE>& backwardCacheBottom = _bottomCacheExtended.backward(x,y,z);
            bottomEval = _bottomExtended.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
            iterate = (topEval / bottomEval);
          }
          else
            iterate = topEval;

          //if (fetestexcept(FE_OVERFLOW))
          //  _overflowed(x,y,z) = 1;

          iterate *= _expScaling;
          magnitude = iterate.magnitude();

          if (isinf(magnitude) || isnan(magnitude))
            _overflowed(x,y,z) = 1;

          totalIterations++;
        }

        fractal(x,y,z) = log(magnitude);

        if (isnan(fractal(x,y,z)))
        {
          fractal(x,y,z) = 0;
          stillNans(x,y,z) = 1;
        }
      }
  }

  cout << "\t" << stillNans.totalNans() << " NaNs persist after adaptive pass " << flush;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogScaledPowerRationalMapCached(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();
  //FIELDVIEW3D(_weightField);
  //exit(0);

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  if (verbose)
    cout << " Computing powerized rational Julia set ... " << flush;

  _topCache = POLYNOMIAL_CACHE(xRes, yRes, zRes);
  _bottomCache = POLYNOMIAL_CACHE(xRes, yRes, zRes);

  _overflowed = FIELD_3D(xRes, yRes, zRes);
  _overflowed = (Real)0.0;

  //FIELD_3D peekTop = _overflowed;
  //FIELD_3D peekBottom = _overflowed;

  // compute a ground truth solution
#if OMP_ENABLED
#pragma omp parallel
#pragma omp for schedule(dynamic)
#endif
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
      
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        int totalIterations = 0;

        //Real magnitudeTrial = 0;
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          vector<QUATERNION>& forwardCacheTop = _topCache.forward(x,y,z);
          vector<QUATERNION>& backwardCacheTop = _topCache.backward(x,y,z);

          // DEBUG
          feclearexcept(FE_OVERFLOW);

#if EXPLICIT_LEADING_ROOT
          QUATERNION topEval = top.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
#else
          QUATERNION topEval = iterate * top.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
#endif

          QUATERNION bottomEval;
          
          if (bottom.totalRoots() > 0)
          {
            vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(x,y,z);
            vector<QUATERNION>& backwardCacheBottom = _bottomCache.backward(x,y,z);
#if EXPLICIT_LEADING_ROOT
            bottomEval = bottom.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
#else
            bottomEval = iterate * bottom.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
#endif
            iterate = (topEval / bottomEval);

          }
          else
            iterate = topEval;

          // DEBUG
          if (fetestexcept(FE_OVERFLOW))
          {
          /*
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << " OVERFLOW DETECTED: " << endl;

            cout << " Top multiply series: " << flush;
            for (int d = 0; d < forwardCacheTop.size(); d++)
              cout << forwardCacheTop[d] << " " << flush;
            cout << endl;
            exit(0);
            */
            _overflowed(x,y,z) = (Real)1.0;
            //peekTop(x,y,z) = topEval.magnitude();
            //peekBottom(x,y,z) = bottomEval.magnitude();
          }

          //magnitudeTrial = iterate.magnitude();

          iterate *= _expScaling;
          magnitude = iterate.magnitude();

          totalIterations++;
        }

        //Real logKnown = log(magnitude);
        //Real logTrial = log(magnitudeTrial) + log(_expScaling);

        fractal(x,y,z) = log(magnitude);
      }
    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;

/*
  // DEBUG
  if (_overflowed.sum() > 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OVERFLOW DETECTED " << endl;

    FIELDVIEW3D(_overflowed);
    //FIELDVIEW3D(peekTop);
    //FIELDVIEW3D(peekBottom);
    exit(0);
  }
  */
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogPowerRationalMap(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  if (verbose)
    cout << " Computing powerized rational Julia set ... " << flush;

  // compute a ground truth solution
#if OMP_ENABLED
#pragma omp parallel
#pragma omp for schedule(static)
#endif
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
    
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        //Real magnitudeTrial = 0;
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
#if EXPLICIT_LEADING_ROOT
          QUATERNION topEval = top.evaluatePowerFactored(iterate);
#else
          QUATERNION topEval = iterate * top.evaluatePowerFactored(iterate);
#endif
          QUATERNION bottomEval;
          
          if (bottom.totalRoots() > 0)
          {
#if EXPLICIT_LEADING_ROOT
            bottomEval = bottom.evaluatePowerFactored(iterate);
#else
            bottomEval = iterate * bottom.evaluatePowerFactored(iterate);
#endif
            iterate = (topEval / bottomEval);
          }
          else
            iterate = topEval;

          //magnitudeTrial = iterate.magnitude();

          iterate *= _expScaling;
          magnitude = iterate.magnitude();

          totalIterations++;
        }

        fractal(x,y,z) = log(magnitude);
      }

    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogRationalMapDistanceEstimated(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

/*
  vector<QUATERNION> bottomRoots;
  for (int x = 0; x < 5; x++)
    //bottomRoots.push_back(QUATERNION(-2,0,0,0));
    bottomRoots.push_back(QUATERNION(-1,0,0,0));
  _bottom = POLYNOMIAL_4D(bottomRoots);
  */

  Real escape = 20.0;
  if (verbose)
    cout << " Computing distance estimated rational Julia set ... " << flush;

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        QUATERNION prime(1,1,1,1);
        
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;
 
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          QUATERNION topEval = iterate * top.evaluateFactored(iterate);
          
          if (bottom.totalRoots() > 0)
          {
            QUATERNION bottomEval = iterate * bottom.evaluateFactored(iterate);
            iterate = _expScaling * (topEval / bottomEval);
          }
          else
            iterate = _expScaling * topEval;

          prime = 2.0 * prime * iterate;

          magnitude = iterate.magnitude();
          totalIterations++;
        }

        //Real primeMagnitude = prime.magnitude();

/*
        Real final = log(magnitude);
        Real sign = final < 0.0 ? -1 : 1;
        final = fabs(final);
        final = sign * log(fabs(log(fabs(log(final)))));
        */

        fractal(x,y,z) = log(magnitude);
        //fractal(x,y,z) = final; 
      }

    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;

  //FIELDVIEW3D(fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogRationalMap(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //const Real eps = 1e-8;
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

/*
  vector<QUATERNION> bottomRoots;
  for (int x = 0; x < 5; x++)
    //bottomRoots.push_back(QUATERNION(-2,0,0,0));
    bottomRoots.push_back(QUATERNION(-1,0,0,0));
  _bottom = POLYNOMIAL_4D(bottomRoots);
  */

  Real escape = 20.0;
  if (verbose)
    cout << " Computing rational Julia set ... " << flush;

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;
 
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          QUATERNION topEval = iterate * top.evaluateFactored(iterate);
          
          if (bottom.totalRoots() > 0)
          {
            QUATERNION bottomEval = iterate * bottom.evaluateFactored(iterate);
            iterate = _expScaling * (topEval / bottomEval);
          }
          else
            iterate = _expScaling * topEval;

          magnitude = iterate.magnitude();
          totalIterations++;
        }

        fractal(x,y,z) = log(magnitude);
      }

    if (verbose)
      if (y % (int)(yRes / 10) == 0)
        cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  if (verbose)
    cout << " done. " << endl;

  //FIELDVIEW3D(fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeKittyMap()
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

  //const Real eps = 1e-8;
  const Real coeff = exp(-108.0);
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real escape = 20.0;
  //int maxIterations = 20;
  cout << " Computing Julia set ... " << flush;

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        //QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        QUATERNION iterate(xReal, yReal, zReal, -0.1);
        
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;
 
        while (magnitude < escape && totalIterations < _maxIterations) 
        {
          // iterate = polynomial.evaluateFactored(iterate);
          iterate = coeff * iterate * _top.evaluateFactored(iterate);
          magnitude = iterate.magnitude();
          totalIterations++;
        }
        _fractal(x,y,z) = (totalIterations == _maxIterations || magnitude < escape) ? -1.0 : 1.0;
        //_fractal(x,y,z) = totalIterations;
      }

    if (y % (int)(yRes / 10) == 0)
      cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  cout << " done. " << endl;

  //FIELDVIEW3D(_fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogPolynomial()
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

  //const Real eps = 1e-8;
  //const Real coeff = exp(-108.0);
#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  //Real escape = 20.0;
  //int maxIterations = 100;
  //int maxIterations = 20;
  cout << " Computing Julia set ... " << flush;

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
  {
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;
 
        //while (magnitude < escape && totalIterations < _maxIterations) 
        {
          // iterate = polynomial.evaluateFactored(iterate);
          iterate = iterate * _top.evaluateFactoredPositive(iterate);
          magnitude = iterate.magnitude();
          //totalIterations++;
        }
        //_fractal(x,y,z) = (totalIterations == _maxIterations || magnitude < escape) ? -1.0 : 1.0;
        _fractal(x,y,z) = log(magnitude);

#ifdef QUAD_PRECISION
        if (isinfq(_fractal(x,y,z)))
#else
        if (isinf(_fractal(x,y,z)))
#endif
          _fractal(x,y,z) = 0.001234;
          //_fractal(x,y,z) = magnitude;
        //_fractal(x,y,z) = totalIterations;
      }

    if (z % (int)(zRes / 10) == 0)
      cout << 100 * ((Real)z / zRes) << "% " << flush;
  }
  cout << " done. " << endl;

  //FIELDVIEW3D(_fractal);
  //exit(0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D OPTIMIZE_3D::computeKittyMap2D(int xRes, int yRes)
{
  TIMER functionTimer(__FUNCTION__);

  //const Real eps = 1e-8;
  const Real coeff = exp(-108.0);
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  //Real zMin = -_lengths[2] * 0.5;

  Real escape = 20.0;
  //int maxIterations = 100;
  cout << " Computing Julia set ... " << flush;

  FIELD_2D final(xRes, yRes);

  // compute a ground truth solution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
  {
    for (int x = 0; x < xRes; x++)
    {
      Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
      Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
      Real zReal = 0;

      QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
      
      //Real distanceTravelled = 0.0;
      Real magnitude = iterate.magnitude();
      int totalIterations = 0;

      while (magnitude < escape && totalIterations < _maxIterations) 
      {
        // iterate = polynomial.evaluateFactored(iterate);
        iterate = coeff * iterate * _top.evaluateFactored(iterate);
        magnitude = iterate.magnitude();
        totalIterations++;
      }
      final(x,y) = (totalIterations == _maxIterations || magnitude < escape) ? -1.0 : 1.0;
      //_fractal(x,y,z) = (totalIterations == _maxIterations || magnitude < escape) ? -1.0 : 1.0;
      //_fractal(x,y,z) = totalIterations;
    }
    if (y % (int)(yRes / 10) == 0)
      cout << 100 * ((Real)y / yRes) << "% " << flush;
  }
  cout << " done. " << endl;

  return final;  
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::initMap()
{
  vector<QUATERNION> topRoots;

  //topRoots.push_back(QUATERNION(-0.9, 0.2, -0.8, 0.1));
  topRoots.push_back(QUATERNION(0.2, 0.2, 0.3, 0.1));
  topRoots.push_back(QUATERNION(0.123, 0.133, 0.182, 0.1));
  topRoots.push_back(QUATERNION(0.23, 0.33, 0.12, 0));

  _top = POLYNOMIAL_4D(topRoots);

  vector<QUATERNION> bottomRoots;
  bottomRoots.push_back(QUATERNION(0.1, 0.1, 0.15, 0.05));
  bottomRoots.push_back(QUATERNION(0.1, 0.1, 0.15, 0.05));
  //bottomRoots.push_back(QUATERNION(-0.1, 0.2, -0.15, 0.12));
  
  _bottom = POLYNOMIAL_4D(bottomRoots);
}

///////////////////////////////////////////////////////////////////////
// compute the score of _fractal versus _distanceField
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBandedScore(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal)
{
  TIMER functionTimer(__FUNCTION__);
  //static bool firstCall = true;

  //int binaryBandwidth = 10;
  //int curvatureBandwidth = 2;
  //int maxBandwidth = (binaryBandwidth > curvatureBandwidth) ? binaryBandwidth : curvatureBandwidth;

  FIELD_3D aux = fractal;
  computeMap(polynomial, fractal, aux);

  Real binaryScore = computeBinaryBandedScore(fractal);
  Real curvatureScore = computeCurvatureBandedScore(fractal);

  return (binaryScore + curvatureScore) * 0.5;
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " CPU BINARY ONLY " << endl;
  //return binaryScore; 
}

///////////////////////////////////////////////////////////////////////
// compute the score of _fractal versus _distanceField
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeScore(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal)
{
  TIMER functionTimer(__FUNCTION__);
  //static bool firstCall = true;

  //int binaryBandwidth = 10;
  //int curvatureBandwidth = 2;
  /*
  if (firstCall)
  {
    cout << " Using banded binary plus banded curvature scoring function " << endl;
    cout << "   binary bandwidth: " << binaryBandwidth << endl; 
    cout << "   curvature bandwidth: " << curvatureBandwidth << endl; 
    firstCall = false;
  }
  */

  //int maxBandwidth = (_binaryBandwidth > _curvatureBandwidth) ? _binaryBandwidth : _curvatureBandwidth;

  FIELD_3D aux = fractal;
  computeMap(polynomial, fractal, aux);

  Real binaryScore = computeBinaryScore(fractal);

  //return binaryScore;
  Real curvatureScore = computeCurvatureScore(fractal);
  //return curvatureScore;

/*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "================================================= " << endl;
  cout << " CPU binary score: " << binaryScore << endl;
  cout << " CPU curvature score: " << curvatureScore << endl;
  */

  return (binaryScore + curvatureScore) * 0.5;
}

///////////////////////////////////////////////////////////////////////
// compute the score with respect to a distance field
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBinaryScore(const FIELD_3D& fractal) const
{
  //assert(_distanceField.totalCells() == fractal.totalCells());

  FIELD_3D product = fractal;
  FIELD_3D threshold = fractal;

  TIMER functionTimer(__FUNCTION__);
  // compute the score
  Real score = 0;
  for (int x = 0; x < fractal.totalCells(); x++)
  {
    VEC3F position = fractal.cellCenter(x);
    double distanceThresholded = (_distanceField(position) < 0) ? -1 : 1;
    double fractalThresholded = (fractal[x] > 0) ? -1 : 1;
    product[x] = fractalThresholded * distanceThresholded;
    threshold[x] = distanceThresholded;
    score += fractalThresholded * distanceThresholded;
  }

  return score / fractal.totalCells();
}

///////////////////////////////////////////////////////////////////////
// compute the score with respect to a distance field
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeInverseDistanceScore(const FIELD_3D& fractal)
{
  //assert(_distanceField.totalCells() == fractal.totalCells());

  FIELD_3D product = fractal;
  FIELD_3D threshold = fractal;

  TIMER functionTimer(__FUNCTION__);

  // make sure that inverse field indeed exists
  if (_inverseDistanceField.totalCells() != _distanceField.totalCells())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Recomputing the inverse distance field ... " << endl;
    _inverseDistanceField = _distanceField.invertedDistance();

    FIELD_3D temp = _inverseDistanceField;
    temp.absoluteValue();
    _inverseFieldSum = temp.sum();

    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //FIELDVIEW3D(_inverseDistanceField);
    //exit(0);
  }

  // compute the score
  Real score = 0;
  for (int x = 0; x < fractal.totalCells(); x++)
  {
    VEC3F position = fractal.cellCenter(x); 
    //double inverseDistance = _inverseDistanceField[x];
    double inverseDistance = _inverseDistanceField(position);
    double fractalThresholded = (fractal[x] > 0) ? -1 : 1;
    product[x] = fractalThresholded * inverseDistance;
    threshold[x] = inverseDistance;
    score += fractalThresholded * inverseDistance;
  }

  return score / _inverseFieldSum;
}

///////////////////////////////////////////////////////////////////////
// compute the score with respect to a distance field
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBinaryBandedScore(const FIELD_3D& fractal) const
{
  assert(_distanceField.totalCells() == fractal.totalCells());

  FIELD_3D product = fractal;
  FIELD_3D threshold = fractal;
  int maxRes = (_xRes > _yRes) ? _xRes : _yRes;
  maxRes = (maxRes > _zRes) ? maxRes : _zRes;

  TIMER functionTimer(__FUNCTION__);
  // compute the score
  Real score = 0;
  int totalBanded = 0;
  for (int x = 0; x < fractal.totalCells(); x++)
  {
    double sign = (_distanceField[x] < 0) ? 1 : -1; 
    int cellsInside = fabs(_distanceField[x] * maxRes);

    if (sign < 0 && cellsInside > _binaryBandwidth)
      continue;
    if (sign > 0 && cellsInside > _binaryBandwidth)
      continue;

    double distanceThresholded = (_distanceField[x] < 0) ? -1 : 1;
    double fractalThresholded = (fractal[x] > 0) ? -1 : 1;
    product[x] = fractalThresholded * distanceThresholded;
    threshold[x] = distanceThresholded;
    score += fractalThresholded * distanceThresholded;
    totalBanded++;
  }

  return score / totalBanded;
}

///////////////////////////////////////////////////////////////////////
// compute the center of the current distance field
///////////////////////////////////////////////////////////////////////
VEC3F OPTIMIZE_3D::computeCenter(const FIELD_3D& field)
{
  assert(field.totalCells() > 0);

  int zRes = field.zRes();
  int yRes = field.yRes();
  int xRes = field.xRes();

  VEC3F sum;
  int totalSummed = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        if (field(x,y,z) > 0) continue;

        sum += field.cellCenter(x,y,z);
        totalSummed++;
      }
  sum *= 1.0 / totalSummed;

  return sum;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeCenterAndRadius(const FIELD_3D& field, VEC3F& center, Real& radius)
{
  assert(field.totalCells() > 0);

  center = computeCenter(field);

  int index = field.cellIndex(center);
  VEC3I middle;
  middle[2] = index / field.slabSize();
  index -= middle[2] * field.slabSize();
  middle[1] = index / field.xRes(); 
  middle[0] = index % field.xRes();
  int zRes = field.zRes();
  int yRes = field.yRes();
  int xRes = field.xRes();

  // get the radii in the x direction
  int xPlus = 0;
  for (int x = middle[0]; x < xRes; x++)
  {
    if (field(x, middle[1], middle[2]) > 0) break;
    xPlus++;
  }
  int xMinus = 0;
  for (int x = middle[0]; x >= 0; x--)
  {
    if (field(x, middle[1], middle[2]) > 0) break;
    xMinus++;
  }
   
  // get the radii in the y direction
  int yPlus = 0;
  for (int y = middle[1]; y < yRes; y++)
  {
    if (field(middle[0], y, middle[2]) > 0) break;
    yPlus++;
  }
  int yMinus = 0;
  for (int y = middle[1]; y >= 0; y--)
  {
    if (field(middle[0], y, middle[2]) > 0) break;
    yMinus++;
  }

  // get the radii in the z direction
  int zPlus = 0;
  for (int z = middle[2]; z < zRes; z++)
  {
    if (field(middle[0], middle[1], z) > 0) break;
    zPlus++;
  }
  int zMinus = 0;
  for (int z = middle[2]; z >= 0; z--)
  {
    if (field(middle[0], middle[1], z) > 0) break;
    zMinus++;
  }

  cout << " x plus minus: " << xMinus << " " << xPlus << endl;
  cout << " y plus minus: " << yMinus << " " << yPlus << endl;
  cout << " z plus minus: " << zMinus << " " << zPlus << endl;

  radius = (xMinus > xPlus) ? xMinus : xPlus;
  radius = (radius > yMinus) ? radius : yMinus;
  radius = (radius > yPlus) ? radius : yPlus;
  radius = (radius > zMinus) ? radius : zMinus;
  radius = (radius > zPlus) ? radius : zPlus;
 
  int maxRes = (xRes > yRes) ? xRes : yRes;
  maxRes = (maxRes > zRes) ? maxRes : zRes;

  radius *= 1.0 / maxRes;
  radius *= field.maxLength();
}

///////////////////////////////////////////////////////////////////////
// add to a real using Kahan summation
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::kahanAdd(Real& original, const Real& add, Real& error) const
{
  Real y = add - error;
  Real t = original + y;
  error = (t - original) - y;
  original = t;
}

///////////////////////////////////////////////////////////////////////
// compute the score with respect to a curvature field
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeCurvatureScore(const FIELD_3D& fractal) const
{
  assert(_curvatureField.totalCells() == fractal.totalCells());

  TIMER functionTimer(__FUNCTION__);
  int maxRes = (_xRes > _yRes) ? _xRes : _yRes;
  maxRes = (maxRes > _zRes) ? maxRes : _zRes;

  // compute the score
  Real score = 0;
  Real scoreKahan = 0;
  Real sum = 0;
  Real sumKahan = 0;
  FIELD_3D product = fractal;
  for (int x = 0; x < fractal.totalCells(); x++)
  {
    //Real sign = (_distanceField[x] < 0) ? -1 : 1; 
    Real sign = (_distanceField[x] < 0) ? 1 : -1; 
    //int cellsInside = fabs(_distanceField[x] * maxRes);

    // remember the fractal is the inverted sign, not the distance field realted stuff
    Real adding = sign * _curvatureField[x];
    product[x] = adding * fractal[x];
    //score += adding * fractal[x];
    //sum += fabs(adding);
   
    kahanAdd(score, adding * fractal[x], scoreKahan);
    kahanAdd(sum, fabs(adding), sumKahan);
  }

  //FIELDVIEW3D(product);
  //FIELDVIEW3D(_curvatureField);
  //product.write("cpuResult.field3d");
  //exit(0);

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " Curvature field sum:        " << sum << endl;

  return score / sum;
  //return sum;
}

///////////////////////////////////////////////////////////////////////
// compute the score with respect to a curvature field
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeCurvatureBandedScore(const FIELD_3D& fractal) const
{
  assert(_curvatureField.totalCells() == fractal.totalCells());

  TIMER functionTimer(__FUNCTION__);
  int maxRes = (_xRes > _yRes) ? _xRes : _yRes;
  maxRes = (maxRes > _zRes) ? maxRes : _zRes;

  // compute the score
  Real score = 0;
  Real scoreKahan = 0;
  Real sum = 0;
  Real sumKahan = 0;
  FIELD_3D product = fractal;
  product.clear();
  for (int x = 0; x < fractal.totalCells(); x++)
  {
    //Real sign = (_distanceField[x] < 0) ? -1 : 1; 
    Real sign = (_distanceField[x] < 0) ? 1 : -1; 
    int cellsInside = fabs(_distanceField[x] * maxRes);

    if (sign < 0 && cellsInside > _curvatureBandwidth)
      continue;
    if (sign > 0 && cellsInside > _curvatureBandwidth)
      continue;

    // remember the fractal is the inverted sign, not the distance field realted stuff
    Real adding = sign * _curvatureField[x];
    product[x] = adding * fractal[x];
    //score += adding * fractal[x];
    //sum += fabs(adding);
   
    kahanAdd(score, adding * fractal[x], scoreKahan);
    kahanAdd(sum, fabs(adding), sumKahan);
  }

  //FIELDVIEW3D(product);
  //FIELDVIEW3D(_curvatureField);
  //product.write("cpuResult.field3d");
  //exit(0);

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " Curvature field sum:        " << sum << endl;

  return score / sum;
  //return sum;
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeCenteredGradientCUDA()
{
#ifdef USING_CUDA
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  Real dx = _lengths[0] / _xRes;
  Real eps = dx;

  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D top(_top);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      top = _top;

      QUATERNION newRoot = top.roots()[x];
      newRoot[y] += eps;
      top.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeScoreCUDA(top, fractal);

      // reset the polynomial
      top = _top;
      newRoot = top.roots()[x];
      newRoot[y] -= eps;
      top.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeScoreCUDA(top, fractal);
      int entry = 4 * x + y;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " left: " << leftScore << " right: " << rightScore << " diff: " << (rightScore - leftScore) << " final: " << gradient[entry] << endl;
    }
  }

  return gradient;
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT SUPPORTED " << endl;
  return VECTOR();
#endif
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeBandedCenteredGradientCUDA()
{
#ifdef USING_CUDA
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  Real dx = _lengths[0] / _xRes;
  Real eps = dx;

  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D top(_top);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
      {
        // reset the polynomial
        top = _top;

        QUATERNION newRoot = top.roots()[x];
        newRoot[y] += eps;
        top.modifyRoot(x, newRoot);

        // compute its new score
        Real rightScore = computeBandedScoreCUDA(top, fractal);

        // reset the polynomial
        top = _top;
        newRoot = top.roots()[x];
        newRoot[y] -= eps;
        top.modifyRoot(x, newRoot);
        
        // compute its new score
        Real leftScore = computeBandedScoreCUDA(top, fractal);
        int entry = 4 * x + y;
        gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
        //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        //cout << " left: " << leftScore << " right: " << rightScore << " diff: " << (rightScore - leftScore) << " final: " << gradient[entry] << endl;
      }
  }

  return gradient;
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT SUPPORTED " << endl;
  return VECTOR();
#endif
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeBandedCenteredGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

  //int maxThreads = omp_get_max_threads();

//#pragma omp parallel
//#pragma omp for  schedule(static)
  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _top;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeBandedScore(polynomial, fractal);

      // reset the polynomial
      polynomial = _top;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeBandedScore(polynomial, fractal);
      int entry = 4 * x + y;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " left: " << leftScore << " right: " << rightScore << " diff: " << (rightScore - leftScore) << " final: " << gradient[entry] << endl;
    }

    fractal.clear();
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeCenteredGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

  //int maxThreads = omp_get_max_threads();

//#pragma omp parallel
//#pragma omp for  schedule(static)
  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _top;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeScore(polynomial, fractal);

      // reset the polynomial
      polynomial = _top;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeScore(polynomial, fractal);
      int entry = 4 * x + y;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " left: " << leftScore << " right: " << rightScore << " diff: " << (rightScore - leftScore) << " final: " << gradient[entry] << endl;
    }

    fractal.clear();
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeLogCenteredGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

  //int maxThreads = omp_get_max_threads();

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _top;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeLogScore(polynomial, fractal);

      // reset the polynomial
      polynomial = _top;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeLogScore(polynomial, fractal);
      int entry = 4 * x + y;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " left: " << leftScore << " right: " << rightScore << " diff: " << (rightScore - leftScore) << " final: " << gradient[entry] << endl;
    }

    fractal.clear();
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeLogRationalCenteredGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  int topRoots = _top.totalRoots();
  int bottomRoots = _bottom.totalRoots();
  int totalRoots = topRoots + bottomRoots;
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < topRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _top;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeLogRationalScore(polynomial, _bottom, fractal);

      // reset the polynomial
      polynomial = _top;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeLogRationalScore(polynomial, _bottom, fractal);
      int entry = 4 * x + y;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
    }
    fractal.clear();
  }

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < bottomRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_bottom);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _bottom;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeLogRationalScore(_top, polynomial, fractal);

      // reset the polynomial
      polynomial = _bottom;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeLogRationalScore(_top, polynomial, fractal);
      int entry = 4 * x + y + 4 * topRoots;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
    }
    fractal.clear();
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeLogRationalBottomCenteredGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  int bottomRoots = _bottom.totalRoots();
  int totalRoots = bottomRoots;
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

//#pragma omp parallel
//#pragma omp for  schedule(static)
  for (int x = 0; x < bottomRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_bottom);
    FIELD_3D fractal(_fractal);
//#pragma omp parallel
//#pragma omp for  schedule(static)
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _bottom;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeLogRationalBottomScore(polynomial, fractal);

      // reset the polynomial
      polynomial = _bottom;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeLogRationalBottomScore(polynomial, fractal);
      int entry = 4 * x + y;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
    }
    fractal.clear();
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computePowerCenteredGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  int topRoots = _top.totalRoots();
  int bottomRoots = _bottom.totalRoots();
  int totalRoots = topRoots + bottomRoots;
  VECTOR gradient(totalRoots);

  //const Real dx = _lengths[0] / _xRes;
  //const Real dx = 0.5;
  const Real dx = _powerDx;
  const Real eps = dx;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < topRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D fractal(_fractal);

    Real power = polynomial.powers()[x];

    // compute its new score
    polynomial.changePower(x, power + eps);
    Real rightScore = computePowerScore(polynomial, _bottom, fractal);
    
    // compute its new score
    polynomial.changePower(x, power - eps);
    Real leftScore = computePowerScore(polynomial, _bottom, fractal);
    
    gradient[x] = (rightScore - leftScore) / (2.0 * eps);
  }

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < bottomRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_bottom);
    FIELD_3D fractal(_fractal);

    Real power = polynomial.powers()[x];

    // compute its new score
    polynomial.changePower(x, power + eps);
    Real rightScore = computePowerScore(_top, polynomial, fractal);

    // compute its new score
    polynomial.changePower(x, power - eps);
    Real leftScore = computePowerScore(_top, polynomial, fractal);
    gradient[topRoots + x] = (rightScore - leftScore) / (2.0 * eps);
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// energy derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeLogBandedCenteredGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

  //int maxThreads = omp_get_max_threads();

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D fractal(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _top;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      Real rightScore = computeLogBandedScore(polynomial, fractal);

      // reset the polynomial
      polynomial = _top;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      Real leftScore = computeLogBandedScore(polynomial, fractal);
      int entry = 4 * x + y;
      gradient[entry] = (rightScore - leftScore) / (2.0 * eps);
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " left: " << leftScore << " right: " << rightScore << " diff: " << (rightScore - leftScore) << " final: " << gradient[entry] << endl;
    }

    fractal.clear();
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// self derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeCenteredSelfGradientCUDA()
{
#ifdef USING_CUDA
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

  //int maxThreads = omp_get_max_threads();

//#pragma omp parallel
//#pragma omp for  schedule(static)
  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D auxiliary(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _top;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      FIELD_3D fractalRight(_fractal);
      computeMapCUDA(polynomial, fractalRight, auxiliary);

      // reset the polynomial
      polynomial = _top;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      FIELD_3D fractalLeft(_fractal);
      computeMapCUDA(polynomial, fractalLeft, auxiliary);

      // get the diff
      FIELD_3D diff = fractalRight - fractalLeft;
      Real diffNorm = sqrt(diff.sumSq()) / diff.totalCells();

      int entry = 4 * x + y;
      gradient[entry] = diffNorm / (2.0 * eps);
      
      //if (fabs(gradient[entry]) < 1e-8)
      //  cout << " Dud DOF at root " << x << ": " << _top.roots()[x] << " entry: " << entry << endl;
    }
  }

  return gradient;
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT SUPPORTED " << endl;
  return VECTOR();
#endif
}

///////////////////////////////////////////////////////////////////////
// self derivatives
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeCenteredSelfGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  int totalRoots = _top.totalRoots();
  int totalDOFs = totalRoots * 4;
  VECTOR gradient(totalDOFs);

  const Real dx = _lengths[0] / _xRes;
  const Real eps = dx;

  //int maxThreads = omp_get_max_threads();

//#pragma omp parallel
//#pragma omp for  schedule(static)
  for (int x = 0; x < totalRoots; x++)
  {
    POLYNOMIAL_4D polynomial(_top);
    FIELD_3D auxiliary(_fractal);
    for (int y = 0; y < 4; y++)
    {
      // reset the polynomial
      polynomial = _top;

      QUATERNION newRoot = polynomial.roots()[x];
      newRoot[y] += eps;
      polynomial.modifyRoot(x, newRoot);

      // compute its new score
      FIELD_3D fractalRight(_fractal);
      computeMap(polynomial, fractalRight, auxiliary);

      // reset the polynomial
      polynomial = _top;
      newRoot = polynomial.roots()[x];
      newRoot[y] -= eps;
      polynomial.modifyRoot(x, newRoot);
      
      // compute its new score
      FIELD_3D fractalLeft(_fractal);
      computeMap(polynomial, fractalLeft, auxiliary);

      // get the diff
      FIELD_3D diff = fractalRight - fractalLeft;
      Real diffNorm = sqrt(diff.sumSq()) / diff.totalCells();

      int entry = 4 * x + y;
      gradient[entry] = diffNorm / (2.0 * eps);
      
      //if (fabs(gradient[entry]) < 1e-8)
      //  cout << " Dud DOF at root " << x << ": " << _top.roots()[x] << " entry: " << entry << endl;
    }
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// write current state to file
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::write(const string& filename) const
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OPTIMIZE_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Writing file " << filename.c_str() << " ... " << flush;

  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);

  _center.write(file);
  _lengths.write(file);

  _top.write(file);
  _bottom.write(file);

  double expScaling = _expScaling;
  fwrite((void*)&expScaling, sizeof(double), 1, file);

  //fwrite((void*)&_flowRadius, sizeof(Real), 1, file);
  //_flowCenter.write(file);

  int size = _historyTop.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  // TODO: Make Eigen write to a file
  //for (int x = 0; x < size; x++)
  //  _historyTop[x].write(file);

  _fractal.write(file);
  _distanceField.write(file);
  _curvatureField.write(file);

  fclose(file);
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// write current state to file
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::read(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " OPTIMIZE_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Reading file " << filename.c_str() << " ... " << flush;

  // read dimensions
  fread((void*)&_xRes, sizeof(int), 1, file);
  fread((void*)&_yRes, sizeof(int), 1, file);
  fread((void*)&_zRes, sizeof(int), 1, file);

  cout << " Read in res: " << _xRes << " " << _yRes << " " << _zRes << endl;

  _center.read(file);
  _lengths.read(file);

  _top.read(file);
  _bottom.read(file);

  double expScaling;
  fread((void*)&expScaling, sizeof(double), 1, file);
  _expScaling = expScaling;

  //fread((void*)&_flowRadius, sizeof(Real), 1, file);
  //_flowCenter.read(file);

  int size;
  fread((void*)&size, sizeof(int), 1, file);
  _historyTop.resize(size);
  // TODO: Have Eigen read from a file
  //for (int x = 0; x < size; x++)
  //  _historyTop[x].read(file);

  _fractal.read(file);
  _auxiliary = _fractal;
  _distanceField.read(file);
  _curvatureField.read(file);

  fclose(file);
  cout << " done. " << endl;
#if 1
  //if (_distanceField.totalCells() > 0 && _curvatureField.totalCells() > 0)
  //  setScoreFields(_distanceField, _curvatureField);
  setScoreFields(_distanceField, _distanceField);

  if (_inverseDistanceField.totalCells() <= 0)
  {
    // this can cause precision issues, so punt on this for now
    /*
    // try to read in an external one
    string inverseFile = filename + string(".inverse.field");
    FILE* checkExists = fopen(inverseFile.c_str(), "rb");
  
    // if it's there, read it 
    if (checkExists != NULL)
    {
      fclose(checkExists);
      cout << " External precomputed inverse field found! " << endl;
      _inverseDistanceField.read(inverseFile);
    }
    // give up, recompute
    else
    */
    {
      cout << " No inverse field found! Recomputing ... " << flush;
      _inverseDistanceField = _distanceField.invertedDistanceNormalized();
      FIELD_3D temp = _inverseDistanceField;
      temp.absoluteValue();
      _inverseFieldSum = temp.sum();
      cout << " done. " << endl;
      //_inverseDistanceField.write(inverseFile);
    }
  }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  computeWeightField();

  cout << " inverse distance field sum: " << _inverseFieldSum << endl;

  _lastLoadedFilename = filename;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
#endif
}

///////////////////////////////////////////////////////////////////////
// read the distance field from a file
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::readDistanceField(const string& filename)
{
  _distanceField.read(filename);

  //if (_distanceField.totalCells() > 0 && _curvatureField.totalCells() > 0)
  //  setScoreFields(_distanceField, _curvatureField);
  setScoreFields(_distanceField, _distanceField);

  if (_inverseDistanceField.totalCells() <= 0)
  {
    // this can cause precision issues, so punt on this for now
    /*
    // try to read in an external one
    string inverseFile = filename + string(".inverse.field");
    FILE* checkExists = fopen(inverseFile.c_str(), "rb");
  
    // if it's there, read it 
    if (checkExists != NULL)
    {
      fclose(checkExists);
      cout << " External precomputed inverse field found! " << endl;
      _inverseDistanceField.read(inverseFile);
    }
    // give up, recompute
    else
    */
    {
      cout << " No inverse field found! Recomputing ... " << flush;
      _inverseDistanceField = _distanceField.invertedDistanceNormalized();
      FIELD_3D temp = _inverseDistanceField;
      temp.absoluteValue();
      _inverseFieldSum = temp.sum();
      cout << " done. " << endl;
      //_inverseDistanceField.write(inverseFile);
    }
  }

  computeWeightField();

  cout << " inverse distance field sum: " << _inverseFieldSum << endl;
}

///////////////////////////////////////////////////////////////////////
// set the score fields and extract the field dimensions
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::setScoreFields(const FIELD_3D& distanceField, const FIELD_3D& curvatureField)
{
  assert(distanceField.xRes() == curvatureField.xRes());
  assert(distanceField.yRes() == curvatureField.yRes());
  assert(distanceField.zRes() == curvatureField.zRes());

  _xRes = distanceField.xRes();
  _yRes = distanceField.yRes();
  _zRes = distanceField.zRes();

  // NEW: assume the distance field actually has decent spatial information
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " RESETTING CENTER AND LENGTH BASED ON DISTANCE FILE" << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  _center = distanceField.center();
  _lengths = distanceField.lengths();

  cout << " Using field center: " << distanceField.center() << endl;
  cout << " Using field lengths: " << distanceField.lengths() << endl;

  _distanceField = FIELD_3D(distanceField.dataConst(), _xRes, _yRes, _zRes, _center, _lengths);
  _curvatureField = FIELD_3D(curvatureField.dataConst(), _xRes, _yRes, _zRes, _center, _lengths);

  _fractal = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _auxiliary = _fractal;
  
  //setScoreFieldsCUDA(_distanceField, _curvatureField);
}

///////////////////////////////////////////////////////////////////////
// set the score fields and extract the field dimensions
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::setScoreFields(const FIELD_3D& distanceField, const FIELD_3D& curvatureField, const FIELD_3D& inverseField)
{
  assert(distanceField.xRes() == curvatureField.xRes());
  assert(distanceField.yRes() == curvatureField.yRes());
  assert(distanceField.zRes() == curvatureField.zRes());

  _xRes = distanceField.xRes();
  _yRes = distanceField.yRes();
  _zRes = distanceField.zRes();

  // NEW: assume the distance field actually has decent spatial information
  _center = distanceField.center();
  _lengths = distanceField.lengths();

  cout << " Using field center: " << distanceField.center() << endl;
  cout << " Using field lengths: " << distanceField.lengths() << endl;

  _distanceField = FIELD_3D(distanceField.dataConst(), _xRes, _yRes, _zRes, _center, _lengths);
  _curvatureField = FIELD_3D(curvatureField.dataConst(), _xRes, _yRes, _zRes, _center, _lengths);

  _inverseDistanceField = FIELD_3D(inverseField.dataConst(), _xRes, _yRes, _zRes, _center, _lengths);
  FIELD_3D temp = _inverseDistanceField;
  temp.absoluteValue();
  _inverseFieldSum = temp.sum();

  computeWeightField();

  _fractal = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _auxiliary = _fractal;
  //setScoreFieldsCUDA(_distanceField, _curvatureField);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeScore(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setRoots(roots);
  return computeScore(_top, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogRationalScore(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setRationalRoots(roots);
  return computeLogRationalScore(_top, _bottom, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogRationalBottomScoreVarying(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setBottomRoots(roots);
  return computeLogRationalBottomScore(_bottom, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogRationalBottomScoreVarying()
{
  TIMER functionTimer(__FUNCTION__);
  return computeLogRationalBottomScore(_bottom, _fractal);
}

#define USING_WEIGHT 1

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogPowerRationalScore()
{
  TIMER functionTimer(__FUNCTION__);
  // DEBUG: see if caching is working
  //computeLogPowerRationalMap(_top, _bottom, _fractal);
  computeLogPowerRationalMapCached(_top, _bottom, _fractal);
  //return computeBinaryScore(_fractal);
  
  //FIELD_3D inverseDistanceField = _distanceField.invertedDistance();
  //Real inverseSum = inverseDistanceField.absSum();
 
  return weightScore(_fractal, _weightField);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeScoreAP()
{
  TIMER functionTimer(__FUNCTION__);
  computeRationalMapAP(_top, _bottom, _fractal);
 
  return weightScore(_fractal, _weightField);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogScaledPowerRationalScore()
{
  TIMER functionTimer(__FUNCTION__);
  // DEBUG: see if caching is working
  //computeLogPowerRationalMap(_top, _bottom, _fractal);
  computeLogScaledPowerRationalMapCached(_top, _bottom, _fractal);
  //return computeBinaryScore(_fractal);
  
  //FIELD_3D inverseDistanceField = _distanceField.invertedDistance();
  //Real inverseSum = inverseDistanceField.absSum();
 
  return weightScore(_fractal, _weightField);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeExpGradientAP()
{
  TIMER functionTimer(__FUNCTION__);

  // DEBUG: Assume that this was called by the scoring function
  //computeRationalMapAP(_top, _bottom, _fractal);
 
  Real gradientScore = 0;
  FIELD_3D gradientField(_fractal);
  gradientField = 0;

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // Shouldn't need to do anything special for overflows here.
    // Any cell with an overflow was already recomputed
    /*
    // ADDING OVERFLOW DETECTION
    if (_overflowed[x] > 0.5) continue;

    // TODO AP: Add another overflow path here?
    */

    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real weight = _weightField(cellPosition);

    Real tanhVersionBefore = tanh(_tanhAlpha * (_fractal[x] - _tanhThreshold));
    Real S = tanhVersionBefore;

    Real currentScore = -_tanhAlpha * weight * ((Real)1.0 - S * S);

    gradientScore += currentScore;
    gradientField[x] = currentScore;
  }

  return gradientScore;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogPowerScaledRationalScoreExpGradient()
{
  TIMER functionTimer(__FUNCTION__);
  // DEBUG: see if caching is working
  //computeLogPowerRationalMap(_top, _bottom, _fractal);
  computeLogScaledPowerRationalMapCached(_top, _bottom, _fractal);
 
  return weightScoreConformalGradient(_fractal, _weightField);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogPowerRationalScoreExpGradient()
{
  TIMER functionTimer(__FUNCTION__);
  // DEBUG: see if caching is working
  //computeLogPowerRationalMap(_top, _bottom, _fractal);
  computeLogPowerRationalMapCached(_top, _bottom, _fractal);
 
  return weightScoreConformalGradient(_fractal, _weightField);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogRationalBottomScore(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setBottomRoots(roots);
  return computeLogRationalBottomScore(_bottom, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computePowerScore(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setPowers(roots);
  return computePowerScore(_top, _bottom, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBandedScore(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setRoots(roots);
  return computeBandedScore(_top, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeScoreCUDA(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setRoots(roots);
  return computeScoreCUDA(_top, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBandedScoreCUDA(const VECTOR& roots)
{
  TIMER functionTimer(__FUNCTION__);
  setRoots(roots);
  return computeBandedScoreCUDA(_top, _fractal);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::addToHistory()
{
  VECTOR newHistory = getRoots();
  _historyTop.push_back(newHistory);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeScoreCUDA(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal)
{
#ifdef USING_CUDA
  TIMER functionTimer(__FUNCTION__);
  static bool firstCall = true;

  int binaryBandwidth = 10;
  int curvatureBandwidth = 2;
  /*
  if (firstCall)
  {
    cout << " Using banded binary plus banded curvature scoring function " << endl;
    cout << "   binary bandwidth: " << binaryBandwidth << endl; 
    cout << "   curvature bandwidth: " << curvatureBandwidth << endl; 
    firstCall = false;
  }
  */

  int maxBandwidth = (binaryBandwidth > curvatureBandwidth) ? binaryBandwidth : curvatureBandwidth;

  //const int maxIterations = 1000;
  const int maxIterations = 100;
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //VEC3F recenter = _flowCenter;
  //Real rescale = _flowRadius;

  const Real eps = 1e-8;
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5; 

  Real escape = 10.0;

  Real4* roots = new Real4[polynomial.totalRoots()];
  for (unsigned int x = 0; x < polynomial.totalRoots(); x++)
  {
    roots[x].x = polynomial.roots()[x].x();
    roots[x].y = polynomial.roots()[x].y();
    roots[x].z = polynomial.roots()[x].z();
    roots[x].w = polynomial.roots()[x].w();
  }
  setTopRoots(roots, polynomial.totalRoots());
  delete[] roots;

  setConsts(xRes, yRes, zRes,
            _center[0], _center[1], _center[2],
            _lengths[0], _lengths[1], _lengths[2], 
            escape, _maxIterations);

  return runScore(_distanceField, _curvatureField);
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT SUPPORTED " << endl;
  return 0;
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBandedScoreCUDA(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal)
{
#ifdef USING_CUDA
  TIMER functionTimer(__FUNCTION__);
  static bool firstCall = true;

  int binaryBandwidth = 10;
  int curvatureBandwidth = 2;

  int maxBandwidth = (binaryBandwidth > curvatureBandwidth) ? binaryBandwidth : curvatureBandwidth;

  //const int maxIterations = 1000;
  const int maxIterations = 100;
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //VEC3F recenter = _flowCenter;
  //Real rescale = _flowRadius;

  const Real eps = 1e-8;
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5; 

  Real escape = 10.0;

  Real4* roots = new Real4[polynomial.totalRoots()];
  for (unsigned int x = 0; x < polynomial.totalRoots(); x++)
  {
    roots[x].x = polynomial.roots()[x].x();
    roots[x].y = polynomial.roots()[x].y();
    roots[x].z = polynomial.roots()[x].z();
    roots[x].w = polynomial.roots()[x].w();
  }
  setTopRoots(roots, polynomial.totalRoots());
  delete[] roots;

  setConsts(xRes, yRes, zRes,
            _center[0], _center[1], _center[2],
            _lengths[0], _lengths[1], _lengths[2], 
            escape, _maxIterations, _binaryBandwidth, _curvatureBandwidth);

  return runScoreBanded(_distanceField, _curvatureField);
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT SUPPORTED " << endl;
  return 0;
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeMapCUDA(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal, FIELD_3D& auxiliary)
{
#ifdef USING_CUDA
  TIMER functionTimer(__FUNCTION__);

  //const int maxIterations = 1000;
  const int maxIterations = 100;
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();

  //VEC3F recenter = _flowCenter;
  //Real rescale = _flowRadius;

  const Real eps = 1e-8;
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5; 

  Real escape = 10.0;

  Real4* roots = new Real4[polynomial.totalRoots()];
  for (unsigned int x = 0; x < polynomial.totalRoots(); x++)
  {
    roots[x].x = polynomial.roots()[x].x();
    roots[x].y = polynomial.roots()[x].y();
    roots[x].z = polynomial.roots()[x].z();
    roots[x].w = polynomial.roots()[x].w();
  }
  setTopRoots(roots, polynomial.totalRoots());
  delete[] roots;

  setConsts(xRes, yRes, zRes, 
            _center[0], _center[1], _center[2],
            _lengths[0], _lengths[1], _lengths[2], 
            escape, _maxIterations);

  runMap(fractal);
  //runBinaryScore(_distanceField);
  //runScore(_distanceField, _curvatureField);
  //exit(0);
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugCUDA()
{
#ifdef USING_CUDA
/*
  FIELD_3D cudaFractal(100,100,100);
  FIELD_3D cudaAux(100,100,100);
  computeMapCUDA(cudaFractal, cudaAux);

  FIELD_3D fractal(100,100,100);
  FIELD_3D aux(100,100,100);
  computeMap(fractal, aux, true);

  FIELD_3D diff = cudaFractal - fractal;
  Real diffNorm = diff.sumSq();
  cout << " CUDA sum: " << cudaFractal.sumSq() << endl;
  cout << " CPU sum: " << fractal.sumSq() << endl;
  cout << " Diff: " << diffNorm << endl;

  if (diffNorm > 1e-5)
  {
    FIELDVIEW3D(cudaFractal);
    FIELDVIEW3D(fractal);
    FIELDVIEW3D(cudaFractal - fractal);
  }
  */
  cout << " CUDA ... " << endl;
  VECTOR gradientCUDA = computeCenteredGradientCUDA();
  cout << " CPU ... " << flush;
  VECTOR gradient = computeCenteredGradient();
  VECTOR diff = gradient - gradientCUDA;

  cout << " gradient: " << gradient << endl;
  cout << " CUDA: " << gradientCUDA << endl;
  cout << " diff: " << diff << endl;
  cout << " diff sum: " << diff.norm2() << endl;

  TIMER::printTimings();
#endif
}

///////////////////////////////////////////////////////////////////////
// verify that CUDA is still consistent
///////////////////////////////////////////////////////////////////////
bool OPTIMIZE_3D::verifyCUDA()
{
  bool final = true;

  /*
  computeMap();
  FIELD_3D cpu = _fractal;

  computeMapCUDA();
  FIELD_3D cuda = _fractal;
  FIELD_3D fractalDiff = cpu - cuda;

  Real fractalDiffNorm = fractalDiff.sumSq();
  if (fractalDiffNorm > 1e-8)
    final = false;
  cout << " Diff between CPU and CUDA fractal: " << fractalDiffNorm << endl;

  // if there's no distance field left, doesn't make sense to test scoring
  if (_distanceField.totalCells() <= 0)
    return final;

  // verify basic scoring function
  Real scoreCPU = computeScore();
  Real scoreGPU = computeScoreCUDA();
  Real scoreDiff = fabs(scoreCPU - scoreGPU);
  if (scoreDiff > 1e-8)
    final = false;
  cout << " CPU score: " << scoreCPU << endl;
  cout << " GPU score: " << scoreGPU << endl;
  cout << " Diff between CPU and CUDA scores: " << scoreDiff << endl;

  // verify banded scoring function
  Real bandedScoreCPU = computeBandedScore();
  Real bandedScoreGPU = computeBandedScoreCUDA();
  Real bandedScoreDiff = fabs(bandedScoreCPU - bandedScoreGPU);
  if (bandedScoreDiff > 1e-8)
    final = false;
  cout << " CPU banded score: " << bandedScoreCPU << endl;
  cout << " GPU banded score: " << bandedScoreGPU << endl;
  cout << " Diff between CPU and CUDA banded scores: " << bandedScoreDiff << endl;

  // verify centered gradients
  VECTOR gradientCPU = computeCenteredGradient();
  VECTOR gradientGPU = computeCenteredGradientCUDA();
  Real gradientDiff = (gradientCPU - gradientGPU).norm2();
  cout << " CPU gradient: " << gradientCPU.norm2() << endl;
  cout << " GPU gradient: " << gradientGPU.norm2() << endl;
  cout << " Diff between CPU and CUDA gradients: " << gradientDiff << endl;

  // verify centered banded gradients
  gradientCPU = computeBandedCenteredGradient();
  gradientGPU = computeBandedCenteredGradientCUDA();
  gradientDiff = (gradientCPU - gradientGPU).norm2();
  cout << " CPU banded gradient: " << gradientCPU.norm2() << endl;
  cout << " GPU banded gradient: " << gradientGPU.norm2() << endl;
  cout << " Diff between CPU and CUDA banded gradients: " << gradientDiff << endl;

  // verify the scale and translate score
  VECTOR scaleTranslateCPU = computeScaleTranslateGradient();
  VECTOR scaleTranslateGPU = computeScaleTranslateGradientCUDA();
  Real scaleTranslateDiff = (scaleTranslateCPU - scaleTranslateGPU).norm2();
  cout << " CPU scale/translate gradient: " << scaleTranslateCPU.norm2() << endl;
  cout << " GPU scale/translate gradient: " << scaleTranslateGPU.norm2() << endl;
  cout << " Diff between CPU and CUDA scale/translate gradients: " << scaleTranslateDiff << endl;
  if (scaleTranslateDiff > 1e-5)
    final = false;
    */

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::polynomialViewer(const OPTIMIZE_3D& optimize3D)
{
  optimize3D.write("temp.o3d");
  system("./bin/viewPolynomial3D temp.o3d &");
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::logPolynomialViewer(const OPTIMIZE_3D& optimize3D, const int expPow)
{
  string command("./bin/viewLogPolynomial3D temp.o3d ");

  char buffer[256];
  sprintf(buffer, "%i", expPow);
  command = command + string(buffer) + string(" &");

  optimize3D.write("temp.o3d");
  cout << " Executing " << command.c_str() << endl;
  system(command.c_str());
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::setLengthCenter(const Real& newLength, const VEC3F& newCenter)
{
  _lengths = newLength;
  _center = newCenter;

  _fractal.setLengths(_lengths);
  _fractal.setCenter(_center);
  _distanceField.setLengths(_lengths);
  _distanceField.setCenter(_center);
  _curvatureField.setLengths(_lengths);
  _curvatureField.setCenter(_center);
}

///////////////////////////////////////////////////////////////////////
// compute the derivative with respect to scaling and translation
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeScaleTranslateGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  //int totalRoots = _top.totalRoots();
  //int totalDOFs = totalRoots * 4;
  VECTOR gradient(4);

  const Real dx = _lengths[0] / _xRes;
  //const Real eps = dx;

  const VEC3F lengthsOriginal = _lengths;
  const VEC3F centerOriginal = _center;

  // get the center gradients
  for (int x = 0; x < 3; x++)
  {
    // do the left
    VEC3F center = centerOriginal;
    center[x] -= dx;
    _center = center;
    Real leftScore = computeScore();

    center = centerOriginal;
    center[x] += dx;
    _center = center;
    Real rightScore = computeScore();
    gradient[x] = (rightScore - leftScore) / (2.0 * dx);
  }

  // get the length gradient
  _lengths = lengthsOriginal + VEC3F(-dx, -dx, -dx);
  Real leftScore = computeScore();
  _lengths = lengthsOriginal + VEC3F(dx, dx, dx);
  Real rightScore = computeScore();
  gradient[3] = (rightScore - leftScore) / (2.0 * dx);

  _center = centerOriginal;
  _lengths = lengthsOriginal;

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative with respect to scaling and translation
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeLogScaleTranslateGradient()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  //int totalRoots = _top.totalRoots();
  //int totalDOFs = totalRoots * 4;
  VECTOR gradient(4);

  const Real dx = _lengths[0] / _xRes;
  //const Real eps = dx;

  const VEC3F lengthsOriginal = _lengths;
  const VEC3F centerOriginal = _center;

  // get the center gradients
  for (int x = 0; x < 3; x++)
  {
    // do the left
    VEC3F center = centerOriginal;
    center[x] -= dx;
    _center = center;
    Real leftScore = computeLogScore();

    center = centerOriginal;
    center[x] += dx;
    _center = center;
    Real rightScore = computeLogScore();
    gradient[x] = (rightScore - leftScore) / (2.0 * dx);
  }

  // get the length gradient
  _lengths = lengthsOriginal + VEC3F(-dx, -dx, -dx);
  Real leftScore = computeLogScore();
  _lengths = lengthsOriginal + VEC3F(dx, dx, dx);
  Real rightScore = computeLogScore();
  gradient[3] = (rightScore - leftScore) / (2.0 * dx);

  _center = centerOriginal;
  _lengths = lengthsOriginal;

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative with respect to scaling and translation
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::computeScaleTranslateGradientCUDA()
{
  // output a dot just to show some progress
  cout << "." << flush;
  TIMER functionTimer(__FUNCTION__);
  //cout << " Computing gradient ... " << flush;
  //int totalRoots = _top.totalRoots();
  //int totalDOFs = totalRoots * 4;
  VECTOR gradient(4);

  const Real dx = _lengths[0] / _xRes;
  //const Real eps = dx;

  const VEC3F lengthsOriginal = _lengths;
  const VEC3F centerOriginal = _center;

  // get the center gradients
  for (int x = 0; x < 3; x++)
  {
    // do the left
    VEC3F center = centerOriginal;
    center[x] -= dx;
    _center = center;
    Real leftScore = computeScoreCUDA();

    center = centerOriginal;
    center[x] += dx;
    _center = center;
    Real rightScore = computeScoreCUDA();
    gradient[x] = (rightScore - leftScore) / (2.0 * dx);
  }

  // get the length gradient
  _lengths = lengthsOriginal + VEC3F(-dx, -dx, -dx);
  Real leftScore = computeScoreCUDA();
  _lengths = lengthsOriginal + VEC3F(dx, dx, dx);
  Real rightScore = computeScoreCUDA();
  gradient[3] = (rightScore - leftScore) / (2.0 * dx);

  _center = centerOriginal;
  _lengths = lengthsOriginal;

  return gradient;
}

///////////////////////////////////////////////////////////////////////
// IO for the scale translate optimization
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::writeScaleTranslate(const string& path, int subdivisions, bool includeOrigin) const
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;

  char buffer[1024];
  sprintf(buffer, "scaleTranslate.sub.%i.origin.%i.opt", subdivisions, (int)includeOrigin);
  string filename = path + string(buffer);
  file = fopen(filename.c_str(), "wb");

  if (file == NULL)
  {
    cout << " Couldn't open the file " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  _center.write(file);
  _lengths.write(file);

  fclose(file);
  cout << " Wrote file: " << filename.c_str() << endl;
}

///////////////////////////////////////////////////////////////////////
// IO for the scale translate optimization
///////////////////////////////////////////////////////////////////////
bool OPTIMIZE_3D::readScaleTranslate(const string& path, int subdivisions, bool includeOrigin)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;

  char buffer[1024];
  sprintf(buffer, "scaleTranslate.sub.%i.origin.%i.opt", subdivisions, (int)includeOrigin);
  string filename = path + string(buffer);
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
    return false;
  cout << " Found file: " << filename.c_str() << endl;

  VEC3F center, lengths;
  center.read(file);
  lengths.read(file);

  setLengthCenter(lengths[0], center);

  fclose(file);

  return true;
}

///////////////////////////////////////////////////////////////////////
// compute an error field along the fractal borer
///////////////////////////////////////////////////////////////////////
FIELD_3D OPTIMIZE_3D::errorField()
{
  FIELD_3D final(_fractal);
  final.clear();

  int xRes = _fractal.xRes();
  int yRes = _fractal.yRes();
  int zRes = _fractal.zRes();

  for (int z = 1; z < zRes - 1; z++)
    for (int y = 1; y < yRes - 1; y++)
      for (int x = 1; x < xRes - 1; x++)
      {
        bool border = false;

        if (_fractal(x,y,z) * _fractal(x + 1,y,z) < 0) border = true;
        if (_fractal(x,y,z) * _fractal(x - 1,y,z) < 0) border = true;
        if (_fractal(x,y,z) * _fractal(x,y + 1,z) < 0) border = true;
        if (_fractal(x,y,z) * _fractal(x,y - 1,z) < 0) border = true;
        if (_fractal(x,y,z) * _fractal(x,y,z + 1) < 0) border = true;
        if (_fractal(x,y,z) * _fractal(x,y,z - 1) < 0) border = true;

        if (!border) continue;

        int cellsInside = fabs(_distanceField(x,y,z)) * _fractal.maxRes();

        final(x,y,z) = _curvatureField(x,y,z) * cellsInside;
        //final(x,y,z) = _curvatureField(x,y,z);
        //final(x,y,z) = cellsInside;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute log versions of the score
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogBandedScore(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal)
{
  computeLogBandedMap(polynomial, fractal);
  return computeBinaryBandedScore(fractal);
}

///////////////////////////////////////////////////////////////////////
// compute log versions of the score
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogScore(const POLYNOMIAL_4D& polynomial, FIELD_3D& fractal)
{
  computeLogMap(polynomial, fractal);
  //return computeBinaryScore(fractal);
  return computeInverseDistanceScore(fractal);
}
/*
{
  computeLogMap(polynomial, fractal);
  Real binaryScore = computeBinaryScore(fractal);
  Real rootDistanceScore = computeRootDistanceScore();

  //return binaryScore + rootDistanceScore;
  return 0.9 * binaryScore + 0.1 * rootDistanceScore;
}
*/

///////////////////////////////////////////////////////////////////////
// compute log versions of the score
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computePowerScore(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal)
{
  computeLogPowerRationalMap(top, bottom, fractal);
  //return computeBinaryScore(fractal);
  return computeInverseDistanceScore(fractal);
  //return computeBinaryBandedScore(fractal);
}

///////////////////////////////////////////////////////////////////////
// compute log versions of the score
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogRationalScore(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, FIELD_3D& fractal)
{
  computeLogRationalMap(top, bottom, fractal);
  return computeBinaryScore(fractal);
}

///////////////////////////////////////////////////////////////////////
// compute log versions of the score
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeLogRationalBottomScore(const POLYNOMIAL_4D& bottom, FIELD_3D& fractal)
{
  computeLogRationalMap(_top, bottom, fractal);
  return computeBinaryScore(fractal);
  //return computeBinaryBandedScore(fractal);
}

///////////////////////////////////////////////////////////////////////
// compute a score based on the distance field values
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeRootDistanceScore() const
{
  Real final = 0;
  int totalRoots = _top.totalRoots();
  //cout << " max distance: " << _distanceField.fieldMax() << endl;

  for (int x = 0; x < _top.totalRoots(); x++)
  {
    QUATERNION root = _top.roots()[x];
    VEC3F point(root[0], root[1], root[2]);

    Real distance = _distanceField(point);
    //cout << " distance: " << distance << endl;
    final += fabs(distance);
  }
  //cout << " final score: " << final << endl;
  //exit(0);
  return final / totalRoots;
}

///////////////////////////////////////////////////////////////////////
// scale the entire optimization problem
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::scaleEverything(const Real& scale)
{
  _lengths *= scale;
  _distanceField.scaleLengths(scale);
  _inverseDistanceField.scaleLengths(scale);
  _weightField.scaleLengths(scale);
  _curvatureField.scaleLengths(scale);
  _fractal.scaleLengths(scale);
  _auxiliary.scaleLengths(scale);

  _center *= scale;
  _distanceField.center() *= scale;
  _inverseDistanceField.center() *= scale;
  _weightField.center() *= scale;
  _curvatureField.center() *= scale;
  _fractal.center() *= scale;
  _auxiliary.center() *= scale;

  _distanceField *= scale;

  _top *= scale;
  _bottom *= scale;
}

///////////////////////////////////////////////////////////////////////
// translate the entire optimization problem
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::translateEverything(const VEC3F& translation)
{
  _center += translation;
  _distanceField.center() += translation;
  _inverseDistanceField.center() += translation;
  _weightField.center() += translation;
  _curvatureField.center() += translation;
  _fractal.center() += translation;
  _auxiliary.center() += translation;

  _top += translation;
  _bottom += translation;
}

///////////////////////////////////////////////////////////////////////
// translate the entire optimization problem
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::translateEverythingExceptFirst(const VEC3F& translation)
{
  _center += translation;
  _distanceField.center() += translation;
  _inverseDistanceField.center() += translation;
  _weightField.center() += translation;
  _curvatureField.center() += translation;
  _fractal.center() += translation;
  _auxiliary.center() += translation;

  _top.translateExceptFirst(translation);
  _bottom.translateExceptFirst(translation);
}

///////////////////////////////////////////////////////////////////////
// load up the kitty as the polynomial
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::loadKitty()
{
  int kittySize = 298;
  double xKitty[] = {-1.50616918,-1.504,-1.502385206,-1.501256117,-1.500948054,-1.5,-1.5,-1.499858653,-1.499002019,-1.498501499,-1.498401111,-1.496983317,-1.496393613,-1.496,-1.493457736,-1.493,-1.489936225,-1.488714563,-1.487846747,-1.486641716,-1.479007292,-1.478,-1.474324366,-1.46868558,-1.464143287,-1.463,-1.457697501,-1.455,-1.449126895,-1.442102386,-1.434120908,-1.433,-1.42056555,-1.419666293,-1.418353777,-1.418159757,-1.418,-1.413380118,-1.412526182,-1.411315727,-1.411309131,-1.411,-1.411,-1.410716215,-1.406027838,-1.405315915,-1.4,-1.398744055,-1.394642464,-1.393267299,-1.392226225,-1.386128432,-1.380283606,-1.379950722,-1.373001408,-1.37226157,-1.37098414,-1.362,-1.356869829,-1.349226246,-1.348,-1.335444052,-1.326770339,-1.306694361,-1.299,-1.295623047,-1.280506499,-1.27308128,-1.251316199,-1.250606427,-1.23479669,-1.220946985,-1.220061242,-1.21941706,-1.217954764,-1.215247398,-1.212344887,-1.210888174,-1.210778099,-1.208944221,-1.206843283,-1.201281367,-1.200669319,-1.187083531,-1.18392272,-1.175919998,-1.17386755,-1.173163747,-1.167377104,-1.161934884,-1.147249669,-1.144362353,-1.140190678,-1.131696954,-1.127188804,-1.119626942,-1.118928941,-1.102533024,-1.088373904,-1.088225147,-1.085416351,-1.083607031,-1.082251041,-1.081807492,-1.07699944,-1.076873061,-1.070922763,-1.068872173,-1.06717507,-1.066798596,-1.054008426,-1.052655827,-1.041011501,-1.034441702,-1.009522298,-0.994214746,-0.966863385,-0.935024213,-0.916579949,-0.891811022,-0.858547144,-0.839012539,-0.829532436,-0.804432258,-0.780227507,-0.75598741,-0.740712302,-0.738935744,-0.732068337,-0.718640386,-0.712404039,-0.691755024,-0.68557809,-0.676708173,-0.663010584,-0.661441298,-0.658049369,-0.656155959,-0.647664437,-0.646796704,-0.646350208,-0.64538565,-0.642159424,-0.640068441,-0.63980835,-0.639594589,-0.638493391,-0.635102156,-0.574234779,-0.57064013,-0.549223701,-0.522525227,-0.503787447,-0.490354262,-0.476483561,-0.433959538,-0.40312035,-0.39419133,-0.363684937,-0.340075715,-0.326167816,-0.321732677,-0.292339793,-0.26191209,-0.256080285,-0.230961391,-0.205009984,-0.198192242,-0.169152552,-0.1610428,-0.135588033,-0.119949056,-0.102829361,-0.080250829,-0.072356715,-0.045036222,-0.031394466,-0.016415547,0.008377875,0.012139371,0.042828359,0.046705186,0.072980429,0.088904311,0.101702142,0.12843301,0.131320878,0.156637573,0.175739065,0.20486406,0.206531781,0.231934013,0.235970752,0.253081037,0.262591547,0.284941667,0.299207935,0.304514774,0.321460715,0.336067013,0.35046414,0.378436439,0.397752337,0.425754539,0.448773056,0.449733874,0.468461557,0.483302322,0.494960127,0.503980641,0.507251668,0.507279764,0.531655458,0.540074505,0.564125979,0.580817967,0.594538558,0.603427761,0.605039211,0.606992828,0.625843335,0.649327432,0.667696918,0.681244077,0.692267196,0.743300992,0.786838304,0.804149316,0.839805196,0.8849345,0.899316313,0.9484047,0.982164213,0.992731179,1.041667946,1.086183407,1.106652611,1.115985589,1.119098319,1.128066712,1.128394176,1.129193303,1.135962746,1.142156638,1.151419152,
               1.157965895,1.16099025,1.167600825,1.173359604,1.187323043,1.1992497,1.201621783,1.21102322,1.220854152,1.222849097,1.228886911,1.239658376,1.255,1.255038556,1.258,1.262,1.263388371,1.275395563,1.291440868,1.292969309,1.3,1.307,1.310213444,1.314870291,1.332148621,1.334730671,1.336,1.341900163,1.348141491,1.349,1.350765783,1.359280131,1.363368472,1.375819991,1.377,1.379097225,1.384776625,1.389,1.393870645,1.398097323,1.407,1.412516789,1.415,1.418,1.422,1.422,1.423568098,1.425408062,1.429,1.429,1.432183632,1.433,1.433074507};
  double yKitty[] = {-1.378339323,-1.41,-1.424515502,-1.120969104,-1.46073756,-1.358,-1.518,-1.092893561,-1.347341472,-1.498501499,-1.149963392,-1.318565722,-1.294742664,-1.486,-1.073440366,-1.548,-1.269607036,-1.208125563,-1.242377081,-1.181859124,-1.552574494,-1.596,-1.048335713,-1.571223734,-1.590308825,-1.612,-1.609360906,-1.596,-1.012617787,-1.647432199,-1.635165896,-1.638,-0.98251226,-0.553317743,-0.519429292,-0.593098749,-1.668,-0.624896307,-1.672034949,-0.489404874,-0.458813112,-1.66,-1.682,-0.659681623,-0.422454966,-0.701424503,-1.674,-0.386005635,-0.743531773,-0.355639835,-1.668467834,-0.332562453,-0.791627039,-0.844636294,-0.304598632,-1.664710757,-0.921264964,-1.656,-0.280268326,-1.65462543,-1.656,-0.249799717,-1.644510627,-0.221930165,-1.616,-1.629278481,-0.188753248,-1.600763756,-1.585045831,-0.154010356,-1.561829521,0.967375021,0.946948803,0.982728332,0.923885098,0.93474542,0.912616696,0.894380285,-1.539363074,-0.11303137,1.007831571,1.019269038,0.879186694,1.024407604,0.853418356,-1.52344911,1.036534856,-0.064421622,0.824551573,1.039583302,1.046754779,0.796162656,-1.498198735,-0.013663197,1.045634074,0.763145557,1.051722461,1.037588252,0.722359476,0.069760563,1.029113407,-1.48967689,0.464052977,0.416596674,0.515685963,0.364486741,0.309409903,1.008750363,0.568055755,0.144480418,0.640313921,0.241440173,0.988789142,-1.485576256,0.960340442,-1.482032873,0.926515145,-1.486580552,0.877919283,-1.509928317,-1.522571744,0.744065579,-1.534330778,-1.546089437,-1.555408265,-1.562173188,-1.558871734,-1.051297782,-1.572442675,-1.568213182,-1.173519175,-1.584856612,-1.258833695,-1.587756961,-1.313749878,-0.802599253,-1.564908645,-1.549112991,-1.51529038,-1.505760631,-1.420104034,-1.533425233,-1.359689673,-1.392031487,-1.48996465,-1.449970978,-1.47616568,-1.474743302,0.940134996,1.006146256,1.043733686,1.069123432,1.093593277,-0.787608221,1.117203974,1.149012218,-0.765566843,1.17237499,1.194779106,1.201306421,-0.73406578,1.209172121,1.219699297,1.233025976,-0.724986048,1.242982528,-0.713493527,1.251205547,1.25488212,-0.69586815,1.258192889,-0.67573395,1.261729614,-0.650160459,1.265974532,1.259954551,-0.630948958,1.260290386,-0.617043213,1.258272145,-0.6000796,1.250342189,-0.57985854,1.242075929,-0.556455244,-0.529045208,1.228715499,-0.494338321,1.205916516,-0.472953941,1.196900843,1.187228996,-0.456201327,1.174469223,-0.435854499,-0.414024283,1.141586208,-0.389609373,-0.36098617,-0.327766341,-0.281719251,1.076104074,-0.243399061,-0.220832611,1.007913446,-0.195163751,-0.169838006,-0.142316965,-0.113362113,-0.079664955,0.013561043,-0.042849244,0.965688847,0.067441427,0.105272779,0.142596385,0.181330075,0.222725593,0.326141647,0.268615989,0.403913074,0.464321153,0.524876619,0.979852108,0.604431749,0.666247448,0.725020746,1.024336423,0.784551958,1.062813985,0.826465973,0.8596539,1.096671773,0.889932696,0.915104169,0.941542624,1.15886128,1.476938963,1.495231624,1.493579012,0.963943756,1.498160881,1.495801877,1.499969195,1.497558587,1.44104954,1.493894957,
               0.98348812,1.498856684,1.497853685,1.500599625,0.997557461,1.498772733,1.313329021,1.496807143,1.012173613,1.495640598,1.016,1.495878866,1.014,1.024,1.025896077,1.490205776,1.487388081,1.039329001,1.456,1.046,1.473504141,1.047750458,1.464057455,1.057346383,1.448,1.066663706,1.444437311,1.39,1.078819574,1.093918612,1.42903381,1.406135762,1.404,1.112276442,1.134309893,1.37,1.371563615,1.151502918,1.166,1.329722732,1.348,1.322,1.188,1.3,1.295162996,1.193793656,1.232,1.274,1.252406579,1.248,1.217907006};
  vector<QUATERNION> roots;
  for (int x = 0; x < kittySize; x++)
  {
    QUATERNION root(xKitty[x], yKitty[x],0,0);
    roots.push_back(root);
  }
  POLYNOMIAL_4D top(roots);
  _top = top;
}

///////////////////////////////////////////////////////////////////////
// compute the error on the outside 
// (assuming the convex hull has been built already)
///////////////////////////////////////////////////////////////////////
FIELD_3D OPTIMIZE_3D::computeExteriorErrorField()
{
  // set up dimensions for an error field run
  _fractal = FIELD_3D(_distanceField);
  _auxiliary = FIELD_3D(_distanceField);
  _fractal = 0;
  _auxiliary = 0;
  _maxIterations = 20;
  
  // compute the convex hull (presumably) so it can be diffed
  // against the SDF
  computeLogRationalMap(true);

  FIELD_3D final(_distanceField);
  final = 0;
  for (int x = 0; x < _distanceField.totalCells(); x++)
  {
    if (_distanceField[x] < 0 && _fractal[x] < 0) continue;

    if (_distanceField[x] > 0 && _fractal[x] < 0)
      final[x] = _distanceField[x];
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get a vector of all the root exponents
///////////////////////////////////////////////////////////////////////
VECTOR OPTIMIZE_3D::getPowers() const
{
  const int topPowers = _top.totalRoots();
  const int bottomPowers = _bottom.totalRoots();
  const int totalPowers = topPowers + bottomPowers;
  VECTOR final(totalPowers);

  for (int x = 0; x < topPowers; x++)
    final[x] = _top.powers()[x];
  
  for (int x = 0; x < bottomPowers; x++)
    final[topPowers + x] = _bottom.powers()[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// set of all the root exponents
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::setPowers(const VECTOR& powers)
{
  const int topPowers = _top.totalRoots();
  const int bottomPowers = _bottom.totalRoots();
  //const int totalPowers = topPowers + bottomPowers;

  for (int x = 0; x < topPowers; x++)
    _top.changePower(x, powers[x]);
  
  for (int x = 0; x < bottomPowers; x++)
    _bottom.changePower(x, powers[topPowers + x]);
}

///////////////////////////////////////////////////////////////////////
// compute the fitting score for a given threshold
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeMaximumScaledPowerRationalScore()
{
  Real tanhScore = 0;
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);

    Real weight = _weightField(cellPosition);
    tanhScore += fabs(weight);
  }

  return tanhScore;
}

///////////////////////////////////////////////////////////////////////
// compute the fitting score for a given threshold
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::weightScore(const FIELD_3D& logFractal, const FIELD_3D& weightField)
{
  Real tanhScore = 0;
  for (int x = 0; x < logFractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = logFractal.cellCenter(x);

    Real weight = weightField(cellPosition);
    Real tanhVersionBefore = -tanh(_tanhAlpha * logFractal[x]);
    Real tanhVersionAfter = tanhVersionBefore * weight;
    tanhScore += tanhVersionAfter;
  }

  return tanhScore;
}

///////////////////////////////////////////////////////////////////////
// compute the fitting score for a given threshold
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::inverseDistanceScore(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField, const Real& inverseSum, const Real& threshold)
{
  Real tanhScore = 0;
  //Real scaleSize = inverseDistanceField.lengths().maxElement() /  logFractal.lengths().maxElement();

  for (int x = 0; x < logFractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = logFractal.cellCenter(x);

    /*
    // make sure the dims match up
    cellPosition -= logFractal.center();
    cellPosition *= scaleSize;
    cellPosition += inverseDistanceField.center();
    */

    Real inverseDistance = inverseDistanceField(cellPosition);

    /*
    if (x == 21 + 9 * logFractal.xRes() + 21 * logFractal.slabSize())
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Inverse cell center: " << cellPosition << endl;
      cout << " Inverse distance: " << inverseDistance << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    }
    */

    //Real fractalThresholded = (logFractal[x] > threshold) ? -1 : 1;
    //Real tanhVersionBefore = -tanh(_tanhAlpha * (logFractal[x] - _tanhThreshold));
    Real tanhVersionBefore = -tanh(_tanhAlpha * (logFractal[x] - threshold));

    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    tanhScore += tanhVersionAfter;
  }
  //tanhScore *= 1.0 / inverseSum;

  return tanhScore;
}

///////////////////////////////////////////////////////////////////////
// compute the fitting score for a given threshold
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::weightScoreConformalGradient(const FIELD_3D& logFractal, const FIELD_3D& weightField)
{
  Real gradientScore = 0;
  //Real scaleSize = weightField.lengths().maxElement() /  logFractal.lengths().maxElement();

  FIELD_3D gradientField(logFractal);
  gradientField = 0;

  /*
  // DEBUG
  const int xDebug = 10;
  const int yDebug = 10;
  const int zDebug = 10;
  const int index = xDebug + yDebug * _fractal.xRes() + zDebug * _fractal.slabSize();
  */

  for (int x = 0; x < logFractal.totalCells(); x++)
  {
    // ADDING OVERFLOW DETECTION
    if (_overflowed[x] > 0.5) continue;

    // do a more robust, resolution-independent version here
    VEC3F cellPosition = logFractal.cellCenter(x);
    Real weight = weightField(cellPosition);

    //Real fractalThresholded = (logFractal[x] > threshold) ? -1 : 1;
    Real tanhVersionBefore = tanh(_tanhAlpha * (logFractal[x] - _tanhThreshold));
    Real S = tanhVersionBefore;

    Real currentScore = -_tanhAlpha * weight * ((Real)1.0 - S * S);

    //gradientScore += _tanhAlpha * inverseDistance * (1 - S * S);
    gradientScore += currentScore;
    gradientField[x] = currentScore;

    /*
    if (x == index)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " F: " << logFractal[x] << endl;
      cout << " tanh: " << tanhVersionBefore << endl;
      cout << " weight: " << weight << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    }
    */
  }

  // this actually seems to change the path, probably because the score is very small sometimes,
  // so the significant digits in the 1.0 dominate
  //
  // not needed, since we normalized ...
  //gradientScore *= 1.0 / _inverseFieldSum;

  return gradientScore;
}

///////////////////////////////////////////////////////////////////////
// compute the fitting score for a given threshold
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::inverseDistanceScoreConformalGradient(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField, const Real& inverseSum, const Real& threshold)
{
  Real gradientScore = 0;
  //Real scaleSize = inverseDistanceField.lengths().maxElement() /  logFractal.lengths().maxElement();
  //const Real alpha = 10.0;

  FIELD_3D gradientField(logFractal);
  gradientField = 0;

  for (int x = 0; x < logFractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = logFractal.cellCenter(x);
    Real inverseDistance = inverseDistanceField(cellPosition);

    //Real fractalThresholded = (logFractal[x] > threshold) ? -1 : 1;
    Real tanhVersionBefore = tanh(_tanhAlpha * (logFractal[x] - _tanhThreshold));
    Real S = tanhVersionBefore;

    Real currentScore = -_tanhAlpha * inverseDistance * ((Real)1.0 - S * S);

    //gradientScore += _tanhAlpha * inverseDistance * (1 - S * S);
    gradientScore += currentScore;
    gradientField[x] = currentScore;
  }

  // not needed, since we normalized ...
  gradientScore *= 1.0 / inverseSum;

  return gradientScore;
}

///////////////////////////////////////////////////////////////////////
// compute the fitting score for a given threshold
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::inverseDistanceScoreDebug(const FIELD_3D& logFractal, const FIELD_3D& inverseDistanceField, const Real& inverseSum, const Real& threshold)
{
  //assert(logFractal.totalCells() == inverseDistanceField.totalCells());

  FIELD_3D product(logFractal);
  product = 0;

  Real score = 0;
  for (int x = 0; x < logFractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = logFractal.cellCenter(x);
    Real inverseDistance = inverseDistanceField(cellPosition);
    Real fractalThresholded = (logFractal[x] > threshold) ? -1 : 1;

    Real currentScore = fractalThresholded * inverseDistance;
    product[x] = currentScore;
    score += currentScore;
  }
  score *= 1.0 / inverseSum;

  FIELDVIEW3D(product);

  return score;
}

///////////////////////////////////////////////////////////////////////
// take a first shot at estimating the conformal radius
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::estimateConformalRadius()
{
  TIMER functionTimer(__FUNCTION__);
  
  // get all the inputs lined up
  VEC3F center = _fractal.center();
  VEC3F lengths = _fractal.lengths();

  int fractalRes = _fractal.xRes();
  _fractal   = FIELD_3D(fractalRes,fractalRes,fractalRes, center, lengths);
  _auxiliary = FIELD_3D(fractalRes,fractalRes,fractalRes, center, lengths);
 
  //Real oldScaling = _expScaling;
  Real oldMaxIterations = _maxIterations;

  // prep to look at the raw polynomial
  _expScaling = 1;
  _maxIterations = 1;

  // compute the log map -- not the power version!
  //computeLogMap(true);
  cout << " USING LOG POWER RATIONAL MAP " << flush;
  computeLogPowerRationalMap(true);

  _fractal.clampInfsToMean();
  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Why are Nans appearing here when the computation should be exactly the same as before?" << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  exit(0);
  */

  const FIELD_3D& logFractal = _fractal;
  const FIELD_3D& distanceField = _distanceField;

  FIELD_3D inverseDistanceField = distanceField.invertedDistance();
  Real inverseSum = inverseDistanceField.absSum();

  Real minScore = logFractal.fieldMin();
  Real maxScore = logFractal.fieldMax();
  Real minThreshold = minScore;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " fractal res: " << logFractal.xRes() << endl;
  cout << " distance res: " << _distanceField.xRes() << endl;
  cout << " score min max: " << minScore << " " << maxScore << endl;

  // do multiple scans to converge on a good value
  Real minScoreFound = FLT_MAX;
  for (int i = 0; i < 4; i++)
  {
    Real increment = (maxScore - minScore) / 100;
    cout << " increment " << increment << endl;
    minScoreFound = FLT_MAX;
  
    // do a linear scan to see where it looks there's a fit
    for (Real x = minScore; x < maxScore; x += increment)
    {
      Real currentScore = inverseDistanceScore(logFractal, inverseDistanceField, inverseSum, x);
      if (currentScore < minScoreFound)
      {
        minScoreFound = currentScore;
        minThreshold = x;
      }
    }

    minScore = minThreshold - 5 * increment;
    maxScore = minThreshold + 5 * increment;
    cout << " Power scan returned " << (double)minThreshold << " as the conformal radius with score " << (double)minScoreFound << endl;
  }

  _expScaling = exp(-minThreshold);
  cout << " Final radius: " << -minThreshold << " exponentiated: " << _expScaling << endl;
  _maxIterations = oldMaxIterations;

  // DEBUG: look at what the score field looks like
  //inverseDistanceScoreDebug(logFractal, inverseDistanceField, inverseSum, minThreshold);

  return minScoreFound;
}

///////////////////////////////////////////////////////////////////////
// plot one of the scoring functions with respect to one variable
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::plotExpEnergy1D()
{
  int totalSamples = 100;
  //int totalSamples = 20;
  //Real start = 1.0;
  //Real end = 10.0;
  Real start = 1.0;
  Real end = 200.0;

  Real dx = (end - start) / totalSamples;
  vector<Real> xs;
  vector<Real> ys;
  vector<Real> gradients;

  cout << " Plotting ... " << flush;
  for (int x = 0; x < totalSamples; x++)
  {
    Real current = start + x * dx;

    // var to probe
    //_bottom.powersMutable()[0] = current;
    _expScaling = exp(-current);

    // score to probe
    Real energy = computeLogPowerRationalScore();
    Real gradient = computeLogPowerRationalScoreExpGradient();

    xs.push_back(current);
    ys.push_back(energy);
    gradients.push_back(gradient);
    cout << x << " " << flush;
  }
  cout << endl;

  /*
  VECTOR xsMatlab(xs);
  VECTOR ysMatlab(ys);
  VECTOR gradientsMatlab(gradients);
  xsMatlab.writeMatlab("xExp1D.m", "x");
  ysMatlab.writeMatlab("yExp1D.m", "y");
  gradientsMatlab.writeMatlab("gradientExp1D.m", "gradient");
  */
}

///////////////////////////////////////////////////////////////////////
// plot two variables with respect to a scoring function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::plotEnergy2D(OPTIMIZE_3D& input)
{
  int totalSamples = 10;
  Real startX = 1.0;
  Real endX = 20.0;

  //Real startY = 1.0;
  //Real endY = 200.0;
  Real startY = 90.0;
  Real endY = 130.0;

  Real dx = (endX - startX) / totalSamples;
  Real dy = (endY - startY) / totalSamples;

  FIELD_2D result(totalSamples, totalSamples);
  vector<Real> xs;
  vector<Real> ys;
  for (int y = 0; y < totalSamples; y++)
    ys.push_back(startY + dy * y);
  
  for (int x = 0; x < totalSamples; x++)
    xs.push_back(startX + dx * x);

  cout << " Plotting ... " << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < totalSamples; y++)
  {
    OPTIMIZE_3D optimize = input;
    for (int x = 0; x < totalSamples; x++)
    {
      Real currentX = startX + x * dx;
      Real currentY = startY + y * dy;

      // vars to probe
      optimize._bottom.powersMutable()[0] = currentX;
      optimize._expScaling = exp(-currentY);

      // score to probe
      Real energy = optimize.computeLogPowerRationalScore();

      result(x,y) = energy;
      //cout << " currentX: " << currentX << endl;
      //cout << " currentY: " << currentY << endl;
      //cout << " energy: " << result(x,y) << flush;
    }
    cout << y << " " << flush;
  }
  cout << endl;

  /*
  result.writeMatlab("energies.m", "energy");
  
  VECTOR xsMatlab(xs);
  VECTOR ysMatlab(ys);
  xsMatlab.writeMatlab("xs.m", "x");
  ysMatlab.writeMatlab("ys.m", "y");
  */
}

///////////////////////////////////////////////////////////////////////
// draw the output of a score computation
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::drawInverseDistanceScore()
{
  const FIELD_3D logFractal = _fractal;
  const FIELD_3D inverseDistanceField = _inverseDistanceField;
  const Real inverseSum = _inverseFieldSum;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " inverse dims: " << _inverseDistanceField.center() << " " << _inverseDistanceField.lengths() << endl;
  cout << " distance dims: " << _distanceField.center() << " " << _distanceField.lengths() << endl;
  cout << " fractal dims: " << _fractal.center() << " " << _fractal.lengths() << endl;

  //Real scaleSize = _inverseDistanceField.lengths().maxElement() /  _fractal.lengths().maxElement();

  FIELD_3D scoreField(_fractal);
  FIELD_3D tanField(_fractal);
  FIELD_3D inverseField(_fractal);
  scoreField = 0;
  tanField = 0;
  inverseField = 0;

  Real tanhScore = 0;
  //Real threshold = 0;
  for (int x = 0; x < logFractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = logFractal.cellCenter(x);

    //cellPosition -= _fractal.center();
    //cellPosition *= scaleSize;
    //cellPosition += _inverseDistanceField.center();
    Real inverseDistance = inverseDistanceField(cellPosition);

    //Real fractalThresholded = (logFractal[x] > threshold) ? -1 : 1;
    Real tanhVersionBefore = -tanh(_tanhAlpha * (logFractal[x] - _tanhThreshold));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    tanhScore += tanhVersionAfter;

    scoreField[x] = tanhVersionAfter;
    tanField[x] = tanhVersionBefore;
    inverseField[x] = inverseDistance;
  }
  tanhScore *= 1.0 / inverseSum;

  scoreField *= 1.0 / scoreField.max();

  FIELDVIEW3D(scoreField);
  FIELDVIEW3D(tanField);
  FIELDVIEW3D(inverseField);

  cout << " score: " << tanhScore << endl;
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalBottomScoreDerivative(int whichRoot)
{
  //const Real alpha = 10.0;
  //const Real threshold = 0;
  //const Real threshold = 0.01;
  
  int xDebug = 21;
  int yDebug = 9;
  int zDebug = 21;
  const QUATERNION iterate(1.36, -0.56, 1.36, 0);

  QUATERNION topEval    = iterate * _top.evaluatePowerFactored(iterate);
  QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
  QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  cout << " log end R (experm): " << F << endl;

  // get dFdt
  QUATERNION bottomDerivative = _bottom.inversePowerDerivative(iterate, whichRoot);
  QUATERNION Rderivative = topEval / (iterate * bottomDerivative);
  // NOTE: R.magnitude instead of R.dot(R) is a fairly pernicious bug
  //Real dFdt = (1.0 / R.magnitude()) * (R.dot(Rderivative));
  Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));

  // get dSdF
  Real inverseDistance = _inverseDistanceField(xDebug, yDebug, zDebug);
  Real tanhTerm = tanh(_tanhAlpha * inverseDistance * (F - _tanhThreshold));
  Real S = -tanhTerm;
  Real dSdF = -_tanhAlpha * inverseDistance * (1 - tanhTerm * tanhTerm); 
  cout << " final S: " << S << endl; 

  // get the final dSdt
  Real dSdt = dSdF * dFdt;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of total score, BOTTOM root " << whichRoot << " derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " top: " << topEval << endl;
  cout << " bottom: " << bottomEval << endl;
  cout << " iterate: " << iterate << endl;
  cout << " R: " << R << endl;
  cout << " ||R||: " << R.magnitude() << endl;
  cout << " dSdt (analytical solution): " << dSdt << endl;
  
  Real dx = 0.1;
  Real originalPower = _bottom.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] += dx;
    QUATERNION topEvalNew = iterate * _top.evaluatePowerFactored(iterate);
    QUATERNION bottomEvalNew = iterate * _bottom.evaluatePowerFactored(iterate);
    _bottom.powersMutable()[whichRoot] = originalPower;
    
    QUATERNION Rnew = (topEvalNew / bottomEvalNew);
    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);

    /*
    Real tanhVersionBefore = -tanh(_tanhAlpha * (Fnew - _tanhThreshold));
    Real Snew = inverseDistance * tanhVersionBefore;
    */
    // DEBUG: try a simpler S
    Real Snew = -tanh(_tanhAlpha * inverseDistance * (Fnew - _tanhThreshold));

    Real numericalDt = (Snew - S) / dx;

    Real diff = (numericalDt - dSdt) / dSdt;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical dt: " << numericalDt << "\t Snew: " << Snew << "\t S: " << S << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalTopSummedScoreDerivative(int whichRoot)
{
  // compute the summed score via function
  Real groundScore = computeLogPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real topScoreGradient = computeTopPowerGradient(whichRoot);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of summed score, TOP root " << whichRoot << " derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " Top gradient: " << topScoreGradient << endl;
  
  Real dx = 0.1;
  Real originalPower = _top.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] += dx;
    Real newScore = computeLogPowerRationalScore();
    _top.powersMutable()[whichRoot] = originalPower;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - topScoreGradient) / topScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalTopScaledPowerDerivative()
{
  // compute the summed score via function
  Real groundScore = computeLogScaledPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluateScaledPowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluateScaledPowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real topScoreGradient = computeTopScaledPowerGradient();

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of SCALED summed score, TOP root GLOBAL SCALAR derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " Top gradient: " << topScoreGradient << endl;
  
  Real dx = 0.1;
  Real originalPower = _top.powerScalar();
  for (int x = 0; x < 10; x++)
  {
    _top.powerScalar() += dx;
    Real newScore = computeLogScaledPowerRationalScore();
    _top.powerScalar() = originalPower;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - topScoreGradient) / topScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalTopScaledPowerDerivative(int whichRoot)
{
  // compute the summed score via function
  Real groundScore = computeLogScaledPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluateScaledPowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluateScaledPowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real topScoreGradient = computeTopScaledPowerGradient(whichRoot);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of SCALED summed score, TOP root " << whichRoot << " derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " Top gradient: " << topScoreGradient << endl;
  
  Real dx = 0.1;
  Real originalPower = _top.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] += dx;
    Real newScore = computeLogScaledPowerRationalScore();
    _top.powersMutable()[whichRoot] = originalPower;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - topScoreGradient) / topScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalBottomScaledPowerDerivative(int whichRoot)
{
  // compute the summed score via function
  Real groundScore = computeLogScaledPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluateScaledPowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluateScaledPowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real bottomScoreGradient = computeBottomScaledPowerGradient(whichRoot);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of SCALED summed score, BOTTOM root " << whichRoot << " derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " Bottom gradient: " << bottomScoreGradient << endl;
  
  Real dx = 0.1;
  Real originalPower = _bottom.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] += dx;
    Real newScore = computeLogScaledPowerRationalScore();
    _bottom.powersMutable()[whichRoot] = originalPower;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - bottomScoreGradient) / bottomScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalBottomScaledPowerDerivative()
{
  // compute the summed score via function
  Real groundScore = computeLogScaledPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluateScaledPowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluateScaledPowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real bottomScoreGradient = computeBottomScaledPowerGradient();

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of SCALED summed score, BOTTOM root GLOBAL SCALAR derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " Bottom gradient: " << bottomScoreGradient << endl;
  
  Real dx = 0.1;
  Real originalPower = _bottom.powerScalar();
  for (int x = 0; x < 10; x++)
  {
    _bottom.powerScalar() += dx;
    Real newScore = computeLogScaledPowerRationalScore();
    _bottom.powerScalar() = originalPower;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - bottomScoreGradient) / bottomScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalBottomSummedScoreDerivative(int whichRoot)
{
  // compute the summed score via function
  Real groundScore = computeLogPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real bottomScoreGradient = computeBottomPowerGradient(whichRoot);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of summed score, BOTTOM root " << whichRoot << " derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " Bottom gradient: " << bottomScoreGradient << endl;
  
  Real dx = 0.1;
  Real originalPower = _bottom.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] += dx;
    Real newScore = computeLogPowerRationalScore();
    _bottom.powersMutable()[whichRoot] = originalPower;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - bottomScoreGradient) / bottomScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of all bottom powers at once
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalBulkBottomDerivative()
{
  // compute the summed score via function
  Real groundScore = computeLogPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real bottomScoreGradient = computeBottomBulkPowerGradient();
  //for (int x = 0; x < _top.totalRoots(); x++)
  //  topScoreGradient += computeTopPowerGradient(x);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of bulk power, BOTTOM root " << endl;
  cout << " ========================================================== " << endl;

  cout << " Bottom gradient: " << bottomScoreGradient << endl;
  
  Real dx = 0.1;
  for (int x = 0; x < 15; x++)
  {
    vector<Real> originalPowers = _bottom.powers();
    for (int y = 0; y < _bottom.totalRoots(); y++)
      _bottom.powersMutable()[y] += dx;
    Real newScore = computeLogPowerRationalScore();
    _bottom.powersMutable() = originalPowers;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - bottomScoreGradient) / bottomScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of all top powers at once
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalBulkTopDerivative()
{
  // compute the summed score via function
  Real groundScore = computeLogPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif
  Real score = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
        iterate = (topEval / bottomEval);
        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        // get the final fractal value
        Real fractalFinal = log(magnitude);

        // compute its contribution to the score
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real inverseDistance = _inverseDistanceField(cellPosition);
        Real tanhVersionBefore = -tanh(_tanhAlpha * (fractalFinal - _tanhThreshold));
        Real cellScore = tanhVersionBefore * inverseDistance;

        score += cellScore;
      }

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Total score, ground truth:     " << groundScore << endl;
  cout << " Total score, locally computed: " << score << endl;

  Real topScoreGradient = computeTopBulkPowerGradient();
  //for (int x = 0; x < _top.totalRoots(); x++)
  //  topScoreGradient += computeTopPowerGradient(x);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of bulk power, TOP root " << endl;
  cout << " ========================================================== " << endl;

  cout << " Top gradient: " << topScoreGradient << endl;
  
  Real dx = 0.1;
  for (int x = 0; x < 15; x++)
  {
    vector<Real> originalPowers = _top.powers();
    for (int y = 0; y < _top.totalRoots(); y++)
      _top.powersMutable()[y] += dx;
    Real newScore = computeLogPowerRationalScore();
    _top.powersMutable() = originalPowers;
    
    Real numericalGradient = (newScore - groundScore) / dx;

    Real diff = (numericalGradient - topScoreGradient) / topScoreGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalTopScoreDerivative(int whichRoot)
{
  //const Real alpha = 10.0;
  //const Real threshold = 0;
  //const Real threshold = 0.01;
  
  int xDebug = 21;
  int yDebug = 9;
  int zDebug = 21;
  const QUATERNION iterate(1.36, -0.56, 1.36, 0);

  QUATERNION topEval    = iterate * _top.evaluatePowerFactored(iterate);
  QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
  QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  cout << " log end R (experm): " << F << endl;

  // get dFdt
  QUATERNION topDerivative = _top.powerDerivative(iterate, whichRoot);
  QUATERNION Rderivative = (iterate * topDerivative) / bottomEval;
  // NOTE: R.magnitude instead of R.dot(R) is a fairly pernicious bug
  //Real dFdt = (1.0 / R.magnitude()) * (R.dot(Rderivative));
  Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));

  // get dSdF
  Real inverseDistance = _inverseDistanceField(xDebug, yDebug, zDebug);
  Real tanhTerm = tanh(_tanhAlpha * inverseDistance * (F - _tanhThreshold));
  Real S = -tanhTerm;
  Real dSdF = -_tanhAlpha * inverseDistance * (1 - tanhTerm * tanhTerm); 
  cout << " inverse distance: " << inverseDistance << endl;
  cout << " final S: " << S << endl; 

  // get the final dSdt
  Real dSdt = dSdF * dFdt;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of total score, TOP root " << whichRoot << " derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " top: " << topEval << endl;
  cout << " bottom: " << bottomEval << endl;
  cout << " iterate: " << iterate << endl;
  cout << " R: " << R << endl;
  cout << " ||R||: " << R.magnitude() << endl;
  cout << " inverse: " << inverseDistance << endl;
  cout << " dSdt: " << dSdt << endl;
  
  Real dx = 0.1;
  Real originalPower = _top.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] += dx;
    QUATERNION topEvalNew = iterate * _top.evaluatePowerFactored(iterate);
    QUATERNION bottomEvalNew = iterate * _bottom.evaluatePowerFactored(iterate);
    _top.powersMutable()[whichRoot] = originalPower;
    
    QUATERNION Rnew = (topEvalNew / bottomEvalNew);
    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);

    /*
    Real tanhVersionBefore = -tanh(_tanhAlpha * (Fnew - threshold));
    Real Snew = inverseDistance * tanhVersionBefore;
    */
    // DEBUG: try a simpler S
    Real Snew = -tanh(_tanhAlpha * inverseDistance * (Fnew - _tanhThreshold));

    Real numericalDt = (Snew - S) / dx;

    Real diff = (numericalDt - dSdt) / dSdt;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical dt: " << numericalDt << "\t Snew: " << Snew << "\t S: " << S << endl;
    dx *= 0.1;
  }
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalLogTopDerivative(int whichRoot)
{
  //QUATERNION iterate(1,0.1,0.2,0.3);
  QUATERNION iterate(1.36, -0.56, 1.36, 0);
  //iterate *= 0.01;
  QUATERNION topEval    = iterate * _top.evaluatePowerFactored(iterate);
  QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
  QUATERNION R = (topEval / bottomEval);
  Real logR = log(R.magnitude());

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of TOP root " << whichRoot << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " ========================================================== " << endl;

  cout << " top: " << topEval << endl;
  cout << " bottom: " << bottomEval << endl;
  cout << " iterate: " << iterate << endl;
  cout << " R: " << R << endl;
  cout << " ||R||: " << R.magnitude() << endl;

  // try for the bottom root
  //QUATERNION bottomDerivative = _bottom.powerDerivative(iterate, 0);
  //QUATERNION Rderivative = topEval / (iterate * bottomDerivative);

  // try for the top root
  QUATERNION topDerivative = _top.powerDerivative(iterate, whichRoot);
  QUATERNION Rderivative = (iterate * topDerivative) / bottomEval;

  Real logDerivative = (1.0 / R.dot(R)) * (R.dot(Rderivative));

  cout << " log(||R||) derivative: " << logDerivative << endl;
  
  Real dx = 0.1;
  //Real originalPower = _bottom.powers()[0];
  Real originalPower = _top.powers()[0];
  for (int x = 0; x < 10; x++)
  {
    //_bottom.powersMutable()[0] += dx;
    //QUATERNION topEvalNew = iterate * _top.evaluatePowerFactored(iterate);
    //QUATERNION bottomEvalNew = iterate * _bottom.evaluatePowerFactored(iterate);
    //_bottom.powersMutable()[0] = originalPower;
    
    _top.powersMutable()[whichRoot] += dx;
    QUATERNION topEvalNew = iterate * _top.evaluatePowerFactored(iterate);
    QUATERNION bottomEvalNew = iterate * _bottom.evaluatePowerFactored(iterate);
    _top.powersMutable()[whichRoot] = originalPower;
    
    QUATERNION Rnew = (topEvalNew / bottomEvalNew);
    Real logNew = log(Rnew.magnitude());
    Real logDiff = (logNew - logR) / dx;

    Real diff = (logDiff - logDerivative) / logDerivative;

    //cout << " original: " << R.magnitude() << endl;
    //cout << " perturbed: " << Rnew.magnitude() << endl;
    //cout << " analytical: " << logDerivative << endl;
    //cout << " numerical: " << Rdiff << endl;
    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical: " << logDiff << "\t perturbed: " << logNew << "\t fixed: " << logR << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalTopDerivative(int whichRoot)
{
  //QUATERNION iterate(1,0.1,0.2,0.3);
  QUATERNION iterate(1.36, -0.56, 1.36, 0);
  //iterate *= 0.01;
  QUATERNION topEval    = iterate * _top.evaluatePowerFactored(iterate);
  QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
  QUATERNION R = (topEval / bottomEval);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of TOP root " << whichRoot << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " ========================================================== " << endl;

  cout << " top: " << topEval << endl;
  cout << " bottom: " << bottomEval << endl;
  cout << " iterate: " << iterate << endl;
  cout << " R: " << R << endl;
  cout << " ||R||: " << R.magnitude() << endl;

  // try for the bottom root
  //QUATERNION bottomDerivative = _bottom.powerDerivative(iterate, 0);
  //QUATERNION Rderivative = topEval / (iterate * bottomDerivative);

  // try for the top root
  QUATERNION topDerivative = _top.powerDerivative(iterate, whichRoot);
  QUATERNION Rderivative = (iterate * topDerivative) / bottomEval;

  Real R2derivative = (1.0 / R.magnitude()) * (R.dot(Rderivative));

  cout << " ||R|| derivative: " << R2derivative << endl;
  
  Real dx = 0.1;
  //Real originalPower = _bottom.powers()[0];
  Real originalPower = _top.powers()[0];
  for (int x = 0; x < 10; x++)
  {
    //_bottom.powersMutable()[0] += dx;
    //QUATERNION topEvalNew = iterate * _top.evaluatePowerFactored(iterate);
    //QUATERNION bottomEvalNew = iterate * _bottom.evaluatePowerFactored(iterate);
    //_bottom.powersMutable()[0] = originalPower;
    
    _top.powersMutable()[whichRoot] += dx;
    QUATERNION topEvalNew = iterate * _top.evaluatePowerFactored(iterate);
    QUATERNION bottomEvalNew = iterate * _bottom.evaluatePowerFactored(iterate);
    _top.powersMutable()[whichRoot] = originalPower;
    
    QUATERNION Rnew = (topEvalNew / bottomEvalNew);
    Real Rdiff = (Rnew.magnitude() - R.magnitude()) / dx;

    Real diff = (Rdiff - R2derivative) / R2derivative;

    //cout << " original: " << R.magnitude() << endl;
    //cout << " perturbed: " << Rnew.magnitude() << endl;
    //cout << " analytical: " << R2derivative << endl;
    //cout << " numerical: " << Rdiff << endl;
    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical: " << Rdiff << "\t perturbed: " << Rnew.magnitude() << "\t fixed: " << R.magnitude() << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a rational function
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugRationalBottomDerivative(int whichRoot)
{
  //QUATERNION iterate(1,0.1,0.2,0.3);
  QUATERNION iterate(10,-20,30,-40);
  //iterate *= 0.01;
  QUATERNION topEval    = iterate * _top.evaluatePowerFactored(iterate);
  QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
  QUATERNION R = (topEval / bottomEval);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of BOTTOM root " << whichRoot << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " ========================================================== " << endl;

  cout << " top: " << topEval << endl;
  cout << " bottom: " << bottomEval << endl;
  cout << " iterate: " << iterate << endl;
  cout << " R: " << R << endl;
  cout << " ||R||: " << R.magnitude() << endl;

  // try for the bottom root
  QUATERNION bottomDerivative = _bottom.inversePowerDerivative(iterate, 0);
  QUATERNION Rderivative = (topEval / (iterate * bottomDerivative));

  Real R2derivative = (1.0 / R.magnitude()) * (R.dot(Rderivative));

  cout << " ||R|| derivative: " << R2derivative << endl;
  
  Real dx = 0.1;
  Real originalPower = _bottom.powers()[0];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[0] += dx;
    QUATERNION topEvalNew    = iterate * _top.evaluatePowerFactored(iterate);
    QUATERNION bottomEvalNew = iterate * _bottom.evaluatePowerFactored(iterate);
    _bottom.powersMutable()[0] = originalPower;
    
    QUATERNION Rnew = (topEvalNew / bottomEvalNew);
    Real Rdiff = (Rnew.magnitude() - R.magnitude()) / dx;

    Real diff = (Rdiff - R2derivative) / R2derivative;

    cout << " diff: " << diff << " dx: " << dx << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopConformalHessian(const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  //const int xDebug = 10;
  //const int yDebug = 10;
  //const int zDebug = 10;
  //const int index = xDebug + yDebug * _fractal.xRes() + zDebug * _fractal.slabSize();

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    Real xMin = -_lengths[0] * 0.5;
    Real yMin = -_lengths[1] * 0.5; 
    Real zMin = -_lengths[2] * 0.5;
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / _fractal.xRes() * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / _fractal.yRes() * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / _fractal.zRes() * _lengths[2] + zMin + _center[2];

    //Real weight = _weightField(cellPosition);
    Real weight = _weightField(_fractal.cellCenter(x));

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], 0);
    const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION left(1,0,0,0);
    QUATERNION right(1,0,0,0);
    const Real powerScalar = _top.powerScalar();
  
    for (int y = 0; y < whichRoot; y++)
      left *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);
    for (int y = whichRoot + 1; y < _top.totalRoots(); y++)
      right *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    QUATERNION base = (iterate - _top.roots()[whichRoot]);
    QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
    QUATERNION middle = powerScalar * power * base.log();
    QUATERNION topDerivative = left * middle * right;
    QUATERNION Rderivative = (topDerivative) / bottomEval;
    Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    final += conformalConformal * dFdt;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopConformalHessianFast(const int whichRoot) const
{
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;

  const bool lastRootTop = !(whichRoot < (_top.totalRoots() - 1));
  const Real powerScalar = _top.powerScalar();

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    Real xMin = -_lengths[0] * 0.5;
    Real yMin = -_lengths[1] * 0.5; 
    Real zMin = -_lengths[2] * 0.5;
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / _fractal.xRes() * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / _fractal.yRes() * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / _fractal.zRes() * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));
    const vector<QUATERNION>& forwardCacheTop = _topCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& backwardCacheTop = _topCache.backward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(xIndex,yIndex,zIndex);

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], 0);
    const QUATERNION topEval    = forwardCacheTop.back();
    const QUATERNION bottomEval = forwardCacheBottom.back();
    //const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    //const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION left = forwardCacheTop[whichRoot];
    QUATERNION right = (!lastRootTop) ? backwardCacheTop[whichRoot + 1] : QUATERNION(1,0,0,0);
    /*
    QUATERNION left(1,0,0,0);
    QUATERNION right(1,0,0,0);
    for (int y = 0; y < whichRoot; y++)
      left *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);
    for (int y = whichRoot + 1; y < _top.totalRoots(); y++)
      right *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);
      */

    QUATERNION base = (iterate - _top.roots()[whichRoot]);
#if 1
    QUATERNION middle = powerScalar * base.log();
#else
    QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
    QUATERNION middle = powerScalar * power * base.log();
#endif
    QUATERNION topDerivative = left * middle * right;
    QUATERNION Rderivative = (topDerivative) / bottomEval;
    Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    final += conformalConformal * dFdt;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopDiagonalHessianFast(const int whichRoot) const
{
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;
  
  const bool lastRootTop = !(whichRoot < (_top.totalRoots() - 1));
  const Real powerScalar = _top.powerScalar();

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));
    const vector<QUATERNION>& forwardCacheTop = _topCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& backwardCacheTop = _topCache.backward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(xIndex,yIndex,zIndex);
    
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = forwardCacheTop.back();
    const QUATERNION bottomEval = forwardCacheBottom.back();
    const QUATERNION R = (topEval / bottomEval);

    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION left = forwardCacheTop[whichRoot];
    QUATERNION right = (!lastRootTop) ? backwardCacheTop[whichRoot + 1] : QUATERNION(1,0,0,0);

    //const vector<QUATERNION>& logCache = _topLogCache(xIndex,yIndex,zIndex);
    //QUATERNION logTerm = logCache[whichRoot];
    QUATERNION base = (iterate - _top.roots()[whichRoot]);
    QUATERNION logTerm = base.log();

    QUATERNION middle = powerScalar * logTerm;
    QUATERNION middle2 = powerScalar * powerScalar * logTerm * logTerm;

    QUATERNION Rderivative2 = (left * middle2 * right) / bottomEval;

    QUATERNION topDerivative = left * middle * right;
    QUATERNION Rderivative = (topDerivative) / bottomEval;
    Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    Real dSdcdt= conformalConformal * dFdt;
    Real dSdFdt = dSdcdt; // i.e. the conformal-top derivative

    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (Rderivative.dot(Rderivative) + R.dot(Rderivative2)) * RdotR;
    Real term2 = R.dot(Rderivative) * (2.0 * R.dot(Rderivative));
    Real d2Fdt2 = (term1 - term2) / (RdotRSq);

    final += dSdFdt * dFdt + dSdF * d2Fdt2;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopDiagonalHessian(const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  //const int xDebug = 10;
  //const int yDebug = 10;
  //const int zDebug = 10;
  //const int index = xDebug + yDebug * _fractal.xRes() + zDebug * _fractal.slabSize();
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    //Real weight = _weightField(cellPosition);
    Real weight = _weightField(_fractal.cellCenter(x));

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION left(1,0,0,0);
    QUATERNION right(1,0,0,0);
    const Real powerScalar = _top.powerScalar();
  
    for (int y = 0; y < whichRoot; y++)
      left *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);
    for (int y = whichRoot + 1; y < _top.totalRoots(); y++)
      right *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    QUATERNION base = (iterate - _top.roots()[whichRoot]);
    QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
    QUATERNION logTerm = base.log();
    QUATERNION middle = powerScalar * power * logTerm;

    QUATERNION middle2 = powerScalar * powerScalar * power * logTerm * logTerm;
    QUATERNION Rderivative2 = (left * middle2 * right) / bottomEval;

    QUATERNION topDerivative = left * middle * right;
    QUATERNION Rderivative = (topDerivative) / bottomEval;
    Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    Real dSdcdt= conformalConformal * dFdt;
    Real dSdFdt = dSdcdt; // i.e. the conformal-top derivative

    //Real RdotRderivative = R.dot(Rderivative);
    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (Rderivative.dot(Rderivative) + R.dot(Rderivative2)) * RdotR;
    Real term2 = R.dot(Rderivative) * (2.0 * R.dot(Rderivative));
    Real d2Fdt2 = (term1 - term2) / (RdotRSq);

    final += dSdFdt * dFdt + dSdF * d2Fdt2;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeMixedTopHessian(const int firstRoot, const int secondRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
    const Real powerScalar = _top.powerScalar();

    QUATERNION leftFirst(1,0,0,0);
    QUATERNION rightFirst(1,0,0,0);
    for (int y = 0; y < firstRoot; y++)
      leftFirst *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);
    for (int y = firstRoot + 1; y < _top.totalRoots(); y++)
      rightFirst *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    QUATERNION leftSecond(1,0,0,0);
    QUATERNION rightSecond(1,0,0,0);
    for (int y = 0; y < secondRoot; y++)
      leftSecond *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);
    for (int y = secondRoot + 1; y < _top.totalRoots(); y++)
      rightSecond *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    const QUATERNION baseFirst = (iterate - _top.roots()[firstRoot]);
    const QUATERNION powerFirst = baseFirst.pow(powerScalar * _top.powers()[firstRoot]);
    const QUATERNION middleFirst = powerScalar * powerFirst * baseFirst.log();
    const QUATERNION topDerivativeFirst = leftFirst * middleFirst * rightFirst;
    QUATERNION RderivativeFirst = (topDerivativeFirst) / bottomEval;

    const QUATERNION baseSecond = (iterate - _top.roots()[secondRoot]);
    const QUATERNION powerSecond = baseSecond.pow(powerScalar * _top.powers()[secondRoot]);
    const QUATERNION middleSecond = powerScalar * powerSecond * baseSecond.log();
    const QUATERNION topDerivativeSecond = leftSecond * middleSecond * rightSecond;
    QUATERNION RderivativeSecond = (topDerivativeSecond) / bottomEval;

    QUATERNION T_dti_dtj(1,0,0,0);
    for (int y = 0; y <= firstRoot; y++)
      T_dti_dtj *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    T_dti_dtj *= powerScalar * baseFirst.log();

    for (int y = firstRoot + 1; y <= secondRoot; y++)
      T_dti_dtj *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    T_dti_dtj *= powerScalar * baseSecond.log();

    for (int y = secondRoot + 1; y < _top.totalRoots(); y++)
      T_dti_dtj *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    QUATERNION R_dti_dtj = T_dti_dtj / bottomEval;

    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (RderivativeSecond.dot(RderivativeFirst) + R.dot(R_dti_dtj)) * RdotR;
    Real term2 = R.dot(RderivativeFirst) * (2.0 * R.dot(RderivativeSecond));
    Real d2F_dti_dtj = (term1 - term2) / (RdotRSq);
 
    Real dFdtj = (1.0 / R.dot(R)) * (R.dot(RderivativeSecond));
    Real dFdti = (1.0 / R.dot(R)) * (R.dot(RderivativeFirst));
    Real dSdFdtj = conformalConformal * dFdtj;
    Real analyticHessian = dSdFdtj * dFdti + dSdF * d2F_dti_dtj;

    final += analyticHessian;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeMixedTopHessianFast(const int firstRoot, const int secondRoot) const
{
  Real final = 0;
  //const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;
  //const bool lastRootFirst  = !(firstRoot  < (_top.totalRoots() - 1));

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    //Real weight = _weightField(_fractal.cellCenter(x));
    const vector<QUATERNION>& forwardCacheTop = _topCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& backwardCacheTop = _topCache.backward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(xIndex,yIndex,zIndex);

    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = forwardCacheTop.back();
    const QUATERNION bottomEval = forwardCacheBottom.back();
    const QUATERNION bottomEvalInverse = bottomEval.inverse();
    const QUATERNION R = topEval * bottomEvalInverse;

    /*
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
    */
    Real dSdF = _dSdF[x];
    Real conformalConformal = _conformalConformal[x];
    const Real powerScalar = _top.powerScalar();

    //QUATERNION leftFirst = forwardCacheTop[firstRoot];
    //QUATERNION rightFirst = (!lastRootFirst) ? backwardCacheTop[firstRoot + 1] : QUATERNION(1,0,0,0);

    QUATERNION leftSecond = forwardCacheTop[secondRoot];
    QUATERNION rightSecond = backwardCacheTop[secondRoot + 1];

    const vector<QUATERNION>& logCache = _topLogCache(xIndex,yIndex,zIndex);
    const QUATERNION baseFirstLog = logCache[firstRoot];

    /*
    const QUATERNION baseFirst = (iterate - _top.roots()[firstRoot]);
    const QUATERNION baseFirstLog = baseFirst.log();
    //const QUATERNION baseFirstLog = logCache[firstRoot];
    const QUATERNION middleFirst = powerScalar * baseFirstLog;
    const QUATERNION topDerivativeFirst = leftFirst * middleFirst * rightFirst;
    //QUATERNION RderivativeFirst = topDerivativeFirst * bottomEvalInverse;
    */

    const vector<QUATERNION>& topDerivatives = _topDerivativeCache(xIndex, yIndex, zIndex);
    //QUATERNION RderivativeFirst = topDerivatives[firstRoot] * bottomEvalInverse;
    QUATERNION RderivativeFirst = topDerivatives[firstRoot];

    const QUATERNION baseSecond = (iterate - _top.roots()[secondRoot]);
    const QUATERNION baseSecondLog = baseSecond.log();
    //const QUATERNION baseSecondLog = logCache[secondRoot];
    const QUATERNION middleSecond = powerScalar * baseSecondLog;
    const QUATERNION topDerivativeSecond = leftSecond * middleSecond * rightSecond;
    QUATERNION RderivativeSecond = topDerivativeSecond * bottomEvalInverse;

    QUATERNION T_dti_dtj = forwardCacheTop[firstRoot];
    T_dti_dtj *= powerScalar * baseFirstLog;
    
    for (int y = firstRoot + 1; y <= secondRoot; y++)
      T_dti_dtj *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    T_dti_dtj *= powerScalar * baseSecondLog;
    T_dti_dtj *= backwardCacheTop[secondRoot + 1];

    QUATERNION R_dti_dtj = T_dti_dtj * bottomEvalInverse;

    Real RdotR = _RdotR[x];
    Real RdotRSq = RdotR * RdotR;
    Real RdotRInv = 1.0 / RdotR;

    Real term1 = (RderivativeSecond.dot(RderivativeFirst) + R.dot(R_dti_dtj)) * RdotR;
    Real term2 = R.dot(RderivativeFirst) * (2.0 * R.dot(RderivativeSecond));
    Real d2F_dti_dtj = (term1 - term2) / (RdotRSq);
 
    Real dFdtj = RdotRInv * (R.dot(RderivativeSecond));
    Real dFdti = RdotRInv * (R.dot(RderivativeFirst));
    Real dSdFdtj = conformalConformal * dFdtj;
    Real analyticHessian = dSdFdtj * dFdti + dSdF * d2F_dti_dtj;
    final += analyticHessian;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeMixedBottomHessian(const int firstRoot, const int secondRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
    const Real powerScalar = _bottom.powerScalar();

    QUATERNION leftFirst(1,0,0,0);
    QUATERNION rightFirst(1,0,0,0);
    for (int y = 0; y < firstRoot; y++)
      leftFirst *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = firstRoot + 1; y < _bottom.totalRoots(); y++)
      rightFirst *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    QUATERNION leftSecond(1,0,0,0);
    QUATERNION rightSecond(1,0,0,0);
    for (int y = 0; y < secondRoot; y++)
      leftSecond *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = secondRoot + 1; y < _bottom.totalRoots(); y++)
      rightSecond *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    const QUATERNION baseFirst = (iterate - _bottom.roots()[firstRoot]);
    const QUATERNION powerFirst = baseFirst.pow(powerScalar * _bottom.powers()[firstRoot]);
    const QUATERNION middleFirst = (1.0 / powerScalar) * powerFirst * baseFirst.log().inverse();
    const QUATERNION bottomDerivativeFirst = leftFirst * middleFirst * rightFirst;
    QUATERNION RderivativeFirst = topEval / bottomDerivativeFirst;

    const QUATERNION baseSecond = (iterate - _bottom.roots()[secondRoot]);
    const QUATERNION powerSecond = baseSecond.pow(powerScalar * _bottom.powers()[secondRoot]);
    const QUATERNION middleSecond = (1.0 / powerScalar) * powerSecond * baseSecond.log().inverse();
    const QUATERNION bottomDerivativeSecond = leftSecond * middleSecond * rightSecond;
    QUATERNION RderivativeSecond = topEval / bottomDerivativeSecond;

    QUATERNION B_dbi_dbj(1,0,0,0);
    for (int y = 0; y <= firstRoot; y++)
      B_dbi_dbj *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    B_dbi_dbj *= (1.0 / powerScalar) * baseFirst.log().inverse();

    for (int y = firstRoot + 1; y <= secondRoot; y++)
      B_dbi_dbj *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    B_dbi_dbj *= (1.0 / powerScalar) * baseSecond.log().inverse();

    for (int y = secondRoot + 1; y < _bottom.totalRoots(); y++)
      B_dbi_dbj *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    QUATERNION R_dbi_dbj = topEval / B_dbi_dbj;

    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (RderivativeSecond.dot(RderivativeFirst) + R.dot(R_dbi_dbj)) * RdotR;
    Real term2 = R.dot(RderivativeFirst) * (2.0 * R.dot(RderivativeSecond));
    Real d2F_dbi_dbj = (term1 - term2) / (RdotRSq);
 
    Real dFdtj = (1.0 / R.dot(R)) * (R.dot(RderivativeSecond));
    Real dFdti = (1.0 / R.dot(R)) * (R.dot(RderivativeFirst));
    Real dSdFdtj = conformalConformal * dFdtj;
    Real analyticHessian = dSdFdtj * dFdti + dSdF * d2F_dbi_dbj;

    final += analyticHessian;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeMixedBottomHessianFast(const int firstRoot, const int secondRoot) const
{
  TIMER functionTimer(__FUNCTION__);

  //assert(firstRoot < secondRoot);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  const bool lastRootFirst  = !(firstRoot  < (_bottom.totalRoots() - 1));
  const bool lastRootSecond = !(secondRoot < (_bottom.totalRoots() - 1));

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));

    const vector<QUATERNION>& forwardCacheTop = _topCache.forward(xIndex,yIndex,zIndex);
    //const vector<QUATERNION>& backwardCacheTop = _topCache.backward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& backwardCacheBottom = _bottomCache.backward(xIndex,yIndex,zIndex);

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    //const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    //const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION topEval    = forwardCacheTop.back();
    const QUATERNION bottomEval = forwardCacheBottom.back();
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
    const Real powerScalar = _bottom.powerScalar();

#if 1
    QUATERNION leftFirst = forwardCacheBottom[firstRoot];
    //QUATERNION leftFirst = forwardCacheBottom[firstRoot - 1];
    QUATERNION rightFirst = !lastRootFirst ? backwardCacheBottom[firstRoot + 1] : QUATERNION(1,0,0,0);
#else
    QUATERNION leftFirst(1,0,0,0);
    QUATERNION rightFirst(1,0,0,0);
    for (int y = 0; y < firstRoot; y++)
      leftFirst *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = firstRoot + 1; y < _bottom.totalRoots(); y++)
      rightFirst *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
#endif

#if 1
    QUATERNION leftSecond = forwardCacheBottom[secondRoot];
    //QUATERNION leftSecond = forwardCacheBottom[secondRoot - 1];
    //QUATERNION rightSecond = backwardCacheBottom[secondRoot + 1];
    QUATERNION rightSecond = !lastRootSecond ? backwardCacheBottom[secondRoot + 1] : QUATERNION(1,0,0,0);
#else
    QUATERNION leftSecond(1,0,0,0);
    QUATERNION rightSecond(1,0,0,0);
    for (int y = 0; y < secondRoot; y++)
      leftSecond *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = secondRoot + 1; y < _bottom.totalRoots(); y++)
      rightSecond *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
#endif

    const vector<QUATERNION>& logCache = _bottomLogCache(xIndex,yIndex,zIndex);
    //const QUATERNION baseFirst = (iterate - _bottom.roots()[firstRoot]);
    //const QUATERNION baseFirstLogInverse = baseFirst.log().inverse();
    const QUATERNION baseFirstLogInverse = logCache[firstRoot];

#if 1
    const QUATERNION middleFirst = (1.0 / powerScalar) * baseFirstLogInverse;
#else
    const QUATERNION powerFirst = baseFirst.pow(powerScalar * _bottom.powers()[firstRoot]);
    const QUATERNION middleFirst = (1.0 / powerScalar) * powerFirst * baseFirst.log().inverse();
#endif
    const QUATERNION bottomDerivativeFirst = leftFirst * middleFirst * rightFirst;
    QUATERNION RderivativeFirst = topEval / bottomDerivativeFirst;

    //const QUATERNION baseSecond = (iterate - _bottom.roots()[secondRoot]);
    //const QUATERNION baseSecondLogInverse = baseSecond.log().inverse();
    const QUATERNION baseSecondLogInverse = logCache[secondRoot];

#if 1
    const QUATERNION middleSecond = (1.0 / powerScalar) * baseSecondLogInverse;
#else
    const QUATERNION powerSecond = baseSecond.pow(powerScalar * _bottom.powers()[secondRoot]);
    const QUATERNION middleSecond = (1.0 / powerScalar) * powerSecond * baseSecond.log().inverse();
#endif
    const QUATERNION bottomDerivativeSecond = leftSecond * middleSecond * rightSecond;
    QUATERNION RderivativeSecond = topEval / bottomDerivativeSecond;

#if 1
    QUATERNION B_dbi_dbj = forwardCacheBottom[firstRoot];
    //QUATERNION B_dbi_dbj = (!lastRootFirst) ? forwardCacheBottom[firstRoot] : QUATERNION(1,0,0,0);
#else
    QUATERNION B_dbi_dbj(1,0,0,0);
    for (int y = 0; y <= firstRoot; y++)
      B_dbi_dbj *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
#endif

    B_dbi_dbj *= (1.0 / powerScalar) * baseFirstLogInverse;

    for (int y = firstRoot + 1; y <= secondRoot; y++)
      B_dbi_dbj *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    B_dbi_dbj *= (1.0 / powerScalar) * baseSecondLogInverse;
   
#if 1
    //B_dbi_dbj *= backwardCacheBottom[secondRoot + 1];
    B_dbi_dbj *= (!lastRootSecond) ? backwardCacheBottom[secondRoot + 1] : QUATERNION(1,0,0,0);
#else
    for (int y = secondRoot + 1; y < _bottom.totalRoots(); y++)
      B_dbi_dbj *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
#endif

    QUATERNION R_dbi_dbj = topEval / B_dbi_dbj;

    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (RderivativeSecond.dot(RderivativeFirst) + R.dot(R_dbi_dbj)) * RdotR;
    Real term2 = R.dot(RderivativeFirst) * (2.0 * R.dot(RderivativeSecond));
    Real d2F_dbi_dbj = (term1 - term2) / (RdotRSq);
 
    Real dFdtj = (1.0 / R.dot(R)) * (R.dot(RderivativeSecond));
    Real dFdti = (1.0 / R.dot(R)) * (R.dot(RderivativeFirst));
    Real dSdFdtj = conformalConformal * dFdtj;
    Real analyticHessian = dSdFdtj * dFdti + dSdF * d2F_dbi_dbj;

    final += analyticHessian;

#if 0
    if (isnan(analyticHessian))
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " NaN found! first: " << firstRoot << " second: " << secondRoot << endl;
      exit(0);
    }
#endif
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopBottomHessian(const int topRoot, const int bottomRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
    const Real powerScalarTop = _top.powerScalar();
    const Real powerScalarBottom = _bottom.powerScalar();

    QUATERNION leftBottom(1,0,0,0);
    QUATERNION rightBottom(1,0,0,0);
    for (int y = 0; y < bottomRoot; y++)
      leftBottom *= (iterate - _bottom.roots()[y]).pow(powerScalarBottom * _bottom.powers()[y]);
    for (int y = bottomRoot + 1; y < _bottom.totalRoots(); y++)
      rightBottom *= (iterate - _bottom.roots()[y]).pow(powerScalarBottom * _bottom.powers()[y]);

    QUATERNION leftTop(1,0,0,0);
    QUATERNION rightTop(1,0,0,0);
    for (int y = 0; y < topRoot; y++)
      leftTop *= (iterate - _top.roots()[y]).pow(powerScalarTop * _top.powers()[y]);
    for (int y = topRoot + 1; y < _top.totalRoots(); y++)
      rightTop *= (iterate - _top.roots()[y]).pow(powerScalarTop * _top.powers()[y]);

    const QUATERNION baseBottom = (iterate - _bottom.roots()[bottomRoot]);
    const QUATERNION powerBottom = baseBottom.pow(powerScalarBottom * _bottom.powers()[bottomRoot]);
    const QUATERNION middleBottom = QUATERNION(-1,0,0,0) * (1.0 / powerScalarBottom) * powerBottom * baseBottom.log().inverse();
    const QUATERNION bottomDerivative = leftBottom * middleBottom * rightBottom;
    QUATERNION RderivativeBottom = topEval / bottomDerivative;

    const QUATERNION baseTop = (iterate - _top.roots()[topRoot]);
    const QUATERNION powerTop = baseTop.pow(powerScalarTop * _top.powers()[topRoot]);
    const QUATERNION middleTop = powerScalarTop * powerTop * baseTop.log();
    const QUATERNION topDerivative = leftTop * middleTop * rightTop;
    QUATERNION RderivativeTop = topDerivative / bottomEval;

    QUATERNION R_dbi_dtj = topDerivative / bottomDerivative;

    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (RderivativeTop.dot(RderivativeBottom) + R.dot(R_dbi_dtj)) * RdotR;
    Real term2 = R.dot(RderivativeBottom) * (2.0 * R.dot(RderivativeTop));
    Real d2F_dbi_dtj = (term1 - term2) / (RdotRSq);
 
    Real dFdtj = (1.0 / R.dot(R)) * (R.dot(RderivativeTop));
    Real dFdbi = (1.0 / R.dot(R)) * (R.dot(RderivativeBottom));
    Real dSdFdtj = conformalConformal * dFdtj;
    Real analyticHessian = dSdFdtj * dFdbi + dSdF * d2F_dbi_dtj;

    final += analyticHessian;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopBottomHessianFast(const int topRoot, const int bottomRoot) const
{
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  const Real powerScalarTop = _top.powerScalar();
  const Real powerScalarBottom = _bottom.powerScalar();
  const Real powerScalarBottomInv = 1.0 / powerScalarBottom;

  const bool lastRootBottom = !(bottomRoot < (_bottom.totalRoots() - 1));
  const bool lastRootTop    = !(topRoot < (_top.totalRoots() - 1));

  const vector<vector<QUATERNION> >& forwardsTop = _topCache.forwards();
  const vector<vector<QUATERNION> >& backwardsTop = _topCache.backwards();
  const vector<vector<QUATERNION> >& forwardsBottom = _bottomCache.forwards();
  const vector<vector<QUATERNION> >& backwardsBottom = _bottomCache.backwards();

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));

    // 20% of time spent in retrival here somehow?
    //const vector<QUATERNION>& forwardCacheTop = _topCache.forward(xIndex,yIndex,zIndex);
    //const vector<QUATERNION>& backwardCacheTop = _topCache.backward(xIndex,yIndex,zIndex);
    //const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(xIndex,yIndex,zIndex);
    //const vector<QUATERNION>& backwardCacheBottom = _bottomCache.backward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& forwardCacheTop = forwardsTop[x];
    const vector<QUATERNION>& backwardCacheTop = backwardsTop[x];
    const vector<QUATERNION>& forwardCacheBottom = forwardsBottom[x];
    const vector<QUATERNION>& backwardCacheBottom = backwardsBottom[x];

    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = forwardCacheTop.back();
    const QUATERNION bottomEval = forwardCacheBottom.back();
    const QUATERNION bottomEvalInverse = bottomEval.inverse();
    const QUATERNION R = topEval * bottomEvalInverse;
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION leftBottom = forwardCacheBottom[bottomRoot];
    QUATERNION rightBottom = !lastRootBottom ? backwardCacheBottom[bottomRoot + 1] : QUATERNION(1,0,0,0);

    QUATERNION leftTop = forwardCacheTop[topRoot];
    QUATERNION rightTop = !lastRootTop   ? backwardCacheTop[topRoot + 1] : QUATERNION(1,0,0,0);

    const QUATERNION baseBottom = (iterate - _bottom.roots()[bottomRoot]);
    // 10% of running time spent in the log
    const QUATERNION middleBottom = QUATERNION(-1,0,0,0) * powerScalarBottomInv * baseBottom.log().inverse();
    const QUATERNION bottomDerivative = leftBottom * middleBottom * rightBottom;
    const QUATERNION bottomDerivativeInverse = bottomDerivative.inverse();
    QUATERNION RderivativeBottom = topEval * bottomDerivativeInverse;

    const QUATERNION baseTop = (iterate - _top.roots()[topRoot]);
    // 10% of running time spent in the log
    const QUATERNION middleTop = powerScalarTop * baseTop.log();
    const QUATERNION topDerivative = leftTop * middleTop * rightTop;
    QUATERNION RderivativeTop = topDerivative * bottomEvalInverse;
    QUATERNION R_dbi_dtj = topDerivative * bottomDerivativeInverse;

    // 12% of running time below here
    Real RdotR = R.dot(R);
    Real RdotRInv = 1.0 / RdotR;
    Real RdotRSqInv = RdotRInv * RdotRInv;

    Real term1 = (RderivativeTop.dot(RderivativeBottom) + R.dot(R_dbi_dtj)) * RdotR;
    Real term2 = R.dot(RderivativeBottom) * (2.0 * R.dot(RderivativeTop));
    Real d2F_dbi_dtj = (term1 - term2) * RdotRSqInv;
 
    Real dFdtj = RdotRInv * (R.dot(RderivativeTop));
    Real dFdbi = RdotRInv * (R.dot(RderivativeBottom));
    Real dSdFdtj = conformalConformal * dFdtj;
    Real analyticHessian = dSdFdtj * dFdbi + dSdF * d2F_dbi_dtj;

    final += analyticHessian;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomDiagonalHessian(const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    //Real weight = _weightField(cellPosition);
    Real weight = _weightField(_fractal.cellCenter(x));

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION left(1,0,0,0);
    QUATERNION right(1,0,0,0);
    const Real powerScalar = _bottom.powerScalar();
  
    for (int y = 0; y < whichRoot; y++)
      left *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = whichRoot + 1; y < _bottom.totalRoots(); y++)
      right *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
    QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
    QUATERNION logTermInverse = base.log().inverse();
    QUATERNION middle = powerScalar * power * logTermInverse;

    QUATERNION middle2 = 1.0 / (powerScalar * powerScalar) * power * logTermInverse * logTermInverse;
    QUATERNION Rderivative2 = topEval / (left * middle2 * right);

    QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;
    QUATERNION Rderivative = topEval / bottomDerivative;
    Real dFdb = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    Real dSdcdb= conformalConformal * dFdb;
    Real dSdFdb = dSdcdb; // i.e. the conformal-top derivative

    //Real RdotRderivative = R.dot(Rderivative);
    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (Rderivative.dot(Rderivative) + R.dot(Rderivative2)) * RdotR;
    Real term2 = R.dot(Rderivative) * (2.0 * R.dot(Rderivative));
    Real d2Fdb2 = (term1 - term2) / (RdotRSq);

    final += dSdFdb * dFdb + dSdF * d2Fdb2;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomDiagonalHessianFast(const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  const bool lastRootBottom = !(whichRoot < (_bottom.totalRoots() - 1));
  const Real powerScalar = _bottom.powerScalar();

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / xRes * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / yRes * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / zRes * _lengths[2] + zMin + _center[2];

    //Real weight = _weightField(cellPosition);
    Real weight = _weightField(_fractal.cellCenter(x));
    const vector<QUATERNION>& forwardCacheTop = _topCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& backwardCacheBottom = _bottomCache.backward(xIndex,yIndex,zIndex);

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], _quaternionSlice);
    //const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    //const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION topEval    = forwardCacheTop.back();
    const QUATERNION bottomEval = forwardCacheBottom.back();
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real dSdF = -_tanhAlpha * weight * sechTermSq;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    //QUATERNION left = (whichRoot != 0) ? forwardCacheBottom[whichRoot - 1] : QUATERNION(1,0,0,0);
    QUATERNION left = forwardCacheBottom[whichRoot];
    QUATERNION right = (!lastRootBottom) ? backwardCacheBottom[whichRoot + 1] : QUATERNION(1,0,0,0);
    /*
    QUATERNION left(1,0,0,0);
    QUATERNION right(1,0,0,0);
  
    for (int y = 0; y < whichRoot; y++)
      left *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = whichRoot + 1; y < _bottom.totalRoots(); y++)
      right *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
      */

    QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
    QUATERNION logTermInverse = base.log().inverse();
   
#if 1
    QUATERNION middle = powerScalar * logTermInverse;
    QUATERNION middle2 = 1.0 / (powerScalar * powerScalar) * logTermInverse * logTermInverse;
#else 
    QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
    QUATERNION middle = powerScalar * power * logTermInverse;
    QUATERNION middle2 = 1.0 / (powerScalar * powerScalar) * power * logTermInverse * logTermInverse;
#endif

    QUATERNION Rderivative2 = topEval / (left * middle2 * right);

    QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;
    QUATERNION Rderivative = topEval / bottomDerivative;
    Real dFdb = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    Real dSdcdb= conformalConformal * dFdb;
    Real dSdFdb = dSdcdb; // i.e. the conformal-top derivative

    Real RdotR = R.dot(R);
    Real RdotRSq = RdotR * RdotR;

    Real term1 = (Rderivative.dot(Rderivative) + R.dot(Rderivative2)) * RdotR;
    Real term2 = R.dot(Rderivative) * (2.0 * R.dot(Rderivative));
    Real d2Fdb2 = (term1 - term2) / (RdotRSq);

    final += dSdFdb * dFdb + dSdF * d2Fdb2;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomConformalHessian(const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    Real xMin = -_lengths[0] * 0.5;
    Real yMin = -_lengths[1] * 0.5; 
    Real zMin = -_lengths[2] * 0.5;
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / _fractal.xRes() * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / _fractal.yRes() * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / _fractal.zRes() * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], 0);
    const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION left(1,0,0,0);
    QUATERNION right(1,0,0,0);
    const Real powerScalar = _bottom.powerScalar();
  
    for (int y = 0; y < whichRoot; y++)
      left *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = whichRoot + 1; y < _bottom.totalRoots(); y++)
      right *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
    QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
    QUATERNION middle = (1.0 / powerScalar) * power * base.log().inverse();
    QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;
    QUATERNION Rderivative = topEval / bottomDerivative;
    Real dFdb = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    final += conformalConformal * dFdb;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomConformalHessianFast(const int whichRoot) const
{
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  
  const bool lastRootBottom = !(whichRoot < (_bottom.totalRoots() - 1));
  const Real powerScalar = _bottom.powerScalar();

  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    int zIndex = x / _fractal.slabSize();
    int yIndex = (x - zIndex * _fractal.slabSize()) / _fractal.xRes();
    int xIndex = x - zIndex * _fractal.slabSize() - yIndex * _fractal.xRes();
    Real xMin = -_lengths[0] * 0.5;
    Real yMin = -_lengths[1] * 0.5; 
    Real zMin = -_lengths[2] * 0.5;
    VEC3F cellPosition;
    cellPosition[0] = (Real)xIndex / _fractal.xRes() * _lengths[0] + xMin + _center[0];
    cellPosition[1] = (Real)yIndex / _fractal.yRes() * _lengths[1] + yMin + _center[1];
    cellPosition[2] = (Real)zIndex / _fractal.zRes() * _lengths[2] + zMin + _center[2];

    Real weight = _weightField(_fractal.cellCenter(x));
    const vector<QUATERNION>& forwardCacheTop = _topCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(xIndex,yIndex,zIndex);
    const vector<QUATERNION>& backwardCacheBottom = _bottomCache.backward(xIndex,yIndex,zIndex);

    vector<QUATERNION> forwardDummy, backwardDummy;
    const QUATERNION iterate(cellPosition[0], cellPosition[1], cellPosition[2], 0);
    //const QUATERNION topEval    = _top.evaluateScaledPowerFactored(iterate);
    //const QUATERNION bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION topEval    = forwardCacheTop.back();
    const QUATERNION bottomEval = forwardCacheBottom.back();
    const QUATERNION R = (topEval / bottomEval);
    Real F = log(R.magnitude()) + log(_expScaling);
    Real interior = _tanhAlpha * (F - _tanhThreshold);
    Real tanhTerm = tanh(interior);

    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

    QUATERNION left = forwardCacheBottom[whichRoot];
    QUATERNION right = (!lastRootBottom) ? backwardCacheBottom[whichRoot + 1] : QUATERNION(1,0,0,0);
    /*
    QUATERNION left(1,0,0,0);
    QUATERNION right(1,0,0,0);
    for (int y = 0; y < whichRoot; y++)
      left *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    for (int y = whichRoot + 1; y < _bottom.totalRoots(); y++)
      right *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
      */

    QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
#if 1
    QUATERNION middle = (1.0 / powerScalar) * base.log().inverse();
#else
    QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
    QUATERNION middle = (1.0 / powerScalar) * power * base.log().inverse();
#endif
    QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;
    QUATERNION Rderivative = topEval / bottomDerivative;
    Real dFdb = (1.0 / R.dot(R)) * (R.dot(Rderivative));

    final += conformalConformal * dFdb;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute entries in the Hessian
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeConformalConformalHessian() const
{
  TIMER functionTimer(__FUNCTION__);
  Real final = 0;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);

    Real weight = _weightField(cellPosition);
    Real tanhTerm = tanh(_tanhAlpha * _fractal[x]);
    Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
    final += 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial c \partial c}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianConformalConformal()
{
  int xDebug = 80;
  int yDebug = 54;
  int zDebug = 48;
  VEC3F point = _weightField.cellCenter(xDebug, yDebug, zDebug);

  //QUATERNION iterate(1,0.1,0.2,0.3);
  QUATERNION iterate(point[0], point[1], point[2], 0);
  //iterate *= 0.01;
  QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  QUATERNION R = (topEval / bottomEval);
  
  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  Real sechTermSq = 1 - tanhTerm * tanhTerm;
  Real weight = _weightField(xDebug, yDebug, zDebug);

  Real originalScore = weight * (-tanhTerm);
  Real derivative = -_tanhAlpha * weight * sechTermSq;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of conformal-conformal Hessian" << endl;
  cout << " ========================================================== " << endl;
  cout << " original score: " << originalScore << endl;
  cout << " derivative: " << derivative << endl;
  
  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    Real Fnew = log(R.magnitude()) + log(_expScaling) + dx;
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScoreNew = weight * (-tanhTermNew);

    Real finiteDiff = (newScoreNew - originalScore) / dx;
    Real relative = (finiteDiff - derivative) / derivative;
    cout << "analytic: " << derivative << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  // get the analytic Hessian
  Real hessian = 2.0 * _tanhAlpha * _tanhAlpha * weight * sechTermSq * tanhTerm;

  cout << " Now the Hessian: " << endl;
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    Real Fnew = log(R.magnitude()) + log(_expScaling) + dx;
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real sechTermSqNew = 1 - tanhTermNew * tanhTermNew;
    Real newDerivative = -_tanhAlpha * weight * sechTermSqNew;

    Real finiteDiff = (newDerivative - derivative) / dx;
    Real relative = (hessian - finiteDiff) / hessian;
    cout << "analytic: " << hessian << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  // get the Hessian over the entire field
  cout << " Now the Hessian over the entire field: " << endl;
  Real fullHessian = computeConformalConformalHessian();
  Real originalGradient = computeLogPowerScaledRationalScoreExpGradient();
  dx = 0.1;
  const Real backup = _expScaling;
  for (int x = 0; x < 10; x++)
  {
    _expScaling = backup + dx;
    Real newGradient = computeLogPowerScaledRationalScoreExpGradient();

    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (finiteDiff - fullHessian) / fullHessian;
    cout << "analytic: " << fullHessian << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial c \partial t_i}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianBottomConformal(const int whichRoot)
{
  //int xDebug = 80;
  //int yDebug = 54;
  //int zDebug = 48;
  int xDebug = 0;
  int yDebug = 0;
  int zDebug = 0;
  //const VEC3F point = _weightField.cellCenter(xDebug, yDebug, zDebug);
  const VEC3F point = _fractal.cellCenter(xDebug, yDebug, zDebug);
  const QUATERNION iterate(point[0], point[1], point[2], 0);
  const QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  const QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  const QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  Real weight = _weightField(point);
  Real originalScore = weight * (-tanhTerm);
  cout << " original: " << originalScore << endl;
  cout << " Weight: " << weight << endl;

  Real analyticGradient;
  QUATERNION left(1,0,0,0);
  QUATERNION right(1,0,0,0);
  const Real powerScalar = _bottom.powerScalar();

  for (int x = 0; x < whichRoot; x++)
    left *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);
  for (int x = whichRoot + 1; x < _bottom.totalRoots(); x++)
    right *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);

  QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
  QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
  QUATERNION middle = (1.0 / powerScalar) * power * base.log().inverse();
  QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;

  QUATERNION Rderivative = topEval / bottomDerivative;
  Real dFdb = (1.0 / R.dot(R)) * (R.dot(Rderivative));
  Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
  analyticGradient = dSdF * dFdb;
  Real dSdb = dSdF * dFdb;
  cout << " analytic gradient: " << analyticGradient << endl;

  Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
  Real analyticHessian = conformalConformal * dFdb;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of bottom-conformal Hessian, root " << whichRoot << endl;
  cout << " ========================================================== " << endl;

  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  Real originalPower = _bottom.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] = originalPower + dx;
    QUATERNION bottomEvalNew    = _bottom.evaluatePowerFactored(iterate);
    QUATERNION Rnew = topEval/ bottomEvalNew;

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScore = weight * (-tanhTermNew);

    Real finiteDiff = (newScore - originalScore) / dx;
    Real relative = (finiteDiff - analyticGradient) / analyticGradient;
    cout << "analytic: " << analyticGradient << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  cout << " now the Hessian: " << endl;
  _bottom.powersMutable()[whichRoot] = originalPower;
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    Real Fnew = log(R.magnitude()) + log(_expScaling) + dx;
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real dSdFnew = -_tanhAlpha * weight * (1 - tanhTermNew * tanhTermNew);
    Real dSdbNew = dSdFnew * dFdb;

    Real finiteDiff = (dSdbNew - dSdb) / dx;

    Real relative = (finiteDiff - analyticHessian) / analyticHessian;

    cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  cout << " Over the full grid:" << endl;
  dx = 0.1;
  Real fullHessian = computeBottomConformalHessian(whichRoot);
  
  Real originalGradient = computeLogPowerScaledRationalScoreExpGradient();

  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] = originalPower + dx;
    Real newGradient = computeLogPowerScaledRationalScoreExpGradient();
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    //cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[whichRoot] = originalPower;

  cout << " Checking for Hessian symmetry:" << endl;
  dx = 0.1;
  originalGradient = computeBottomScaledPowerGradient(whichRoot);
  const Real backup = _expScaling;
  for (int x = 0; x < 10; x++)
  {
    _expScaling = backup + dx;

    Real newGradient = computeBottomScaledPowerGradient(whichRoot);
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[whichRoot] = originalPower;
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial t_i \partial t_i}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianTopDiagonal(const int whichRoot)
{
  int xDebug = 0;
  int yDebug = 0;
  int zDebug = 0;
 
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  VEC3F point = _fractal.cellCenter(xDebug, yDebug, zDebug);
  point[0] = (Real)xDebug / xRes * _lengths[0] + xMin + _center[0];
  point[1] = (Real)yDebug / yRes * _lengths[1] + yMin + _center[1];
  point[2] = (Real)zDebug / zRes * _lengths[2] + zMin + _center[2];

  const QUATERNION iterate(point[0], point[1], point[2], 0);
  const QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  const QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  const QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  //Real weight = _weightField(point);
  Real weight = _weightField(_fractal.cellCenter(xDebug, yDebug, zDebug));
  Real originalScore = weight * (-tanhTerm);
  cout << " original: " << originalScore << endl;
  cout << " Weight: " << weight << endl;

  Real analyticGradient;
  QUATERNION left(1,0,0,0);
  QUATERNION right(1,0,0,0);
  const Real powerScalar = _top.powerScalar();

  for (int x = 0; x < whichRoot; x++)
    left *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);
  for (int x = whichRoot + 1; x < _top.totalRoots(); x++)
    right *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);

  const QUATERNION base = (iterate - _top.roots()[whichRoot]);
  const QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
  const QUATERNION middle = powerScalar * power * base.log();
  const QUATERNION topDerivative = left * middle * right;
  QUATERNION Rderivative = (topDerivative) / bottomEval;
  Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
  Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
  analyticGradient = dSdF * dFdt;
  Real dSdt = dSdF * dFdt;
  cout << " analytic gradient: " << analyticGradient << endl;

  Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
  Real dSdcdt= conformalConformal * dFdt;
  Real dSdFdt = dSdcdt; // i.e. the conformal-top derivative

  //Real RdotRderivative = R.dot(Rderivative);
  Real RdotR = R.dot(R);
  Real RdotRSq = RdotR * RdotR;

  QUATERNION logTerm = base.log();
  QUATERNION middle2 = powerScalar * powerScalar * power * logTerm * logTerm;
  QUATERNION Rderivative2 = (left * middle2 * right) / bottomEval;

  Real term1 = (Rderivative.dot(Rderivative) + R.dot(Rderivative2)) * RdotR;
  Real term2 = R.dot(Rderivative) * (2.0 * R.dot(Rderivative));

  Real d2Fdt2 = (term1 - term2) / (RdotRSq);
  Real analyticHessian = dSdFdt * dFdt + dSdF * d2Fdt2;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of top-top Hessian, root " << whichRoot << endl;
  cout << " ========================================================== " << endl;

  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  Real originalPower = _top.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] = originalPower + dx;
    QUATERNION topEvalNew    = _top.evaluateScaledPowerFactored(iterate);
    QUATERNION Rnew = (topEvalNew / bottomEval);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScore = weight * (-tanhTermNew);

    Real finiteDiff = (newScore - originalScore) / dx;
    Real relative = (finiteDiff - analyticGradient) / analyticGradient;
    cout << "analytic: " << analyticGradient << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  cout << " now the Hessian: " << endl;
  _top.powersMutable()[whichRoot] = originalPower;
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] = originalPower + dx;
    const QUATERNION topEvalNew    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEvalNew = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION Rnew = (topEvalNew / bottomEvalNew);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    //Real originalScoreNew = weight * (-tanhTermNew);

    //Real analyticGradientNew;

    QUATERNION baseNew = (iterate - _top.roots()[whichRoot]);
    QUATERNION powerNew = baseNew.pow(powerScalar * _top.powers()[whichRoot]);
    QUATERNION middleNew = powerScalar * powerNew * baseNew.log();
    QUATERNION topDerivativeNew = left * middleNew * right;
    QUATERNION RderivativeNew = (topDerivativeNew) / bottomEval;
    Real dFdtNew = (1.0 / Rnew.dot(Rnew)) * (Rnew.dot(RderivativeNew));
    Real dSdFNew = -_tanhAlpha * weight * (1 - tanhTermNew * tanhTermNew);
    //analyticGradientNew = dSdFNew * dFdtNew;
    Real dSdtNew = dSdFNew * dFdtNew;
    Real finiteDiff = (dSdtNew - dSdt) / dx;
    Real relative = (finiteDiff - analyticHessian) / analyticHessian;
    cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[whichRoot] = originalPower;

  cout << " Over the full grid:" << endl;
  dx = 0.1;
  Real fullHessian = computeTopDiagonalHessian(whichRoot);
  Real originalGradient = computeTopScaledPowerGradient(whichRoot);

  cout << " original: " << originalGradient << " original single sample: " << analyticGradient << " diff: " << originalGradient - analyticGradient << endl;
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] = originalPower + dx;
    refreshPolynomialCache(_topCache, _top, _fractal);    

    Real newGradient = computeTopScaledPowerGradient(whichRoot);
    //cout << " dx: " << dx << endl;
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    //cout << " new gradient: " << newGradient << endl;

    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[whichRoot] = originalPower;
  refreshPolynomialCache(_topCache, _top, _fractal);    
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial t_i \partial t_j}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianMixedTop(const int firstRoot, const int secondRoot)
{
  assert(firstRoot < secondRoot);
  assert(firstRoot != secondRoot);

  int xDebug = 0;
  int yDebug = 0;
  int zDebug = 0;
 
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  VEC3F point = _fractal.cellCenter(xDebug, yDebug, zDebug);
  point[0] = (Real)xDebug / xRes * _lengths[0] + xMin + _center[0];
  point[1] = (Real)yDebug / yRes * _lengths[1] + yMin + _center[1];
  point[2] = (Real)zDebug / zRes * _lengths[2] + zMin + _center[2];

  const QUATERNION iterate(point[0], point[1], point[2], 0);
  const QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  const QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  const QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  //Real weight = _weightField(point);
  Real weight = _weightField(_fractal.cellCenter(xDebug, yDebug, zDebug));
  Real originalScore = weight * (-tanhTerm);
  cout << " original: " << originalScore << endl;
  cout << " Weight: " << weight << endl;

  QUATERNION lefti(1,0,0,0);
  QUATERNION righti(1,0,0,0);
  QUATERNION leftj(1,0,0,0);
  QUATERNION rightj(1,0,0,0);
  const Real powerScalar = _top.powerScalar();

  for (int x = 0; x < firstRoot; x++)
    lefti *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);
  for (int x = firstRoot + 1; x < _top.totalRoots(); x++)
    righti *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);
  
  for (int x = 0; x < secondRoot; x++)
    leftj *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);
  for (int x = secondRoot + 1; x < _top.totalRoots(); x++)
    rightj *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);

  const QUATERNION basei = (iterate - _top.roots()[firstRoot]);
  const QUATERNION poweri = basei.pow(powerScalar * _top.powers()[firstRoot]);
  const QUATERNION middlei = powerScalar * poweri * basei.log();
  const QUATERNION topDerivativei = lefti * middlei * righti;
  QUATERNION Rderivativei = (topDerivativei) / bottomEval;

  const QUATERNION basej = (iterate - _top.roots()[secondRoot]);
  const QUATERNION powerj = basej.pow(powerScalar * _top.powers()[secondRoot]);
  const QUATERNION middlej = powerScalar * powerj * basej.log();
  const QUATERNION topDerivativej = leftj * middlej * rightj;
  QUATERNION Rderivativej = (topDerivativej) / bottomEval;

  Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);

  Real dFdti = (1.0 / R.dot(R)) * (R.dot(Rderivativei));
  Real dSdti = dSdF * dFdti;
  cout << " dS / dt_i: " << dSdti << endl;
  
  Real dFdtj = (1.0 / R.dot(R)) * (R.dot(Rderivativej));
  Real dSdtj = dSdF * dFdtj;
  cout << " dS / dt_j: " << dSdtj << endl;

  Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
  //Real dSdcdti = conformalConformal * dFdti;
  //Real dSdFdti = dSdcdti; // i.e. the conformal-top_i derivative

  Real dSdcdtj = conformalConformal * dFdtj;
  Real dSdFdtj = dSdcdtj; // i.e. the conformal-top_j derivative

  Real RdotR = R.dot(R);
  Real RdotRSq = RdotR * RdotR;

  QUATERNION T_dti_dtj(1,0,0,0);
  for (int x = 0; x <= firstRoot; x++)
    T_dti_dtj *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);

  T_dti_dtj *= powerScalar * basei.log();

  for (int x = firstRoot + 1; x <= secondRoot; x++)
    T_dti_dtj *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);

  T_dti_dtj *= powerScalar * basej.log();

  for (int x = secondRoot + 1; x < _top.totalRoots(); x++)
    T_dti_dtj *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);

  QUATERNION R_dti_dtj = T_dti_dtj / bottomEval;

  Real term1 = (Rderivativej.dot(Rderivativei) + R.dot(R_dti_dtj)) * RdotR;
  Real term2 = R.dot(Rderivativei) * (2.0 * R.dot(Rderivativej));

  Real d2F_dti_dtj = (term1 - term2) / (RdotRSq);
  Real analyticHessian = dSdFdtj * dFdti + dSdF * d2F_dti_dtj;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of mixed top-top Hessian, roots " << firstRoot << ", " << secondRoot << endl;
  cout << " ========================================================== " << endl;

  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  Real originalPoweri = _top.powers()[firstRoot];
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[firstRoot] = originalPoweri + dx;
    QUATERNION topEvalNew    = _top.evaluateScaledPowerFactored(iterate);
    QUATERNION Rnew = (topEvalNew / bottomEval);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScore = weight * (-tanhTermNew);

    Real finiteDiff = (newScore - originalScore) / dx;
    Real relative = (finiteDiff - dSdti) / dSdti;
    cout << "analytic: " << dSdti << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[firstRoot] = originalPoweri;

  cout << " now the Hessian: " << endl;
  Real originalPowerj = _top.powers()[secondRoot];
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[secondRoot] = originalPowerj + dx;
    const QUATERNION topEvalNew    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEvalNew = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION Rnew = (topEvalNew / bottomEvalNew);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);

    //Real analyticGradientNew;

    QUATERNION baseNew = (iterate - _top.roots()[firstRoot]);
    //QUATERNION powerNew = baseNew.pow(powerScalar * _top.powers()[firstRoot]);
    //QUATERNION middleNew = powerScalar * powerNew * baseNew.log();
    
    lefti = QUATERNION(1,0,0,0);
    for (int y = 0; y < firstRoot; y++)
      lefti *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);
    righti = QUATERNION(1,0,0,0);
    for (int y = firstRoot + 1; y < _top.totalRoots(); y++)
      righti *= (iterate - _top.roots()[y]).pow(powerScalar * _top.powers()[y]);

    QUATERNION topDerivativeNew = lefti * middlei * righti;
    QUATERNION RderivativeNew = (topDerivativeNew) / bottomEvalNew;
    Real dFdtNew = (1.0 / Rnew.dot(Rnew)) * (Rnew.dot(RderivativeNew));
    Real dSdFNew = -_tanhAlpha * weight * (1 - tanhTermNew * tanhTermNew);
    //analyticGradientNew = dSdFNew * dFdtNew;
    Real dSdtNew = dSdFNew * dFdtNew;

    Real finiteDiff = (dSdtNew - dSdti) / dx;
    Real relative = (finiteDiff - analyticHessian) / analyticHessian;
    cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[secondRoot] = originalPowerj;

  cout << " Over the full grid:" << endl;
  dx = 0.1;
  Real fullHessian = computeMixedTopHessian(firstRoot, secondRoot);

  Real originalGradient = computeTopScaledPowerGradient(firstRoot);

  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[secondRoot] = originalPowerj + dx;
    refreshPolynomialCache(_topCache, _top, _fractal);

    Real newGradient = computeTopScaledPowerGradient(firstRoot);
    //cout << " dx: " << dx << endl;
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    //cout << " new gradient: " << newGradient << endl;

    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[secondRoot] = originalPowerj;
  refreshPolynomialCache(_topCache, _top, _fractal);
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial b_i \partial b_j}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianMixedBottom(const int firstRoot, const int secondRoot)
{
  assert(firstRoot < secondRoot);
  assert(firstRoot != secondRoot);

  int xDebug = 0;
  int yDebug = 0;
  int zDebug = 0;
 
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  VEC3F point = _fractal.cellCenter(xDebug, yDebug, zDebug);
  point[0] = (Real)xDebug / xRes * _lengths[0] + xMin + _center[0];
  point[1] = (Real)yDebug / yRes * _lengths[1] + yMin + _center[1];
  point[2] = (Real)zDebug / zRes * _lengths[2] + zMin + _center[2];

  const QUATERNION iterate(point[0], point[1], point[2], 0);
  const QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  const QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  const QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  Real weight = _weightField(_fractal.cellCenter(xDebug, yDebug, zDebug));
  Real originalScore = weight * (-tanhTerm);

  QUATERNION lefti(1,0,0,0);
  QUATERNION righti(1,0,0,0);
  QUATERNION leftj(1,0,0,0);
  QUATERNION rightj(1,0,0,0);
  const Real powerScalar = _bottom.powerScalar();

  for (int x = 0; x < firstRoot; x++)
    lefti *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);
  for (int x = firstRoot + 1; x < _bottom.totalRoots(); x++)
    righti *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);
  
  for (int x = 0; x < secondRoot; x++)
    leftj *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);
  for (int x = secondRoot + 1; x < _bottom.totalRoots(); x++)
    rightj *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);

  const QUATERNION basei = (iterate - _bottom.roots()[firstRoot]);
  const QUATERNION poweri = basei.pow(powerScalar * _bottom.powers()[firstRoot]);
  const QUATERNION middlei = QUATERNION(-1,0,0,0) * (1.0 / powerScalar) * poweri * basei.log().inverse();
  const QUATERNION bottomDerivativei = lefti * middlei * righti;
  QUATERNION Rderivativei = topEval / bottomDerivativei;

  const QUATERNION basej = (iterate - _bottom.roots()[secondRoot]);
  const QUATERNION powerj = basej.pow(powerScalar * _bottom.powers()[secondRoot]);
  const QUATERNION middlej = QUATERNION(-1,0,0,0) * (1.0 / powerScalar) * powerj * basej.log().inverse();
  const QUATERNION bottomDerivativej = leftj * middlej * rightj;
  QUATERNION Rderivativej = topEval / bottomDerivativej;

  Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);

  Real dFdbi = (1.0 / R.dot(R)) * (R.dot(Rderivativei));
  Real dSdbi = dSdF * dFdbi;
  cout << " dS / db_i: " << dSdbi << endl;
  
  Real dFdbj = (1.0 / R.dot(R)) * (R.dot(Rderivativej));
  Real dSdbj = dSdF * dFdbj;
  cout << " dS / db_j: " << dSdbj << endl;

  Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
  //Real dSdcdbi = conformalConformal * dFdbi;
  //Real dSdFdbi = dSdcdbi; // i.e. the conformal-top_i derivative

  Real dSdcdbj = conformalConformal * dFdbj;
  Real dSdFdbj = dSdcdbj; // i.e. the conformal-top_j derivative

  Real RdotR = R.dot(R);
  Real RdotRSq = RdotR * RdotR;

  QUATERNION B_dbi_dbj(1,0,0,0);
  for (int x = 0; x <= firstRoot; x++)
    B_dbi_dbj *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);

  B_dbi_dbj *= (1.0 / powerScalar) * basei.log().inverse();

  for (int x = firstRoot + 1; x <= secondRoot; x++)
    B_dbi_dbj *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);

  B_dbi_dbj *= (1.0 / powerScalar) * basej.log().inverse();

  for (int x = secondRoot + 1; x < _bottom.totalRoots(); x++)
    B_dbi_dbj *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);

  QUATERNION R_dbi_dbj = topEval / B_dbi_dbj;

  Real term1 = (Rderivativej.dot(Rderivativei) + R.dot(R_dbi_dbj)) * RdotR;
  Real term2 = R.dot(Rderivativei) * (2.0 * R.dot(Rderivativej));

  Real d2F_dbi_dbj = (term1 - term2) / (RdotRSq);
  Real analyticHessian = dSdFdbj * dFdbi + dSdF * d2F_dbi_dbj;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of mixed bottom-bottom Hessian, roots " << firstRoot << ", " << secondRoot << endl;
  cout << " ========================================================== " << endl;

  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  Real originalPoweri = _bottom.powers()[firstRoot];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[firstRoot] = originalPoweri + dx;
    QUATERNION bottomEvalNew    = _bottom.evaluateScaledPowerFactored(iterate);
    QUATERNION Rnew = (topEval / bottomEvalNew);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScore = weight * (-tanhTermNew);

    Real finiteDiff = (newScore - originalScore) / dx;
    Real relative = (finiteDiff - dSdbi) / dSdbi;
    cout << "analytic: " << dSdbi << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[firstRoot] = originalPoweri;

  cout << " now the Hessian: " << endl;
  Real originalPowerj = _bottom.powers()[secondRoot];
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[secondRoot] = originalPowerj + dx;
    const QUATERNION topEvalNew    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEvalNew = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION Rnew = (topEvalNew / bottomEvalNew);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);

    //Real analyticGradientNew;

    lefti = QUATERNION(1,0,0,0);
    for (int y = 0; y < firstRoot; y++)
      lefti *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);
    righti = QUATERNION(1,0,0,0);
    for (int y = firstRoot + 1; y < _bottom.totalRoots(); y++)
      righti *= (iterate - _bottom.roots()[y]).pow(powerScalar * _bottom.powers()[y]);

    //QUATERNION bottomDerivativeNew = lefti * middleNew * righti;
    QUATERNION bottomDerivativeNew = lefti * middlei * righti;
    QUATERNION RderivativeNew = topEvalNew / bottomDerivativeNew;
    Real dFdbNew = (1.0 / Rnew.dot(Rnew)) * (Rnew.dot(RderivativeNew));
    Real dSdFNew = -_tanhAlpha * weight * (1.0 - tanhTermNew * tanhTermNew);
    //analyticGradientNew = dSdFNew * dFdbNew;
    Real dSdbNew = dSdFNew * dFdbNew;

    Real finiteDiff = (dSdbNew - dSdbi) / dx;
    Real relative = (finiteDiff - analyticHessian) / analyticHessian;
    cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[secondRoot] = originalPowerj;

  cout << " Over the full grid:" << endl;
  dx = 0.1;
  Real fullHessian = computeMixedBottomHessian(firstRoot, secondRoot);

  Real originalGradient = computeBottomScaledPowerGradient(firstRoot);

  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[secondRoot] = originalPowerj + dx;
    refreshPolynomialCache(_bottomCache, _bottom, _fractal);

    Real newGradient = computeBottomScaledPowerGradient(firstRoot);
    //cout << " dx: " << dx << endl;
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    //cout << " new gradient: " << newGradient << endl;

    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[secondRoot] = originalPowerj;
  refreshPolynomialCache(_bottomCache, _bottom, _fractal);
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial b_i \partial t_j}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianTopBottom(const int topRoot, const int bottomRoot)
{
  int xDebug = 0;
  int yDebug = 0;
  int zDebug = 0;
 
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  VEC3F point = _fractal.cellCenter(xDebug, yDebug, zDebug);
  point[0] = (Real)xDebug / xRes * _lengths[0] + xMin + _center[0];
  point[1] = (Real)yDebug / yRes * _lengths[1] + yMin + _center[1];
  point[2] = (Real)zDebug / zRes * _lengths[2] + zMin + _center[2];

  const QUATERNION iterate(point[0], point[1], point[2], 0);
  const QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  const QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  const QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  Real weight = _weightField(_fractal.cellCenter(xDebug, yDebug, zDebug));
  Real originalScore = weight * (-tanhTerm);

  QUATERNION leftBottom(1,0,0,0);
  QUATERNION rightBottom(1,0,0,0);
  QUATERNION leftTop(1,0,0,0);
  QUATERNION rightTop(1,0,0,0);
  const Real powerScalarBottom = _bottom.powerScalar();
  const Real powerScalarTop    = _top.powerScalar();

  for (int x = 0; x < bottomRoot; x++)
    leftBottom *= (iterate - _bottom.roots()[x]).pow(powerScalarBottom * _bottom.powers()[x]);
  for (int x = bottomRoot + 1; x < _bottom.totalRoots(); x++)
    rightBottom *= (iterate - _bottom.roots()[x]).pow(powerScalarBottom * _bottom.powers()[x]);
  
  for (int x = 0; x < topRoot; x++)
    leftTop *= (iterate - _top.roots()[x]).pow(powerScalarTop * _top.powers()[x]);
  for (int x = topRoot + 1; x < _top.totalRoots(); x++)
    rightTop *= (iterate - _top.roots()[x]).pow(powerScalarTop * _top.powers()[x]);

  const QUATERNION baseBottom = (iterate - _bottom.roots()[bottomRoot]);
  const QUATERNION powerBottom = baseBottom.pow(powerScalarBottom * _bottom.powers()[bottomRoot]);
  const QUATERNION middleBottom = QUATERNION(-1,0,0,0) * (1.0 / powerScalarBottom) * powerBottom * baseBottom.log().inverse();
  const QUATERNION bottomDerivative = leftBottom * middleBottom * rightBottom;
  QUATERNION RderivativeBottom = topEval / bottomDerivative;

  const QUATERNION baseTop = (iterate - _top.roots()[topRoot]);
  const QUATERNION powerTop = baseTop.pow(powerScalarTop * _top.powers()[topRoot]);
  const QUATERNION middleTop = powerScalarTop * powerTop * baseTop.log();
  const QUATERNION topDerivative = leftTop * middleTop * rightTop;
  QUATERNION RderivativeTop = topDerivative / bottomEval;

  Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);

  Real dFdbi = (1.0 / R.dot(R)) * (R.dot(RderivativeBottom));
  Real dSdbi = dSdF * dFdbi;
  cout << " dS / db_i: " << dSdbi << endl;
  
  Real dFdtj = (1.0 / R.dot(R)) * (R.dot(RderivativeTop));
  Real dSdtj = dSdF * dFdtj;
  cout << " dS / dt_j: " << dSdtj << endl;

  Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

  Real dSdcdtj = conformalConformal * dFdtj;
  Real dSdFdtj = dSdcdtj; // i.e. the conformal-top_j derivative

  Real RdotR = R.dot(R);
  Real RdotRSq = RdotR * RdotR;

  QUATERNION R_dbi_dtj = topDerivative / bottomDerivative;

  Real term1 = (RderivativeTop.dot(RderivativeBottom) + R.dot(R_dbi_dtj)) * RdotR;
  Real term2 = R.dot(RderivativeBottom) * (2.0 * R.dot(RderivativeTop));

  Real d2F_dbi_dtj = (term1 - term2) / (RdotRSq);
  Real analyticHessian = dSdFdtj * dFdbi + dSdF * d2F_dbi_dtj;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of mixed top-bottom Hessian, roots " << bottomRoot << ", " << topRoot << endl;
  cout << " ========================================================== " << endl;

  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  Real originalPowerBottom = _bottom.powers()[bottomRoot];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[bottomRoot] = originalPowerBottom + dx;
    QUATERNION bottomEvalNew = _bottom.evaluateScaledPowerFactored(iterate);
    QUATERNION Rnew = (topEval / bottomEvalNew);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScore = weight * (-tanhTermNew);

    Real finiteDiff = (newScore - originalScore) / dx;
    Real relative = (finiteDiff - dSdbi) / dSdbi;
    cout << "analytic: " << dSdbi << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[bottomRoot] = originalPowerBottom;

  cout << " now the Hessian: " << endl;
  Real originalPowerTop = _top.powers()[topRoot];
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[topRoot] = originalPowerTop + dx;
    const QUATERNION topEvalNew    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION Rnew = (topEvalNew / bottomEval);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);

    //Real analyticGradientNew;
    QUATERNION RderivativeNew = topEvalNew / bottomDerivative;
    Real dFdbNew = (1.0 / Rnew.dot(Rnew)) * (Rnew.dot(RderivativeNew));
    Real dSdFNew = -_tanhAlpha * weight * (1.0 - tanhTermNew * tanhTermNew);
    //analyticGradientNew = dSdFNew * dFdbNew;
    Real dSdbNew = dSdFNew * dFdbNew;

    Real finiteDiff = (dSdbNew - dSdbi) / dx;
    Real relative = (finiteDiff - analyticHessian) / analyticHessian;
    cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[topRoot] = originalPowerTop;

  cout << " Over the full grid:" << endl;
  dx = 0.1;
  Real fullHessian = computeTopBottomHessian(topRoot, bottomRoot);

  Real originalGradient = computeBottomScaledPowerGradient(bottomRoot);
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[topRoot] = originalPowerTop + dx;
    refreshPolynomialCache(_topCache, _top, _fractal);

    Real newGradient = computeBottomScaledPowerGradient(bottomRoot);
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    cout  << "analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[topRoot] = originalPowerTop;
  refreshPolynomialCache(_topCache, _top, _fractal);

  cout << " Checking for symmetry:" << endl;
  dx = 0.1;
  originalGradient = computeTopScaledPowerGradient(topRoot);
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[bottomRoot] = originalPowerBottom + dx;
    refreshPolynomialCache(_bottomCache, _bottom, _fractal);

    Real newGradient = computeTopScaledPowerGradient(topRoot);
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    cout  << "analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[bottomRoot] = originalPowerBottom;
  refreshPolynomialCache(_bottomCache, _bottom, _fractal);
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial b_i \partial b_i}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianBottomDiagonal(const int whichRoot)
{
  int xDebug = 0;
  int yDebug = 0;
  int zDebug = 0;
 
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;

  VEC3F point = _fractal.cellCenter(xDebug, yDebug, zDebug);
  point[0] = (Real)xDebug / xRes * _lengths[0] + xMin + _center[0];
  point[1] = (Real)yDebug / yRes * _lengths[1] + yMin + _center[1];
  point[2] = (Real)zDebug / zRes * _lengths[2] + zMin + _center[2];

  const QUATERNION iterate(point[0], point[1], point[2], 0);
  const QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  const QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  const QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  Real weight = _weightField(_fractal.cellCenter(xDebug, yDebug, zDebug));
  Real originalScore = weight * (-tanhTerm);
  cout << " original: " << originalScore << endl;
  cout << " Weight: " << weight << endl;

  Real analyticGradient;
  QUATERNION left(1,0,0,0);
  QUATERNION right(1,0,0,0);
  const Real powerScalar = _bottom.powerScalar();

  for (int x = 0; x < whichRoot; x++)
    left *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);
  for (int x = whichRoot + 1; x < _bottom.totalRoots(); x++)
    right *= (iterate - _bottom.roots()[x]).pow(powerScalar * _bottom.powers()[x]);

  const QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
  const QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
  const QUATERNION middle = (1.0 / powerScalar) * power * base.log().inverse();
  const QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;
  QUATERNION Rderivative = topEval / bottomDerivative;
  Real dFdb = (1.0 / R.dot(R)) * (R.dot(Rderivative));
  Real dSdF = -_tanhAlpha * weight * (1.0 - tanhTerm * tanhTerm);
  analyticGradient = dSdF * dFdb;
  Real dSdb = dSdF * dFdb;
  cout << " analytic gradient: " << analyticGradient << endl;

  Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
  Real dSdcdb= conformalConformal * dFdb;
  Real dSdFdb = dSdcdb;

  Real RdotR = R.dot(R);
  Real RdotRSq = RdotR * RdotR;

  QUATERNION logTermInverse = base.log().inverse();
  QUATERNION middle2 = 1.0 / (powerScalar * powerScalar) * power * logTermInverse * logTermInverse;
  QUATERNION Rderivative2 = topEval / (left * middle2 * right);

  Real term1 = (Rderivative.dot(Rderivative) + R.dot(Rderivative2)) * RdotR;
  Real term2 = R.dot(Rderivative) * (2.0 * R.dot(Rderivative));

  Real d2Fdb2 = (term1 - term2) / (RdotRSq);
  Real analyticHessian = dSdFdb * dFdb + dSdF * d2Fdb2;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of bottom-bottom Hessian, root " << whichRoot << endl;
  cout << " ========================================================== " << endl;

  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  Real originalPower = _bottom.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] = originalPower + dx;
    QUATERNION bottomEvalNew    = _bottom.evaluateScaledPowerFactored(iterate);
    QUATERNION Rnew = (topEval / bottomEvalNew);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScore = weight * (-tanhTermNew);

    Real finiteDiff = (newScore - originalScore) / dx;
    Real relative = (finiteDiff - analyticGradient) / analyticGradient;
    cout << "analytic: " << analyticGradient << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  cout << " now the Hessian: " << endl;
  _bottom.powersMutable()[whichRoot] = originalPower;
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] = originalPower + dx;
    const QUATERNION topEvalNew    = _top.evaluateScaledPowerFactored(iterate);
    const QUATERNION bottomEvalNew = _bottom.evaluateScaledPowerFactored(iterate);
    const QUATERNION Rnew = (topEvalNew / bottomEvalNew);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);

    //Real analyticGradientNew;

    QUATERNION baseNew = (iterate - _bottom.roots()[whichRoot]);
    QUATERNION powerNew = baseNew.pow(powerScalar * _bottom.powers()[whichRoot]);
    QUATERNION middleNew = (1.0 / powerScalar) * powerNew * baseNew.log().inverse();
    QUATERNION bottomDerivativeNew = QUATERNION(-1,0,0,0) * left * middleNew * right;
    QUATERNION RderivativeNew = topEval / bottomDerivativeNew;
    Real dFdbNew = (1.0 / Rnew.dot(Rnew)) * (Rnew.dot(RderivativeNew));
    Real dSdFNew = -_tanhAlpha * weight * (1 - tanhTermNew * tanhTermNew);
    //analyticGradientNew = dSdFNew * dFdbNew;
    Real dSdbNew = dSdFNew * dFdbNew;
    Real finiteDiff = (dSdbNew - dSdb) / dx;
    Real relative = (finiteDiff - analyticHessian) / analyticHessian;
    cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[whichRoot] = originalPower;

  cout << " Over the full grid:" << endl;
  dx = 0.1;
  Real fullHessian = computeBottomDiagonalHessian(whichRoot);
  Real originalGradient = computeBottomScaledPowerGradient(whichRoot);

  cout << " original: " << originalGradient << " original single sample: " << analyticGradient << " diff: " << originalGradient - analyticGradient << endl;
  for (int x = 0; x < 10; x++)
  {
    _bottom.powersMutable()[whichRoot] = originalPower + dx;
    refreshPolynomialCache(_bottomCache, _bottom, _fractal);    

    Real newGradient = computeBottomScaledPowerGradient(whichRoot);
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _bottom.powersMutable()[whichRoot] = originalPower;
  refreshPolynomialCache(_bottomCache, _bottom, _fractal);
}

///////////////////////////////////////////////////////////////////////
// debug the Hessian entry for \frac{\partial S}{\partial c \partial t_i}
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugHessianTopConformal(const int whichRoot)
{
  //int xDebug = 80;
  //int yDebug = 54;
  //int zDebug = 48;
  int xDebug = 0;
  int yDebug = 0;
  int zDebug = 0;
  //const VEC3F point = _weightField.cellCenter(xDebug, yDebug, zDebug);
  const VEC3F point = _fractal.cellCenter(xDebug, yDebug, zDebug);
  const QUATERNION iterate(point[0], point[1], point[2], 0);
  const QUATERNION topEval    = _top.evaluatePowerFactored(iterate);
  const QUATERNION bottomEval = _bottom.evaluatePowerFactored(iterate);
  const QUATERNION R = (topEval / bottomEval);

  Real F = log(R.magnitude()) + log(_expScaling);
  Real interior = _tanhAlpha * (F - _tanhThreshold);
  Real tanhTerm = tanh(interior);
  Real weight = _weightField(point);
  Real originalScore = weight * (-tanhTerm);
  cout << " original: " << originalScore << endl;
  cout << " Weight: " << weight << endl;

  Real analyticGradient;
  QUATERNION left(1,0,0,0);
  QUATERNION right(1,0,0,0);
  const Real powerScalar = _top.powerScalar();

  for (int x = 0; x < whichRoot; x++)
    left *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);
  for (int x = whichRoot + 1; x < _top.totalRoots(); x++)
    right *= (iterate - _top.roots()[x]).pow(powerScalar * _top.powers()[x]);

  QUATERNION base = (iterate - _top.roots()[whichRoot]);
  QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
  QUATERNION middle = power * base.log();
  QUATERNION topDerivative = left * middle * right;
  QUATERNION Rderivative = (topDerivative) / bottomEval;
  Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
  Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
  analyticGradient = dSdF * dFdt;
  Real dSdt = dSdF * dFdt;
  cout << " analytic gradient: " << analyticGradient << endl;

  Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  Real conformalConformal = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;
  Real analyticHessian = conformalConformal * dFdt;

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of top-conformal Hessian, root " << whichRoot << endl;
  cout << " ========================================================== " << endl;

  cout << " First the gradient: " << endl;
  Real dx = 0.1;
  Real originalPower = _top.powers()[whichRoot];
  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] = originalPower + dx;
    QUATERNION topEvalNew    = _top.evaluatePowerFactored(iterate);
    QUATERNION Rnew = (topEvalNew / bottomEval);

    Real Fnew = log(Rnew.magnitude()) + log(_expScaling);
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real newScore = weight * (-tanhTermNew);

    Real finiteDiff = (newScore - originalScore) / dx;
    Real relative = (finiteDiff - analyticGradient) / analyticGradient;
    cout << "analytic: " << analyticGradient << " finite: "<< finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  cout << " now the Hessian: " << endl;
  _top.powersMutable()[whichRoot] = originalPower;
  dx = 0.1;
  for (int x = 0; x < 10; x++)
  {
    Real Fnew = log(R.magnitude()) + log(_expScaling) + dx;
    Real interiorNew = _tanhAlpha * (Fnew - _tanhThreshold);
    Real tanhTermNew = tanh(interiorNew);
    Real dSdFnew = -_tanhAlpha * weight * (1 - tanhTermNew * tanhTermNew);
    Real dSdtNew = dSdFnew * dFdt;

    Real finiteDiff = (dSdtNew - dSdt) / dx;

    Real relative = (finiteDiff - analyticHessian) / analyticHessian;

    cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }

  cout << " Over the full grid:" << endl;
  dx = 0.1;
  Real fullHessian = computeTopConformalHessian(whichRoot);
  
  Real originalGradient = computeLogPowerScaledRationalScoreExpGradient();

  for (int x = 0; x < 10; x++)
  {
    _top.powersMutable()[whichRoot] = originalPower + dx;
    Real newGradient = computeLogPowerScaledRationalScoreExpGradient();
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    //cout << "analytic: " << analyticHessian << " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[whichRoot] = originalPower;

  cout << " Checking for Hessian symmetry:" << endl;
  dx = 0.1;
  originalGradient = computeTopScaledPowerGradient(whichRoot);
  const Real backup = _expScaling;
  for (int x = 0; x < 10; x++)
  {
    _expScaling = backup + dx;

    Real newGradient = computeTopScaledPowerGradient(whichRoot);
    Real finiteDiff = (newGradient - originalGradient) / dx;
    Real relative = (fullHessian - finiteDiff) / fullHessian;

    cout  <<" analytic: " << fullHessian <<  " finite diff: " << finiteDiff << " relative: " << relative << " dx: " << dx << endl;
    dx *= 0.1;
  }
  _top.powersMutable()[whichRoot] = originalPower;
}

///////////////////////////////////////////////////////////////////////
// debug the derivative of a fractal conformal radius
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::debugConformalRadiusDerivative()
{
  // compute the summed score via function
  Real groundScore = computeLogPowerRationalScore();
  FIELD_3D groundFractal = _fractal;

  // push the ground fractal through the tanh function
  for (int x = 0; x < _fractal.totalCells(); x++)
  {
    // do a more robust, resolution-independent version here
    VEC3F cellPosition = _fractal.cellCenter(x);
    Real inverseDistance = _inverseDistanceField(cellPosition);

    Real tanhVersionBefore = -tanh(_tanhAlpha * (_fractal[x]));
    Real tanhVersionAfter = tanhVersionBefore * inverseDistance;

    groundFractal[x] = tanhVersionAfter;
  }

  Real conformalGradient = inverseDistanceScoreConformalGradient(_fractal, _inverseDistanceField, _inverseFieldSum, 0);

  cout << " ========================================================== " << endl;
  cout << "  Convergence plot of summed score, CONFORMAL RADIUS derivative " << endl;
  cout << " ========================================================== " << endl;

  cout << " Conformal radius gradient: " << conformalGradient << endl;
  
  Real dx = 0.1;
  Real cOriginal = log(_expScaling);
  for (int x = 0; x < 10; x++)
  {
    _expScaling = exp(cOriginal + dx);
    Real newScore = computeLogPowerRationalScore();
    _expScaling = exp(cOriginal);
    
    Real numericalGradient = (newScore - groundScore) / dx;
    Real diff = (numericalGradient - conformalGradient) / conformalGradient;

    cout << " diff: " << diff << "\t dx: " << dx << "\t numerical grad: " << numericalGradient << "\t Snew: " << newScore << "\t S: " << groundScore << endl;
    dx *= 0.1;
  }
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopPowerGradient(const int whichRoot) const
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        //Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

#if 0
        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
#else
        // VERIFIED
        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
#if EXPLICIT_LEADING_ROOT
        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNION(1,0,0,0);
#else
        QUATERNION topEval = iterate * topForward.back();
        QUATERNION bottomEval = (bottomForward.size() > 0) ? iterate * bottomForward.back() : QUATERNION(1,0,0,0);
#endif
#endif

        QUATERNION R = (topEval / bottomEval);
        Real F = log(R.magnitude()) + log(_expScaling);

#if 0
        QUATERNION topDerivative = _top.powerDerivative(iterate, whichRoot);
#else
        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_top.totalRoots() > 1)
        {  
          // VERIFIED
          //for (int i = 0; i < whichRoot; i++)
          //  left *= (iterate - _top.roots()[i]).pow(_top.powers()[i]);
          if (whichRoot > 0)
            left = topForward[whichRoot - 1];
         
          // VERIFIED
          //for (int i = whichRoot + 1; i < _top.totalRoots(); i++)
          //  right *= (iterate - _top.roots()[i]).pow(_top.powers()[i]);
          if (whichRoot < _top.totalRoots() - 1)
            right = _topCache.backward(x,y,z)[whichRoot + 1];
        }

        // get the derivative of the root in question
        QUATERNION base = (iterate - _top.roots()[whichRoot]);
        QUATERNION power = base.pow(_top.powers()[whichRoot]);
        QUATERNION middle = power * base.log();

        QUATERNION topDerivative = left * middle * right;
#endif

#if EXPLICIT_LEADING_ROOT
        QUATERNION Rderivative = (topDerivative) / bottomEval;
#else
        QUATERNION Rderivative = (iterate * topDerivative) / bottomEval;
#endif
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);

        // TODO: REPLACE HERE
#if USING_WEIGHT
        Real weight = _weightField(cellPosition);
#else
        Real inverseDistance = _inverseDistanceField(cellPosition);
#endif
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
#if USING_WEIGHT
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
#else
        Real dSdF = -_tanhAlpha * inverseDistance * (1 - tanhTerm * tanhTerm); 
#endif
        Real dSdt = dSdF * dFdt;

        final += dSdt;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopScaledPowerGradientDebug(const int whichRoot) const
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
    
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        // VERIFIED
        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);

        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNION(1,0,0,0);

        QUATERNION R = (topEval / bottomEval);

        // should probably track how often this happens, because it suggests we've run out of precision
        Real Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        //if (Rmagnitude <= 0.0) continue;
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        Real F = log(Rmagnitude) + log(_expScaling);

        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_top.totalRoots() > 1)
        { 
          // VERIFIED
          //for (int i = 0; i < whichRoot; i++)
          //  left *= (iterate - _top.roots()[i]).pow(_top.powerScalar() * _top.powers()[i]);
          if (whichRoot > 0)
            left = topForward[whichRoot - 1];
         
          // VERIFIED
          //for (int i = whichRoot + 1; i < _top.totalRoots(); i++)
          //  right *= (iterate - _top.roots()[i]).pow(_top.powerScalar() * _top.powers()[i]);
          if (whichRoot < _top.totalRoots() - 1)
            right = _topCache.backward(x,y,z)[whichRoot + 1];
        }

        // get the derivative of the root in question
        const Real powerScalar = _top.powerScalar();
        QUATERNION base = (iterate - _top.roots()[whichRoot]);
        //QUATERNION power = base.pow(_top.powers()[whichRoot]);
        //QUATERNION middle = power * base.log();
        QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
        QUATERNION middle = powerScalar * power * base.log();

        QUATERNION topDerivative = left * middle * right;

        QUATERNION Rderivative = topDerivative / bottomEval;
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);

        Real weight = _weightField(cellPosition);
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
        Real dSdt = dSdF * dFdt;
        if (isnan(dSdt)) continue;

        final += dSdt;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopGradientAP(const int whichRoot)
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
    
        Real magnitude = iterate.magnitude();
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        // VERIFIED
        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);

        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNION(1,0,0,0);
        QUATERNION R = (topEval / bottomEval);

        // should probably track how often this happens, because it suggests we've run out of precision
        Real Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        Real F = log(Rmagnitude) + log(_expScaling);

        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_top.totalRoots() > 1)
        { 
          if (whichRoot > 0)
            left = topForward[whichRoot - 1];
         
          if (whichRoot < _top.totalRoots() - 1)
            right = _topCache.backward(x,y,z)[whichRoot + 1];
        }

        // get the derivative of the root in question
        const Real powerScalar = _top.powerScalar();
        QUATERNION base = (iterate - _top.roots()[whichRoot]);
        QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
        QUATERNION middle = powerScalar * power * base.log();

        QUATERNION topDerivative = left * middle * right;

        QUATERNION Rderivative = topDerivative / bottomEval;
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);

        Real weight = _weightField(cellPosition);
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
        Real dSdt = dSdF * dFdt;

        if (isnan(dSdt))
        {
          _nans(x,y,z) = 1.0;
          continue;
        }
        if (isinf(dSdt))
        {
          _infs(x,y,z) = 1.0;
          continue;
        }
        
        final += dSdt;
      }

  int problems = _overflowed.sum() + _nans.sum() + _infs.sum();
  if (problems == 0) return final;
  _adaptiveFired = true;
  
  // do an overflow pass
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        if (_overflowed(x,y,z) + _nans(x,y,z) + _infs(x,y,z) < 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        long double xReal = point[0];
        long double yReal = point[1];
        long double zReal = point[2];
#else
        long double xReal = (long double)x / xRes * _lengths[0] + xMin + _center[0];
        long double yReal = (long double)y / yRes * _lengths[1] + yMin + _center[1];
        long double zReal = (long double)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNIONE iterate(xReal, yReal, zReal, _quaternionSlice);
    
        long double magnitude = iterate.magnitude();
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        // make a copy, but don't write it back
        vector<QUATERNIONE> forwardCacheTop = _topCacheExtended.forward(x,y,z);
        vector<QUATERNIONE> backwardCacheTop = _topCacheExtended.backward(x,y,z);
        vector<QUATERNIONE> forwardCacheBottom = _bottomCacheExtended.forward(x,y,z);
        vector<QUATERNIONE> backwardCacheBottom = _bottomCacheExtended.backward(x,y,z);

        // check that the cache was built
        if (_nans(x,y,z) > 0.5 || _infs(x,y,z) > 0.5)
        {
          if (forwardCacheTop.size() == 0)
            QUATERNIONE topEval = _topExtended.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);

          //if (_bottomCacheExtended.forward(x,y,z).size() == 0)
          if (forwardCacheBottom.size() == 0)
            QUATERNIONE bottomEval = _bottomExtended.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
        }

        //const vector<QUATERNIONE>& topForward = _topCacheExtended.forward(x,y,z);
        //const vector<QUATERNIONE>& bottomForward = _bottomCacheExtended.forward(x,y,z);
        const vector<QUATERNIONE>& topForward = forwardCacheTop;
        const vector<QUATERNIONE>& bottomForward = forwardCacheBottom;

        QUATERNIONE topEval = topForward.back();
        QUATERNIONE bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNIONE(1,0,0,0);
        QUATERNIONE R = (topEval / bottomEval);

        // should probably track how often this happens, because it suggests we've run out of precision
        long double Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        long double F = log(Rmagnitude) + log(_expScaling);

        QUATERNIONE left(1,0,0,0);
        QUATERNIONE right(1,0,0,0);
      
        if (_top.totalRoots() > 1)
        { 
          if (whichRoot > 0)
            left = topForward[whichRoot - 1];
         
          if (whichRoot < _top.totalRoots() - 1)
            //right = _topCacheExtended.backward(x,y,z)[whichRoot + 1];
            right = backwardCacheTop[whichRoot + 1];
        }

        // get the derivative of the root in question
        const long double powerScalar = _topExtended.powerScalar();
        QUATERNIONE base = (iterate - _topExtended.roots()[whichRoot]);
        QUATERNIONE power = base.pow(powerScalar * _topExtended.powers()[whichRoot]);
        QUATERNIONE middle = powerScalar * power * base.log();

        QUATERNIONE topDerivative = left * middle * right;

        QUATERNIONE Rderivative = topDerivative / bottomEval;
        long double dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);

        long double weight = _weightField(cellPosition);
  
        long double tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
        long double dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
        long double dSdt = dSdF * dFdt;

        final += dSdt;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopScaledPowerGradient(const int whichRoot) const
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
    
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

#if 0
        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
#else
        // VERIFIED
        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);

// NEW: IMPROVES CONVERGENCE?
#if EXPLICIT_LEADING_ROOT
        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNION(1,0,0,0);
#else
        QUATERNION topEval = iterate * topForward.back();
        QUATERNION bottomEval = (bottomForward.size() > 0) ? iterate * bottomForward.back() : QUATERNION(1,0,0,0);
#endif
#endif

        QUATERNION R = (topEval / bottomEval);

        // should probably track how often this happens, because it suggests we've run out of precision
        Real Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        //if (Rmagnitude <= 0.0) continue;
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        Real F = log(Rmagnitude) + log(_expScaling);
        //Real F = log(R.magnitude()) + log(_expScaling);

#if 0
        QUATERNION topDerivative = _top.powerDerivative(iterate, whichRoot);
#else
        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_top.totalRoots() > 1)
        { 
          // VERIFIED
          //for (int i = 0; i < whichRoot; i++)
          //  left *= (iterate - _top.roots()[i]).pow(_top.powerScalar() * _top.powers()[i]);
          if (whichRoot > 0)
            left = topForward[whichRoot - 1];
         
          // VERIFIED
          //for (int i = whichRoot + 1; i < _top.totalRoots(); i++)
          //  right *= (iterate - _top.roots()[i]).pow(_top.powerScalar() * _top.powers()[i]);
          if (whichRoot < _top.totalRoots() - 1)
            right = _topCache.backward(x,y,z)[whichRoot + 1];
        }

        // get the derivative of the root in question
        const Real powerScalar = _top.powerScalar();
        QUATERNION base = (iterate - _top.roots()[whichRoot]);
        //QUATERNION power = base.pow(_top.powers()[whichRoot]);
        //QUATERNION middle = power * base.log();
        QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
        QUATERNION middle = powerScalar * power * base.log();

        QUATERNION topDerivative = left * middle * right;
#endif

// NEW: IMPROVES CONVERGENCE?
#if EXPLICIT_LEADING_ROOT
        QUATERNION Rderivative = topDerivative / bottomEval;
#else
        QUATERNION Rderivative = (iterate * topDerivative) / bottomEval;
#endif
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);

        // TODO: REPLACE HERE
#if USING_WEIGHT
        Real weight = _weightField(cellPosition);
#else
        Real inverseDistance = _inverseDistanceField(cellPosition);
#endif
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
#if USING_WEIGHT
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
#else
        Real dSdF = -_tanhAlpha * inverseDistance * (1 - tanhTerm * tanhTerm); 
#endif
        Real dSdt = dSdF * dFdt;
        if (isnan(dSdt)) continue;

        final += dSdt;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopBulkScalarGradientAP()
{
  TIMER functionTimer(__FUNCTION__);
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  vector<Real> sums(zRes);

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
  {
    Real& threadSum = sums[z];
    threadSum = 0;
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        Real magnitude = iterate.magnitude();

        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        for (int i = 0; i < _top.totalRoots(); i++)
        {
          int whichRoot = i;
          const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
          const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
          QUATERNION topEval = topForward.back();
          QUATERNION bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNION(1,0,0,0);

          QUATERNION R = (topEval / bottomEval);
          
          // should probably track how often this happens, because it suggests we've run out of precision
          Real Rmagnitude = R.magnitude();
          
          // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
          if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;
          
          Real F = log(Rmagnitude) + log(_expScaling);

          QUATERNION left(1,0,0,0);
          QUATERNION right(1,0,0,0);
        
          if (_top.totalRoots() > 1)
          {  
            if (whichRoot > 0)
              left = topForward[whichRoot - 1];
           
            if (whichRoot < _top.totalRoots() - 1)
              right = _topCache.backward(x,y,z)[whichRoot + 1];
          }

          // get the derivative of the root in question
          const Real powerScalar = _top.powerScalar();
          QUATERNION base = (iterate - _top.roots()[whichRoot]);
          QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
          // this is the only switch: powerScalar -> _top.powers()[whichRoot]
          //  (remember we're in the BULK gradient, not the individual one.)
          //QUATERNION middle = powerScalar * power * base.log();
          QUATERNION middle = _top.powers()[whichRoot] * power * base.log();

          QUATERNION topDerivative = left * middle * right;

          QUATERNION Rderivative = (topDerivative) / bottomEval;
          Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);

          Real weight = _weightField(cellPosition);
    
          Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          Real dSdt = dSdF * dFdt;

          if (isnan(dSdt))
          {
            _nans(x,y,z) = 1.0;
            continue;
          }
          if (isinf(dSdt))
          {
            _infs(x,y,z) = 1.0;
            continue;
          }

          threadSum += dSdt;
        }
      }
  }

  // if there were no issues, all done
  int problems = _overflowed.sum() + _nans.sum() + _infs.sum();
  if (problems == 0)
  { 
    Real final = 0;
    for (int z = 0; z < zRes; z++)
      final += sums[z];

    return final;
  }
  _adaptiveFired = true;

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
  {
    Real& threadSum = sums[z];
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        if (_overflowed(x,y,z) + _nans(x,y,z) + _infs(x,y,z) < 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        long double xReal = point[0];
        long double yReal = point[1];
        long double zReal = point[2];
#else
        long double xReal = (long double)x / xRes * _lengths[0] + xMin + _center[0];
        long double yReal = (long double)y / yRes * _lengths[1] + yMin + _center[1];
        long double zReal = (long double)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNIONE iterate(xReal, yReal, zReal, _quaternionSlice);
        long double magnitude = iterate.magnitude();

        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        // check that the cache was built
        if (_nans(x,y,z) > 0.5 || _infs(x,y,z) > 0.5)
        {
          if (_topCacheExtended.forward(x,y,z).size() == 0)
          {
            vector<QUATERNIONE>& forwardCacheTop = _topCacheExtended.forward(x,y,z);
            vector<QUATERNIONE>& backwardCacheTop = _topCacheExtended.backward(x,y,z);
            QUATERNIONE topEval = _topExtended.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
          }
          if (_bottomCacheExtended.forward(x,y,z).size() == 0)
          {
            vector<QUATERNIONE>& forwardCacheBottom = _bottomCacheExtended.forward(x,y,z);
            vector<QUATERNIONE>& backwardCacheBottom = _bottomCacheExtended.backward(x,y,z);
            QUATERNIONE bottomEval = _bottomExtended.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
          }
        }

        for (int i = 0; i < _top.totalRoots(); i++)
        {
          int whichRoot = i;
          const vector<QUATERNIONE>& topForward = _topCacheExtended.forward(x,y,z);
          const vector<QUATERNIONE>& bottomForward = _bottomCacheExtended.forward(x,y,z);
          QUATERNIONE topEval = topForward.back();
          QUATERNIONE bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNIONE(1,0,0,0);

          QUATERNIONE R = (topEval / bottomEval);
          
          // should probably track how often this happens, because it suggests we've run out of precision
          long double Rmagnitude = R.magnitude();
          
          // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
          if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;
          
          long double F = log(Rmagnitude) + log(_expScaling);

          QUATERNIONE left(1,0,0,0);
          QUATERNIONE right(1,0,0,0);
        
          if (_top.totalRoots() > 1)
          {  
            if (whichRoot > 0)
              left = topForward[whichRoot - 1];
           
            if (whichRoot < _top.totalRoots() - 1)
              right = _topCacheExtended.backward(x,y,z)[whichRoot + 1];
          }

          // get the derivative of the root in question
          const long double powerScalar = _topExtended.powerScalar();
          QUATERNIONE base = (iterate - _topExtended.roots()[whichRoot]);
          QUATERNIONE power = base.pow(powerScalar * _topExtended.powers()[whichRoot]);
          QUATERNIONE middle = _topExtended.powers()[whichRoot] * power * base.log();

          QUATERNIONE topDerivative = left * middle * right;

          QUATERNIONE Rderivative = (topDerivative) / bottomEval;
          long double dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);

          long double weight = _weightField(cellPosition);
    
          long double tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          long double dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          long double dSdt = dSdF * dFdt;

          threadSum += dSdt;
        }
      }
  }

  Real final = 0;
  for (int z = 0; z < zRes; z++)
    final += sums[z];

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopScaledPowerGradient() const
{
  TIMER functionTimer(__FUNCTION__);
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  vector<Real> sums(zRes);

#if OMP_ENABLED
#pragma omp parallel
#pragma omp for schedule(dynamic)
#endif
  for (int z = 0; z < zRes; z++)
  {
    Real& threadSum = sums[z];
    threadSum = 0;
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();

        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        //int totalIterations = 0;

        for (int i = 0; i < _top.totalRoots(); i++)
        {
          int whichRoot = i;
          const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
          const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
#if EXPLICIT_LEADING_ROOT
          QUATERNION topEval = topForward.back();
          QUATERNION bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNION(1,0,0,0);
#else
          QUATERNION topEval = iterate * topForward.back();
          QUATERNION bottomEval = (bottomForward.size() > 0) ? iterate * bottomForward.back() : QUATERNION(1,0,0,0);
#endif

          QUATERNION R = (topEval / bottomEval);
          //Real F = log(R.magnitude()) + log(_expScaling);
          
          // should probably track how often this happens, because it suggests we've run out of precision
          Real Rmagnitude = R.magnitude();
          
          // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
          //if (Rmagnitude <= 0.0) continue;
          if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;
          
          Real F = log(Rmagnitude) + log(_expScaling);

          QUATERNION left(1,0,0,0);
          QUATERNION right(1,0,0,0);
        
          if (_top.totalRoots() > 1)
          {  
            if (whichRoot > 0)
              left = topForward[whichRoot - 1];
           
            if (whichRoot < _top.totalRoots() - 1)
              right = _topCache.backward(x,y,z)[whichRoot + 1];
          }

          // get the derivative of the root in question
          const Real powerScalar = _top.powerScalar();
          QUATERNION base = (iterate - _top.roots()[whichRoot]);
          QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
          // this is the only switch: powerScalar -> _top.powers()[whichRoot]
          //  (remember we're in the BULK gradient, not the individual one.)
          //QUATERNION middle = powerScalar * power * base.log();
          QUATERNION middle = _top.powers()[whichRoot] * power * base.log();

          QUATERNION topDerivative = left * middle * right;

#if EXPLICIT_LEADING_ROOT
          QUATERNION Rderivative = (topDerivative) / bottomEval;
#else
          QUATERNION Rderivative = (iterate * topDerivative) / bottomEval;
#endif
          Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);

          Real weight = _weightField(cellPosition);
    
          Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          Real dSdt = dSdF * dFdt;
          if (isnan(dSdt)) continue;

          threadSum += dSdt;
        }
      }
// inconsistent sum ordering here creates problems
//#pragma omp atomic
//    final += threadSum;
  }
  
  Real final = 0;
  for (int z = 0; z < zRes; z++)
    final += sums[z];

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopScaledPowerGradientDebug() const
{
  TIMER functionTimer(__FUNCTION__);
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;

  vector<Real> sums(zRes);

  for (int z = 0; z < zRes; z++)
  {
    Real& threadSum = sums[z];
    threadSum = 0;
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);

        for (int i = 0; i < _top.totalRoots(); i++)
        {
          int whichRoot = i;
          const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
          const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
          QUATERNION topEval = topForward.back();
          QUATERNION bottomEval = (bottomForward.size() > 0) ? bottomForward.back() : QUATERNION(1,0,0,0);

          QUATERNION R = (topEval / bottomEval);
          //Real F = log(R.magnitude()) + log(_expScaling);
          
          // should probably track how often this happens, because it suggests we've run out of precision
          Real Rmagnitude = R.magnitude();

          // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
          //if (Rmagnitude <= 0.0) continue;
          if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

          Real F = log(Rmagnitude) + log(_expScaling);

          QUATERNION left(1,0,0,0);
          QUATERNION right(1,0,0,0);
        
          if (_top.totalRoots() > 1)
          {  
            if (whichRoot > 0)
              left = topForward[whichRoot - 1];
           
            if (whichRoot < _top.totalRoots() - 1)
              right = _topCache.backward(x,y,z)[whichRoot + 1];
          }

          // get the derivative of the root in question
          const Real powerScalar = _top.powerScalar();
          QUATERNION base = (iterate - _top.roots()[whichRoot]);
          QUATERNION power = base.pow(powerScalar * _top.powers()[whichRoot]);
          QUATERNION middle = _top.powers()[whichRoot] * power * base.log();

          QUATERNION topDerivative = left * middle * right;

          QUATERNION Rderivative = (topDerivative) / bottomEval;
          Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);

          Real weight = _weightField(cellPosition);
    
          Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          Real dSdt = dSdF * dFdt;

          if (isnan(dSdt))
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << " NaN found at: " << x << " " << y << " " << z << endl;
            cout << " top cache dims: " << _topCache.dims() << endl;
            cout << " bottom cache dims: " << _bottomCache.dims() << endl;
            cout << " iterate: " << iterate << endl;
            cout << " dSdt: " << dSdt << endl;
            cout << " dSdF: " << dSdF << endl;
            cout << " dFdt: " << dFdt << endl;
            cout << " tanhTerm: " << tanhTerm << endl;
            cout << " R: " << R << endl;
            cout << " topEval: " << topEval << endl;
            cout << " bottomEval: " << bottomEval << endl;
            cout << " bottom inverse: " << bottomEval.inverse() << endl;
            cout << " F: " << F << endl;
            cout << " left: " << left << endl;
            cout << " middle: " << middle << endl;
            cout << " right: " << right << endl;
            cout << " right inverse: " << right.inverse() << endl;

            cout << " Top multiply series: " << flush;
            for (unsigned int d = 0; d < topForward.size(); d++)
              cout << topForward[d] << " " << flush;
            cout << endl;
            
            cout << " Bottom multiply series: " << flush;
            for (unsigned int d = 0; d < bottomForward.size(); d++)
              cout << bottomForward[d] << " " << flush;
            cout << endl;

            printRootPowers();

            cout << " Trying to reproduce first top Nan: " << endl;
            QUATERNION first = iterate - _top.roots()[0];
            first = first.pow(_top.powerScalar() * _top.powers()[0]);
            cout << first << endl;

            vector<QUATERNION> forward;
            vector<QUATERNION> backward;
            QUATERNION eval = _top.evaluateScaledPowerFactored(iterate,forward, backward);
            cout << " redone eval: " << eval << endl;
            cout << " Redone forward cache: " << endl;
            for (unsigned int d = 0; d < forward.size(); d++)
              cout << forward[d] << " " << flush;
            cout << endl;
            exit(0);
          }
          if (isnan(dSdt)) continue;

          threadSum += dSdt;
        }
      }
  }
 
  cout << " Slices: " << flush;
  Real final = 0;
  for (int z = 0; z < zRes; z++)
  {
    final += sums[z];
    cout << sums[z] << " " << flush;
  }
  cout << endl;

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomPowerGradient(const int whichRoot) const
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        //Real magnitude = iterate.magnitude();
        //int totalIterations = 0;

#if 0
        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
#else        
        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
#if EXPLICIT_LEADING_ROOT
        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = bottomForward.back();
#else
        QUATERNION topEval = iterate * topForward.back();
        QUATERNION bottomEval = iterate * bottomForward.back();
#endif
#endif
        QUATERNION R = (topEval / bottomEval);
        Real F = log(R.magnitude()) + log(_expScaling);

#if 0
        QUATERNION bottomDerivative = _bottom.inversePowerDerivative(iterate, whichRoot);
#else
        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_bottom.totalRoots() > 1)
        {
          // VERIFIED
          if (whichRoot > 0)
            left = bottomForward[whichRoot - 1];
          
          if (whichRoot < _bottom.totalRoots() - 1)
            right = _bottomCache.backward(x,y,z)[whichRoot + 1];
          /*
          // TODO: this can be optimized using the cache, but need more roots
          // to properly debug
          for (int i = 0; i < whichRoot; i++)
            left *= (iterate - _bottom.roots()[i]).pow(_bottom.powers()[i]);
          
          for (int i = whichRoot + 1; i < _bottom.totalRoots(); i++)
            right *= (iterate - _bottom.roots()[i]).pow(_bottom.powers()[i]);
            */
        }

        // get the derivative of the root in question
        QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
        QUATERNION power = base.pow(_bottom.powers()[whichRoot]);

        // DEBUG: where does the inverse come from?
        QUATERNION middle = power * base.log().inverse();
        //QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;
        //QUATERNION middle = power * base.log();
        QUATERNION bottomDerivative = left * middle * right;
#endif

        //QUATERNION Rderivative = topEval / (iterate * bottomDerivative);
#if EXPLICIT_LEADING_ROOT
        QUATERNION Rderivative = QUATERNION(-1,0,0,0) * (topEval / (bottomDerivative));
#else
        QUATERNION Rderivative = QUATERNION(-1,0,0,0) * (topEval / (iterate * bottomDerivative));
#endif
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
#if USING_WEIGHT
        Real weight = _weightField(cellPosition);
#else
        Real inverseDistance = _inverseDistanceField(cellPosition);
#endif
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
#if USING_WEIGHT
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
#else
        Real dSdF = -_tanhAlpha * inverseDistance * (1 - tanhTerm * tanhTerm); 
#endif
        Real dSdt = dSdF * dFdt;

        final += dSdt;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomScaledPowerGradientDebug(const int whichRoot) const
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = bottomForward.back();
        
        QUATERNION R = (topEval / bottomEval);

        // should probably track how often this happens, because it suggests we've run out of precision
        Real Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        //if (Rmagnitude <= 0.0) continue;
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        Real F = log(Rmagnitude) + log(_expScaling);

        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_bottom.totalRoots() > 1)
        {
          // VERIFIED
          if (whichRoot > 0)
            left = bottomForward[whichRoot - 1];
          
          if (whichRoot < _bottom.totalRoots() - 1)
            right = _bottomCache.backward(x,y,z)[whichRoot + 1];
        }

        // get the derivative of the root in question
        const Real powerScalar = _bottom.powerScalar();
        QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
        QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
        QUATERNION middle = (1.0 / powerScalar) * power * base.log().inverse();

        QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;

        QUATERNION Rderivative = topEval / bottomDerivative;
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real weight = _weightField(cellPosition);
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
        Real dSdt = dSdF * dFdt;

        if (isnan(dSdt))
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " NaN found at: " << x << " " << y << " " << z << endl;
          cout << " top cache dims: " << _topCache.dims() << endl;
          cout << " bottom cache dims: " << _bottomCache.dims() << endl;
          cout << " iterate: " << iterate << endl;
          cout << " dSdt: " << dSdt << endl;
          cout << " dSdF: " << dSdF << endl;
          cout << " dFdt: " << dFdt << endl;
          cout << " tanhTerm: " << tanhTerm << endl;
          cout << " R: " << R << endl;
          cout << " R derivative: " << Rderivative << endl;
          cout << " R dot R: " << R.dot(R) << endl;
          cout << " R magnitude: " << R.magnitude() << endl;
          cout << " topEval: " << topEval << endl;
          cout << " bottomEval: " << bottomEval << endl;
          cout << " bottom inverse: " << bottomEval.inverse() << endl;
          cout << " F: " << F << endl;
          cout << " left: " << left << endl;
          cout << " middle: " << middle << endl;
          cout << " right: " << right << endl;
          cout << " right inverse: " << right.inverse() << endl;

          return -1;
        }
        final += dSdt;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomGradientAP(const int whichRoot)
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        Real magnitude = iterate.magnitude();
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = bottomForward.back();
        QUATERNION R = (topEval / bottomEval);
        //Real F = log(R.magnitude()) + log(_expScaling);

        // should probably track how often this happens, because it suggests we've run out of precision
        Real Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        //if (Rmagnitude <= 0.0) continue;
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        Real F = log(Rmagnitude) + log(_expScaling);

        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_bottom.totalRoots() > 1)
        {
          // VERIFIED
          if (whichRoot > 0)
            left = bottomForward[whichRoot - 1];
          
          if (whichRoot < _bottom.totalRoots() - 1)
            right = _bottomCache.backward(x,y,z)[whichRoot + 1];
        }

        // get the derivative of the root in question
        const Real powerScalar = _bottom.powerScalar();

        // if the power scalar for the bottom is zero, the correct thing to do is to zero out this cell
        // (justified using an inverse multiply argument in the working notes)
        if (powerScalar <= 0.0) continue;

        QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
        QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
        QUATERNION middle = (1.0 / powerScalar) * power * base.log().inverse();

        QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;

        QUATERNION Rderivative = topEval / bottomDerivative;
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        Real weight = _weightField(cellPosition);
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
        Real dSdt = dSdF * dFdt;

        if (isnan(dSdt))
        {
          _nans(x,y,z) = 1.0;
          continue;
        }
        if (isinf(dSdt))
        {
          _infs(x,y,z) = 1.0;
          continue;
        }

        final += dSdt;
      }

  int problems = _overflowed.sum() + _nans.sum() + _infs.sum();
  if (problems == 0) return final;
  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " overflows: " << _overflowed.sum() << endl;
  cout << " nans: " << _nans.sum() << endl;
  cout << " infs: " << _infs.sum() << endl;
  */
  _adaptiveFired = true;
  
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        if (_overflowed(x,y,z) + _nans(x,y,z) + _infs(x,y,z) < 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        long double xReal = point[0];
        long double yReal = point[1];
        long double zReal = point[2];
#else
        long double xReal = (long double)x / xRes * _lengths[0] + xMin + _center[0];
        long double yReal = (long double)y / yRes * _lengths[1] + yMin + _center[1];
        long double zReal = (long double)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNIONE iterate(xReal, yReal, zReal, _quaternionSlice);
     
        long double magnitude = iterate.magnitude();
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        vector<QUATERNIONE> forwardCacheTop = _topCacheExtended.forward(x,y,z);
        vector<QUATERNIONE> backwardCacheTop = _topCacheExtended.backward(x,y,z);
        vector<QUATERNIONE> forwardCacheBottom = _bottomCacheExtended.forward(x,y,z);
        vector<QUATERNIONE> backwardCacheBottom = _bottomCacheExtended.backward(x,y,z);

        // check that the cache was built
        if (_nans(x,y,z) > 0.5 || _infs(x,y,z) > 0.5)
        {
          //if (_topCacheExtended.forward(x,y,z).size() == 0)
          if (forwardCacheTop.size() == 0)
            QUATERNIONE topEval = _topExtended.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
          
          //if (_bottomCacheExtended.forward(x,y,z).size() == 0)
          if (forwardCacheBottom.size() == 0)
            QUATERNIONE bottomEval = _bottomExtended.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
        }

        const vector<QUATERNIONE>& topForward = forwardCacheTop;
        const vector<QUATERNIONE>& bottomForward = forwardCacheBottom;
        QUATERNIONE topEval = topForward.back();
        QUATERNIONE bottomEval = bottomForward.back();
        QUATERNIONE R = (topEval / bottomEval);

        long double Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        long double F = log(Rmagnitude) + log(_expScaling);

        QUATERNIONE left(1,0,0,0);
        QUATERNIONE right(1,0,0,0);
      
        if (_bottom.totalRoots() > 1)
        {
          // VERIFIED
          if (whichRoot > 0)
            left = bottomForward[whichRoot - 1];
          
          if (whichRoot < _bottom.totalRoots() - 1)
            right = backwardCacheBottom[whichRoot + 1];
        }

        // get the derivative of the root in question
        const long double powerScalar = _bottomExtended.powerScalar();
        QUATERNIONE base = (iterate - _bottomExtended.roots()[whichRoot]);
        QUATERNIONE power = base.pow(powerScalar * _bottomExtended.powers()[whichRoot]);
        QUATERNIONE middle = (1.0 / powerScalar) * power * base.log().inverse();

        QUATERNIONE bottomDerivative = QUATERNIONE(-1,0,0,0) * left * middle * right;

        QUATERNIONE Rderivative = topEval / bottomDerivative;
        long double dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
        long double weight = _weightField(cellPosition);
  
        long double tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
        long double dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
        long double dSdt = dSdF * dFdt;

        final += dSdt;
      }
  
  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomScaledPowerGradient(const int whichRoot) const
{
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  Real final = 0;
  for (int y = 0; y < yRes; y++)
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        Real magnitude = iterate.magnitude();
        //int totalIterations = 0;
        
        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

#if 0
        QUATERNION topEval = iterate * _top.evaluatePowerFactored(iterate);
        QUATERNION bottomEval = iterate * _bottom.evaluatePowerFactored(iterate);
#else        
        const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
        const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
// NEW: IMPROVES CONVERGENCE?
#if EXPLICIT_LEADING_ROOT
        QUATERNION topEval = topForward.back();
        QUATERNION bottomEval = bottomForward.back();
#else
        QUATERNION topEval = iterate * topForward.back();
        QUATERNION bottomEval = iterate * bottomForward.back();
#endif
#endif
        QUATERNION R = (topEval / bottomEval);
        //Real F = log(R.magnitude()) + log(_expScaling);

        // should probably track how often this happens, because it suggests we've run out of precision
        Real Rmagnitude = R.magnitude();

        // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
        //if (Rmagnitude <= 0.0) continue;
        if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

        Real F = log(Rmagnitude) + log(_expScaling);

#if 0
        QUATERNION bottomDerivative = _bottom.inversePowerDerivative(iterate, whichRoot);
#else
        QUATERNION left(1,0,0,0);
        QUATERNION right(1,0,0,0);
      
        if (_bottom.totalRoots() > 1)
        {
          // VERIFIED
          if (whichRoot > 0)
            left = bottomForward[whichRoot - 1];
          
          if (whichRoot < _bottom.totalRoots() - 1)
            right = _bottomCache.backward(x,y,z)[whichRoot + 1];
          /*
          // TODO: this can be optimized using the cache, but need more roots
          // to properly debug
          for (int i = 0; i < whichRoot; i++)
            left *= (iterate - _bottom.roots()[i]).pow(_bottom.powers()[i]);
          
          for (int i = whichRoot + 1; i < _bottom.totalRoots(); i++)
            right *= (iterate - _bottom.roots()[i]).pow(_bottom.powers()[i]);
            */
        }

        // get the derivative of the root in question
        const Real powerScalar = _bottom.powerScalar();
        QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
        QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
        QUATERNION middle = (1.0 / powerScalar) * power * base.log().inverse();

        QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;
#endif

// NEW: IMPROVES CONVERGENCE?
#if EXPLICIT_LEADING_ROOT
        QUATERNION Rderivative = topEval / bottomDerivative;
#else
        QUATERNION Rderivative = topEval / (iterate * bottomDerivative);
#endif
        Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
       
        VEC3F cellPosition = _fractal.cellCenter(x,y,z);
#if USING_WEIGHT
        Real weight = _weightField(cellPosition);
#else
        Real inverseDistance = _inverseDistanceField(cellPosition);
#endif
  
        //Real threshold = 0;
        //Real alpha = 10;
        Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
#if USING_WEIGHT
        Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
#else
        Real dSdF = -_tanhAlpha * inverseDistance * (1 - tanhTerm * tanhTerm); 
#endif
        Real dSdt = dSdF * dFdt;
        if (isnan(dSdt)) continue;

        final += dSdt;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomBulkScalarGradientAP()
{
  TIMER functionTimer(__FUNCTION__);
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  vector<Real> sums(zRes);

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
  {
    Real& threadSum = sums[z];
    threadSum = 0;

    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        Real magnitude = iterate.magnitude();

        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        for (int i = 0; i < _bottom.totalRoots(); i++)
        {
          // trap the root power = 0 exception -- in this case, the correct
          // answer is to produce a zero, but it will produce a NaN in this code.i
          if (_bottom.powers()[i] < 1e-8)
            continue;

          int whichRoot = i;
          const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
          const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
          QUATERNION topEval = topForward.back();
          QUATERNION bottomEval = bottomForward.back();
          QUATERNION R = (topEval / bottomEval);

          // should probably track how often this happens, because it suggests we've run out of precision
          Real Rmagnitude = R.magnitude();

          // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
          if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

          Real F = log(Rmagnitude) + log(_expScaling);

          QUATERNION left(1,0,0,0);
          QUATERNION right(1,0,0,0);
        
          if (_bottom.totalRoots() > 1)
          {
            if (whichRoot > 0)
              left = bottomForward[whichRoot - 1];
            
            if (whichRoot < _bottom.totalRoots() - 1)
              right = _bottomCache.backward(x,y,z)[whichRoot + 1];
          }

          // get the derivative of the root in question
          const Real powerScalar = _bottom.powerScalar();
          QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
          QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
          QUATERNION middle = (1.0 / _bottom.powers()[whichRoot]) * power * base.log().inverse();

          QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;

          QUATERNION Rderivative = topEval / (bottomDerivative);
          Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);
          Real weight = _weightField(cellPosition);
    
          //Real threshold = 0;
          //Real alpha = 10;
          Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          Real dSdt = dSdF * dFdt;

          if (isnan(dSdt))
          {
            _nans(x,y,z) = 1.0;
            continue;
          }
          if (isinf(dSdt))
          {
            _infs(x,y,z) = 1.0;
            continue;
          }

          threadSum += dSdt;
        }
      }
  }

  // if there were no issues, all done
  int problems = _overflowed.sum() + _nans.sum() + _infs.sum();
  if (problems == 0)
  { 
    Real final = 0;
    for (int z = 0; z < zRes; z++)
      final += sums[z];

    return final;
  }
  _adaptiveFired = true;

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
  {
    Real& threadSum = sums[z];
    threadSum = 0;

    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        if (_overflowed(x,y,z) + _nans(x,y,z) + _infs(x,y,z) < 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        long double xReal = point[0];
        long double yReal = point[1];
        long double zReal = point[2];
#else
        long double xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        long double yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        long double zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNIONE iterate(xReal, yReal, zReal, _quaternionSlice);
     
        long double magnitude = iterate.magnitude();

        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        // check that the cache was built
        if (_nans(x,y,z) > 0.5 || _infs(x,y,z) > 0.5)
        {
          if (_topCacheExtended.forward(x,y,z).size() == 0)
          {
            vector<QUATERNIONE>& forwardCacheTop = _topCacheExtended.forward(x,y,z);
            vector<QUATERNIONE>& backwardCacheTop = _topCacheExtended.backward(x,y,z);
            QUATERNIONE topEval = _topExtended.evaluateScaledPowerFactored(iterate, forwardCacheTop, backwardCacheTop);
          }
          if (_bottomCacheExtended.forward(x,y,z).size() == 0)
          {
            vector<QUATERNIONE>& forwardCacheBottom = _bottomCacheExtended.forward(x,y,z);
            vector<QUATERNIONE>& backwardCacheBottom = _bottomCacheExtended.backward(x,y,z);
            QUATERNIONE bottomEval = _bottomExtended.evaluateScaledPowerFactored(iterate, forwardCacheBottom, backwardCacheBottom);
          }
        }

        for (int i = 0; i < _bottomExtended.totalRoots(); i++)
        {
          // trap the root power = 0 exception -- in this case, the correct
          // answer is to produce a zero, but it will produce a NaN in this code.i
          if (_bottom.powers()[i] < 1e-8)
            continue;

          int whichRoot = i;
          const vector<QUATERNIONE>& topForward = _topCacheExtended.forward(x,y,z);
          const vector<QUATERNIONE>& bottomForward = _bottomCacheExtended.forward(x,y,z);
          QUATERNIONE topEval = topForward.back();
          QUATERNIONE bottomEval = bottomForward.back();
          QUATERNIONE R = (topEval / bottomEval);

          // should probably track how often this happens, because it suggests we've run out of precision
          long double Rmagnitude = R.magnitude();

          // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
          if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

          long double F = log(Rmagnitude) + log(_expScaling);

          QUATERNIONE left(1,0,0,0);
          QUATERNIONE right(1,0,0,0);
        
          if (_bottom.totalRoots() > 1)
          {
            if (whichRoot > 0)
              left = bottomForward[whichRoot - 1];
            
            if (whichRoot < _bottom.totalRoots() - 1)
              right = _bottomCacheExtended.backward(x,y,z)[whichRoot + 1];
          }

          // get the derivative of the root in question
          const long double powerScalar = _bottomExtended.powerScalar();
          QUATERNIONE base = (iterate - _bottomExtended.roots()[whichRoot]);
          QUATERNIONE power = base.pow(powerScalar * _bottomExtended.powers()[whichRoot]);
          QUATERNIONE middle = (1.0 / _bottomExtended.powers()[whichRoot]) * power * base.log().inverse();

          QUATERNIONE bottomDerivative = QUATERNIONE(-1,0,0,0) * left * middle * right;

          QUATERNIONE Rderivative = topEval / (bottomDerivative);
          long double dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);
          long double weight = _weightField(cellPosition);
    
          long double tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          long double dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          long double dSdt = dSdF * dFdt;

          threadSum += dSdt;
        }
      }
  }

  Real final = 0;
  for (int z = 0; z < zRes; z++)
    final += sums[z];

  return final;
}

///////////////////////////////////////////////////////////////////////
// compute the derivative of a root power
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomScaledPowerGradient() const
{
  TIMER functionTimer(__FUNCTION__);
  // compute it here, make sure it matches
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();

#if !RECENTER_ITERATE
  Real xMin = -_lengths[0] * 0.5;
  Real yMin = -_lengths[1] * 0.5; 
  Real zMin = -_lengths[2] * 0.5;
#endif

  vector<Real> sums(zRes);

#if OMP_ENABLED
#pragma omp parallel
#pragma omp for schedule(dynamic)
#endif
  for (int z = 0; z < zRes; z++)
  {
    Real& threadSum = sums[z];
    threadSum = 0;

    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // ADDING OVERFLOW DETECTION
        if (_overflowed(x,y,z) > 0.5) continue;

#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
     
        //Real distanceTravelled = 0.0;
        //Real magnitude = iterate.magnitude();
        //int totalIterations = 0;
        Real magnitude = iterate.magnitude();

        // if we hit the origin, move on
        if (magnitude < 1e-8) continue;

        for (int i = 0; i < _bottom.totalRoots(); i++)
        {
          // trap the root power = 0 exception -- in this case, the correct
          // answer is to produce a zero, but it will produce a NaN in this code.i
          if (_bottom.powers()[i] < 1e-8)
            continue;

          int whichRoot = i;
          const vector<QUATERNION>& topForward = _topCache.forward(x,y,z);
          const vector<QUATERNION>& bottomForward = _bottomCache.forward(x,y,z);
#if EXPLICIT_LEADING_ROOT
          QUATERNION topEval = topForward.back();
          QUATERNION bottomEval = bottomForward.back();
#else
          QUATERNION topEval = iterate * topForward.back();
          QUATERNION bottomEval = iterate * bottomForward.back();
#endif
          QUATERNION R = (topEval / bottomEval);

          // should probably track how often this happens, because it suggests we've run out of precision
          Real Rmagnitude = R.magnitude();

          // if R is zero, just bail. Almost certainly a precision problem, or we've stumbled on the origin.
          //if (Rmagnitude <= 0.0) continue;
          if (Rmagnitude <= 0.0 || isinf(Rmagnitude)) continue;

          Real F = log(Rmagnitude) + log(_expScaling);
          //Real F = log(R.magnitude()) + log(_expScaling);

          QUATERNION left(1,0,0,0);
          QUATERNION right(1,0,0,0);
        
          if (_bottom.totalRoots() > 1)
          {
            if (whichRoot > 0)
              left = bottomForward[whichRoot - 1];
            
            if (whichRoot < _bottom.totalRoots() - 1)
              right = _bottomCache.backward(x,y,z)[whichRoot + 1];
          }

#if 1
          // get the derivative of the root in question
          const Real powerScalar = _bottom.powerScalar();
          QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
          QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
          QUATERNION middle = (1.0 / _bottom.powers()[whichRoot]) * power * base.log().inverse();

          QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;

#if EXPLICIT_LEADING_ROOT
          QUATERNION Rderivative = topEval / (bottomDerivative);
#else
          QUATERNION Rderivative = topEval / (iterate * bottomDerivative);
#endif
          Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);
          Real weight = _weightField(cellPosition);
    
          //Real threshold = 0;
          //Real alpha = 10;
          Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          Real dSdt = dSdF * dFdt;

          if (isnan(dSdt)) continue;

          threadSum += dSdt;
#else
          // get the derivative of the root in question
          const Real powerScalar = _bottom.powerScalar();
          QUATERNION base = (iterate - _bottom.roots()[whichRoot]);
          QUATERNION power = base.pow(powerScalar * _bottom.powers()[whichRoot]);
          //QUATERNION middle = _bottom.powers()[whichRoot] * power * base.log().inverse();
          //Real inversePower = (_bottom.powers()[whichRoot] > 1e-8) ? 1.0 / _bottom.powers()[whichRoot] : 0;
          //QUATERNION middle = inversePower * power * base.log().inverse();
          QUATERNION middle = (1.0 / _bottom.powers()[whichRoot]) * power * base.log().inverse();

          QUATERNION bottomDerivative = QUATERNION(-1,0,0,0) * left * middle * right;

#if EXPLICIT_LEADING_ROOT
          QUATERNION Rderivative = topEval / bottomDerivative;
#else
          QUATERNION Rderivative = topEval / (iterate * bottomDerivative);
#endif
          Real dFdt = (1.0 / R.dot(R)) * (R.dot(Rderivative));
         
          VEC3F cellPosition = _fractal.cellCenter(x,y,z);
          Real weight = _weightField(cellPosition);
    
          //Real threshold = 0;
          //Real alpha = 10;
          Real tanhTerm = tanh(_tanhAlpha * (F - _tanhThreshold));
          Real dSdF = -_tanhAlpha * weight * (1 - tanhTerm * tanhTerm);
          Real dSdt = dSdF * dFdt;
          if (isnan(dSdt)) continue;

          threadSum += dSdt;
          if (isnan(final))
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << " NaN found! " << endl;
            cout << " dSdt: " << dSdt << endl;
            cout << " dSdF: " << dSdF << endl;
            cout << " dFdt: " << dFdt << endl;
            cout << " R dot R: " << R.dot(R) << endl;
            cout << " R dot dR: " << R.dot(Rderivative) << endl;
            cout << " R: " << R << endl;
            cout << " topEval: " << topEval << endl;
            cout << " iterate: " << iterate << endl;
            cout << " bottomDerivative: " << bottomDerivative << endl;
            cout << " left: " << left << endl;
            cout << " middle: " << middle << endl;
            cout << " right: " << right << endl;
            cout << " power: " << power << endl;
            cout << " log inverse: " << base.log().inverse() << endl;
            cout << " bottom power:" << _bottom.powers()[whichRoot] << endl;
            exit(0);
          }
#endif
        }
      }
// inconsistent sum ordering here creates problems
//#pragma omp atomic
//    final += threadSum;
  }

  Real final = 0;
  for (int z = 0; z < zRes; z++)
    final += sums[z];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeTopBulkPowerGradient() const
{
  VECTOR gradients(_top.totalRoots());
#pragma omp parallel
#pragma omp for schedule(static)
  for (int x = 0; x < _top.totalRoots(); x++)
    gradients[x] = computeTopPowerGradient(x);
  return gradients.sum();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real OPTIMIZE_3D::computeBottomBulkPowerGradient() const
{
  VECTOR gradients(_bottom.totalRoots());
#pragma omp parallel
#pragma omp for schedule(static)
  for (int x = 0; x < _bottom.totalRoots(); x++)
    gradients[x] = computeBottomPowerGradient(x);
  return gradients.sum();
}

///////////////////////////////////////////////////////////////////////
// print out the current powers of the roots
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::printRootPowers() const
{
  //VECTOR::printVertical = false;
  //cout << " Top powers: " << VECTOR(_top.powers()) << endl;
  //cout << " Bottom powers: " << VECTOR(_bottom.powers()) << endl;
}

///////////////////////////////////////////////////////////////////////
// print out all the roots
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::printRoots() const
{
  cout << " Top polynomial roots " << endl;
  for (int x = 0; x < _top.totalRoots(); x++)
    cout << "\tPosition: " << _top.roots()[x] << "\tPower: " << _top.powers()[x] << endl;
  cout << " Bottom polynomial roots " << endl;
  for (int x = 0; x < _bottom.totalRoots(); x++)
    cout << "\tPosition: " << _bottom.roots()[x] << "\tPower: " << _bottom.powers()[x] << endl;
}

///////////////////////////////////////////////////////////////////////
// compute some arbitrary weight field
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeWeightField()
{
  _weightField = _inverseDistanceField;
}

///////////////////////////////////////////////////////////////////////
// reweight the score weight field with curvature
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::reweightWithCurvature(const FIELD_3D& curvatureField)
{
  assert(curvatureField.xRes() == _weightField.xRes());
  assert(curvatureField.yRes() == _weightField.yRes());
  assert(curvatureField.zRes() == _weightField.zRes());

  //FIELDVIEW3D(_weightField);
  //FIELDVIEW3D(curvatureField);
  //exit(0);
  for (int x = 0; x < _weightField.totalCells(); x++)
    _weightField[x] *= curvatureField[x];

  // normalize the whole field
  FIELD_3D temp = _weightField;
  temp.absoluteValue();
  _weightField *= 1.0 / temp.sum();

  Real insideSum = 0;
  Real outsideSum = 0;
  for (int x = 0; x < _weightField.totalCells(); x++)
  {
    if (_inverseDistanceField[x] < 0)
      insideSum += _weightField[x];
    else
      outsideSum += _weightField[x];
  }

  cout << " inside sum: " << insideSum << endl;
  cout << " outside sum: " << outsideSum << endl;

  Real insideNormalize = 0.5 / fabs(insideSum);
  Real outsideNormalize = 0.5 / fabs(outsideSum);

  for (int x = 0; x < _weightField.totalCells(); x++)
  {
    if (_inverseDistanceField[x] < 0)
      _weightField[x] *= insideNormalize;
    else
      _weightField[x] *= outsideNormalize;
  }
  
  //FIELDVIEW3D(_weightField);
  //exit(0);
}

///////////////////////////////////////////////////////////////////////
// refresh polynomial evaluation caches
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::refreshPolynomialCache(POLYNOMIAL_CACHE& cache, const POLYNOMIAL_4D& polynomial, const FIELD_3D& fractal)
{
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();
  
#if !RECENTER_ITERATE
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;
#endif
  
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
        
        vector<QUATERNION>& forwardCache = cache.forward(x,y,z);
        vector<QUATERNION>& backwardCache = cache.backward(x,y,z);
        QUATERNION eval = polynomial.evaluatePowerFactored(iterate, forwardCache, backwardCache);
      }
  }
}

///////////////////////////////////////////////////////////////////////
// compute cached quaternion logs for everybody
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeLogCache(QUATERNION_CACHE& cache, const POLYNOMIAL_4D& polynomial, const FIELD_3D& fractal, const bool inverse)
{
  const int xRes = fractal.xRes();
  const int yRes = fractal.yRes();
  const int zRes = fractal.zRes();
  
#if !RECENTER_ITERATE
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;
#endif
  
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
      
        // stomp any previous results  
        vector<QUATERNION>& logCache = cache(x,y,z);
        logCache.clear();

        // make some new ones
        for (int i = 0; i < polynomial.totalRoots(); i++)
        {
          QUATERNION base = (iterate - polynomial.roots()[i]);
          QUATERNION baseLog = base.log();
          if (inverse)
            baseLog = baseLog.inverse();
          logCache.push_back(baseLog);
        }
      }
  }
}

///////////////////////////////////////////////////////////////////////
// compute some cached quantities for use later in the Hessian
//
// Assumes that the _topCache and _bottomCache are up to date!!!!
///////////////////////////////////////////////////////////////////////
void OPTIMIZE_3D::computeHessianCaches()
{
  const int xRes = _fractal.xRes();
  const int yRes = _fractal.yRes();
  const int zRes = _fractal.zRes();
  const int slabSize = xRes * yRes;
  
#if !RECENTER_ITERATE
  const Real xMin = -_lengths[0] * 0.5;
  const Real yMin = -_lengths[1] * 0.5; 
  const Real zMin = -_lengths[2] * 0.5;
#endif

  // make sure everything is the right size 
  if (_RdotR.xRes() != xRes || _RdotR.yRes() != yRes || _RdotR.zRes() != zRes)
  {
    _RdotR = _fractal;
    _RdotR.clear();

    _conformalConformal = _fractal;
    _conformalConformal.clear();
    
    _dSdF = _fractal;
    _dSdF.clear();
    
    _topDerivativeCache = QUATERNION_CACHE(_fractal.xRes(), _fractal.yRes(), _fractal.zRes());
    _bottomDerivativeCache = QUATERNION_CACHE(_fractal.xRes(), _fractal.yRes(), _fractal.zRes());
  }

  const Real tanAlphaSq = _tanhAlpha * _tanhAlpha;
  for (int y = 0; y < yRes; y++)
  {
    for (int z = 0; z < zRes; z++)
      for (int x = 0; x < xRes; x++)
      {
#if RECENTER_ITERATE
        VEC3F point = _fractal.cellCenter(x,y,z);
        Real xReal = point[0];
        Real yReal = point[1];
        Real zReal = point[2];
#else
        Real xReal = (Real)x / xRes * _lengths[0] + xMin + _center[0];
        Real yReal = (Real)y / yRes * _lengths[1] + yMin + _center[1];
        Real zReal = (Real)z / zRes * _lengths[2] + zMin + _center[2];
#endif
        int index = x + y * xRes + z * slabSize;

        QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
      
        Real weight = _weightField(_fractal.cellCenter(index));
        const vector<QUATERNION>& forwardCacheTop = _topCache.forward(x,y,z);
        const vector<QUATERNION>& forwardCacheBottom = _bottomCache.forward(x,y,z);

        const QUATERNION topEval    = forwardCacheTop.back();
        const QUATERNION bottomEval = forwardCacheBottom.back();
        const QUATERNION bottomEvalInverse = bottomEval.inverse();
        const QUATERNION R = topEval * bottomEvalInverse;

        _RdotR[index] = R.dot(R);

        const Real F = log(R.magnitude()) + log(_expScaling);
        const Real interior = _tanhAlpha * (F - _tanhThreshold);
        const Real tanhTerm = tanh(interior);

        const Real sechTermSq = 1.0 - tanhTerm * tanhTerm;
        const Real dSdF = -_tanhAlpha * weight * sechTermSq;
        _dSdF[index] = dSdF;
        _conformalConformal[index] = 2.0 * tanAlphaSq * weight * sechTermSq * tanhTerm;

        const vector<QUATERNION>& backwardCacheTop = _topCache.backward(x,y,z);
        const Real powerScalarTop = _top.powerScalar();
        vector<QUATERNION>& topDerivatives = _topDerivativeCache(x, y, z);
        for (int r = 0; r < _top.totalRoots(); r++)
        {
          const QUATERNION left  = (r != 0) ? forwardCacheTop[r] : QUATERNION(1,0,0,0);
          const QUATERNION right = (r != _top.totalRoots() - 1) ? backwardCacheTop[r + 1] : QUATERNION(1,0,0,0);
          const QUATERNION base    = (iterate - _top.roots()[r]);
          const QUATERNION baseLog = base.log();
          const QUATERNION middle = powerScalarTop * baseLog;
          topDerivatives.push_back(left * middle * right * bottomEvalInverse);
        }
      }
  }
}

///////////////////////////////////////////////////////////////////////
// build the whole Hessian
///////////////////////////////////////////////////////////////////////
MATRIX OPTIMIZE_3D::buildHessian()
{
  TIMER functionTimer(__FUNCTION__);
  int size = 1 + _top.totalRoots() + _bottom.totalRoots();
  int bottomOffset = 1;
  int topOffset = 1 + _bottom.totalRoots();
  MATRIX final(size, size);

  // populate the diagonal first
  final(0,0) = computeConformalConformalHessian();
  TIMER bottomDiagonalTimer("Bottom Diagonal");
  cout << " Computing bottom diagonal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _bottom.totalRoots(); x++)
    final(x + bottomOffset, x + bottomOffset) = computeBottomDiagonalHessian(x);
  cout << " done." << endl;
  bottomDiagonalTimer.stop();

  TIMER topDiagonalTimer("Top Diagonal");
  cout << " Computing top diagonal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _top.totalRoots(); x++)
    final(x + topOffset, x + topOffset) = computeTopDiagonalHessian(x);
  cout << " done." << endl;
  topDiagonalTimer.stop();

  // populate all the conformal-bottom entries
  TIMER bottomConformalTimer("Bottom Conformal");
  cout << " Computing bottom-conformal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _bottom.totalRoots(); x++)
    final(x + bottomOffset, 0) = computeBottomConformalHessian(x);
  cout << " done." << endl;
  bottomConformalTimer.stop();

  // populate all the conformal-top entries
  TIMER topConformalTimer("Top Conformal");
  cout << " Computing top-conformal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _top.totalRoots(); x++)
    final(x + topOffset, 0) = computeTopConformalHessian(x);
  cout << " done." << endl;
  topConformalTimer.stop();

  // populate all the top-top entries
  TIMER topTopTimer("Mixed Top");
  cout << " Computing mixed-top ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _top.totalRoots(); x++)
    for (int y = 0; y < x; y++) // skip (0,0) already got it on the diagonal
    {
      Real mixed = computeMixedTopHessian(x,y);
      final(x + topOffset, y + topOffset) = mixed;
    }
  cout << " done." << endl;
  topTopTimer.stop();
  
  // populate all the bottom-bottom entries
  TIMER bottomBottomTimer("Mixed Bottom");
  cout << " Computing mixed-bottom..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _bottom.totalRoots(); x++)
    for (int y = 0; y < x; y++) // skip (0,0) already got it on the diagonal
    {
      Real mixed = computeMixedBottomHessian(x,y);
      final(x + bottomOffset, y + bottomOffset) = mixed;
    }
  cout << " done." << endl;
  bottomBottomTimer.stop();
  
  // populate all the top-bottom entries
  TIMER topBottomTimer("Top Bottom");
  cout << " Computing top-bottom..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int t = 0; t < _top.totalRoots(); t++)
    for (int b = 0; b < _bottom.totalRoots(); b++)
    {
      Real mixed = computeTopBottomHessian(t,b);
      final(t + topOffset, b + bottomOffset) = mixed;
    }
  cout << " done." << endl;
  topBottomTimer.stop();

  //final.writeMatlab("beforeSym.m", "halfHessian");

  // do the symmetric copy
  for (int x = 0; x < size; x++)
    for (int y = 0; y < x; y++)
      final(y,x) = final(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// build the whole Hessian
///////////////////////////////////////////////////////////////////////
MATRIX OPTIMIZE_3D::buildHessianFast()
{
  TIMER functionTimer(__FUNCTION__);
  int size = 1 + _top.totalRoots() + _bottom.totalRoots();
  int bottomOffset = 1;
  int topOffset = 1 + _bottom.totalRoots();
  MATRIX final(size, size);

  if (_topLogCache.totalCells() != _fractal.totalCells())
    _topLogCache = QUATERNION_CACHE(_fractal.xRes(), _fractal.yRes(), _fractal.zRes());
  if (_bottomLogCache.totalCells() != _fractal.totalCells())
    _bottomLogCache = QUATERNION_CACHE(_fractal.xRes(), _fractal.yRes(), _fractal.zRes());

  cout << " Computing log cache ..." << flush;
  computeLogCache(_topLogCache, _top, _fractal);
  computeLogCache(_bottomLogCache, _bottom, _fractal, true);

  cout << " other caches ..." << flush;
  computeHessianCaches();

  // populate the diagonal first
  final(0,0) = computeConformalConformalHessian();
  TIMER bottomDiagonalTimer("Bottom Diagonal Fast");
  cout << " bottom diagonal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _bottom.totalRoots(); x++)
    final(x + bottomOffset, x + bottomOffset) = computeBottomDiagonalHessianFast(x);
    //final(x + bottomOffset, x + bottomOffset) = computeBottomDiagonalHessian(x);
  bottomDiagonalTimer.stop();

  TIMER topDiagonalTimer("Top Diagonal Fast");
  cout << " top diagonal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _top.totalRoots(); x++)
    final(x + topOffset, x + topOffset) = computeTopDiagonalHessianFast(x);
    //final(x + topOffset, x + topOffset) = computeTopDiagonalHessian(x);
  topDiagonalTimer.stop();

  // populate all the conformal-bottom entries
  TIMER bottomConformalTimer("Bottom Conformal Fast");
  cout << " bottom-conformal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _bottom.totalRoots(); x++)
    final(x + bottomOffset, 0) = computeBottomConformalHessianFast(x);
    //final(x + bottomOffset, 0) = computeBottomConformalHessian(x);
  bottomConformalTimer.stop();

  // populate all the conformal-top entries
  TIMER topConformalTimer("Top Conformal Fast");
  cout << " top-conformal ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _top.totalRoots(); x++)
    final(x + topOffset, 0) = computeTopConformalHessianFast(x);
    //final(x + topOffset, 0) = computeTopConformalHessian(x);
  topConformalTimer.stop();

  // populate all the top-top entries
  TIMER topTopTimer("Mixed Top Fast");
  cout << " mixed-top ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _top.totalRoots(); x++)
  {
    for (int y = 0; y < x; y++)
    {
      Real mixed = computeMixedTopHessianFast(x,y);
      final(x + topOffset, y + topOffset) = mixed;
    }

    //if (x % 10 == 0 && x != 0)
    //  TIMER::printTimings();
  }
  topTopTimer.stop();
  
  // populate all the bottom-bottom entries
  TIMER bottomBottomTimer("Mixed Bottom Fast");
  cout << " mixed-bottom ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int x = 0; x < _bottom.totalRoots(); x++)
    for (int y = 0; y < x; y++)
    {
      Real mixed = computeMixedBottomHessianFast(x,y);
      final(x + bottomOffset, y + bottomOffset) = mixed;
    }
  bottomBottomTimer.stop();
  
  // populate all the top-bottom entries
  TIMER topBottomTimer("Top Bottom Fast");
  cout << " top-bottom ..." << flush;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int t = 0; t < _top.totalRoots(); t++)
    for (int b = 0; b < _bottom.totalRoots(); b++)
    {
      Real mixed = computeTopBottomHessianFast(t,b);
      final(t + topOffset, b + bottomOffset) = mixed;
    }
  cout << " done." << endl;
  topBottomTimer.stop();

  //final.writeMatlab("beforeSym.m", "halfHessian");

  // do the symmetric copy
  for (int x = 0; x < size; x++)
    for (int y = 0; y < x; y++)
      final(y,x) = final(x,y);

  return final;
}
  
