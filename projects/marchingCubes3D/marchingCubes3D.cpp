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
#include <float.h>

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>

#include <iostream>
#include "QUATERNION.h"
#include "POLYNOMIAL_4D.h"
#include "OPTIMIZE_3D.h"
#include "TRIANGLE_MESH.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Meshes and integrators
//////////////////////////////////////////////////////////////////////////////
TRIANGLE_MESH* triangleMesh = NULL;
TRIANGLE_MESH* distanceMesh = NULL;
TRIANGLE_MESH* fractalMesh = NULL;
FIELD_3D toMarch(100, 100, 100);
Real marchingSurface = 0;

OPTIMIZE_3D optimize3D;
string outputFilename("./temp/temp.obj");

bool topPolynomial = true;
bool drawDistanceMesh = false;

int maxIterations = 1;
int res = 97;

///////////////////////////////////////////////////////////////////////
// compute the fitting score for a given threshold
///////////////////////////////////////////////////////////////////////
Real distanceScore(const FIELD_3D& logFractal, const FIELD_3D& distanceField, const Real& threshold)
{
  Real score = 0;
  for (int x = 0; x < logFractal.totalCells(); x++)
  {
    //Real distanceThresholded = (distanceField[x] < 0) ? 1 : -1;
    Real distanceThresholded = (distanceField[x] < 0) ? -1 : 1;
    Real fractalThresholded = (logFractal[x] > threshold) ? -1 : 1;
    score += fractalThresholded  * distanceThresholded;
  }
  score *= 1.0 / logFractal.totalCells();
  return score;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void computeFractal(int res)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing fractal with res " << res << " ... " << flush;

  VEC3I reses(res,res,res);
  VEC3F center = optimize3D.fractal().center();
  VEC3F lengths = optimize3D.fractal().lengths();
  triangleMesh = new TRIANGLE_MESH(center, lengths, reses, optimize3D.top(), optimize3D.bottom(), optimize3D.expScaling(), optimize3D.maxIterations(), optimize3D.slice(), 0);
  cout << " Marching cubes found " << triangleMesh->triangles().size() << " triangles and " << triangleMesh->vertices().size() << " vertices " << endl;

  VEC3F minBox, maxBox;
  triangleMesh->boundingBox(minBox, maxBox);
  cout << " final fractal bounding box: " << minBox << " " << maxBox << endl;

  TIMER::printTimings();
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (sizeof(Real) == sizeof(double))
    cout << " OPTIMIZING USING DOUBLE PRECISION " << endl;
  else if (sizeof(Real) == sizeof(float))
    cout << " OPTIMIZING USING SINGLE PRECISION " << endl;
  else if (sizeof(Real) == sizeof(long double))
    cout << " OPTIMIZING USING EXTENDED PRECISION " << endl;
  else
  {
    cout << " OPTIMIZING USING AN UNDOCUMENTED PRECISION, size: " << sizeof(Real) << endl;
  }

  TIMER functionTimer(__FUNCTION__);
  if (argc < 2)
  {
    cout << " USAGE: " << argv[0] << " <optimize 3D file> <res (optional)> <max iterations (optional)> <output filename (optional)> <headless?>" << endl;
    return 0;
  }

  optimize3D.read(argv[1]);

  if (argc >= 3)
    res = atoi(argv[2]);
  cout << " Using resolution: " << res << endl;

  if (argc >= 4)
    maxIterations = atoi(argv[3]);
  cout << " Using max iterations: " << maxIterations << endl;
  optimize3D.maxIterations() = maxIterations;
 
  if (argc >= 5)
    outputFilename = string(argv[4]);
  cout << " Using output filename: " << outputFilename.c_str() << endl;

  cout << " Conformal radius is: " << optimize3D.expScaling() << endl;
  optimize3D.printRoots();

  VEC3F center = optimize3D.fractal().center();
  VEC3F lengths = optimize3D.fractal().lengths();

  computeFractal(res);
  //Real globalMin = optimize3D.fractal().fieldMin();
  //Real globalMax = optimize3D.fractal().fieldMax();

  // print out the root powers, just to see if it all checks out
  triangleMesh->writeOBJ(outputFilename.c_str());
  //triangleMesh->writeOBJgz(outputFilename.c_str());

  return 0;
}
