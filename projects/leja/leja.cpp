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
#include "TRIANGLE_MESH.h"
#include <MERSENNETWISTER.h>
#include <queue>
#include <algorithm>
#include <float.h>
#include "OPTIMIZE_3D.h"
#include "POLYNOMIAL_4D.h"
#include "SPHERE.h"
#include "VECTOR3_FIELD_3D.h"

#include "SIMPLE_PARSER.h"

// used the squared distance? If so, this is not actually potential anymore ...
#define USING_SQUARED 0

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Meshes and integrators
//////////////////////////////////////////////////////////////////////////////
TRIANGLE_MESH* triangleMesh;
FIELD_3D distanceField;
FIELD_3D colorField;
FIELD_3D offsetColors;
FIELD_3D surfaceShell;
FIELD_3D shellPotential;
float marchingSurface = 0;

// the potential being approximated at an outer shell
vector<int> shellIndices;
vector<VEC3F> shellPositions;

// the surface points being used to approximate the outer shell
vector<int> surfaceIndices;
vector<VEC3F> surfacePositions;

vector<VEC3F> rootPositions;
vector<VEC3F> plyPoints;

//Real relativeErrorThreshold = 0.04;
//Real relativeErrorThreshold = 0.03;
Real relativeErrorThreshold = 0.02;

vector<SPHERE> topSpheres;
vector<SPHERE> bottomSpheres;
vector<Real> topWeights;
vector<Real> bottomWeights;

vector<VEC3F> topPoints;

Real bottomPower = 3.0;

////////////////////////////////////////////////////////////////////////////
// write out an oriented point cloud to PLY
////////////////////////////////////////////////////////////////////////////
bool writePlyPoints(const string& filename, const vector<VEC3F>& points)
{
  FILE* file = NULL;
  file = fopen(filename.c_str(), "w");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Couldn't open file " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  fprintf(file, "ply\n");
  //fprintf(file, "format ascii 1.0\n");
  fprintf(file, "format binary_little_endian 1.0\n");

  int totalVertices = points.size();
  fprintf(file, "element vertex %i\n", totalVertices);
  fprintf(file, "property float x\n");
  fprintf(file, "property float y\n");
  fprintf(file, "property float z\n");
  fprintf(file, "element face 0\n");
  fprintf(file, "property list uchar int vertex_index\n");
  fprintf(file, "end_header\n");

  // close the text file
  fclose(file);

  // reopen in binary for output
  file = fopen(filename.c_str(), "ab");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Couldn't open file " << filename.c_str() << " for binary output!!!" << endl;
    exit(0);
  }

  for (int x = 0; x < totalVertices; x++)
  {
    float vertex[] = {points[x][0], points[x][1], points[x][2]};
    fwrite((void*)&(vertex[0]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[1]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[2]), sizeof(float), 1, file);
  }

  // close the binary file
  fclose(file);

  return true;
}


//////////////////////////////////////////////////////////////////////////////
// Read in a PLY file of points
//////////////////////////////////////////////////////////////////////////////
void readPointsPLY(const string& filename, vector<VEC3F>& points)
{
  points.clear();

  FILE* file = NULL;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Couldn't open file " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  fscanf(file, "ply\n");
  fscanf(file, "format binary_little_endian 1.0\n");

  int totalPoints;
  fscanf(file, "element vertex %i\n", &totalPoints);
  fscanf(file, "property float x\n");
  fscanf(file, "property float y\n");
  fscanf(file, "property float z\n");
  fscanf(file, "element face 0\n");
  fscanf(file, "property list uchar int vertex_index\n");
  fscanf(file, "end_header\n");
  cout << " scanning for " << totalPoints << " points " << endl;

  Real maxMagnitude = 0;
  for (int x = 0; x < totalPoints; x++)
  {
    float vertex[3];
    fread((void*)&(vertex[0]), sizeof(float), 1, file);
    fread((void*)&(vertex[1]), sizeof(float), 1, file);
    fread((void*)&(vertex[2]), sizeof(float), 1, file);

    VEC3F vertexFinal(vertex[0], vertex[1], vertex[2]);

    if (vertexFinal.magnitude() > 1000)
      continue;

    if (vertexFinal.magnitude() > maxMagnitude)
    {
      cout << " max found: " << vertexFinal << endl;
      maxMagnitude = vertexFinal.magnitude();
    }
    
    points.push_back(vertexFinal);
  }

  cout << " max point magnitude found: " << maxMagnitude << endl;

  // close the binary file
  fclose(file);
}

//////////////////////////////////////////////////////////////////////////////
// shuffle a vector in a deterministic way using Mersenne Twister.
// C++11 has this built in, but too late for that now ...
//////////////////////////////////////////////////////////////////////////////
vector<VEC3F> shuffle(vector<VEC3F>& toShuffle)
{
  TIMER functionTimer(__FUNCTION__);
  MERSENNETWISTER twister(123456);
  vector<VEC3F> final(toShuffle.size());
  for (unsigned int x = 0; x < toShuffle.size(); x++)
  {
    int back = toShuffle.size() - 1;

    // pick a random element
    int pick = twister.randInt(back);

    final[x] = toShuffle[pick];
    final[pick] = toShuffle[x];
  }

  return final;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void buildLejaPoints()
{
  // Shuffle the random points
  cout << " Shuffling ... " << flush;
  plyPoints = shuffle(plyPoints);
  cout << endl;

  // pick one at random as the first Leja point
  MERSENNETWISTER twister(123456);

  map<int, bool> picked;
  int first = twister.randInt(plyPoints.size());
  topSpheres.push_back(SPHERE(plyPoints[first], 0.02));
  topPoints.push_back(plyPoints[first]);
  picked[first] = true; 

  const VEC3F& currentLeja = topSpheres.back().center();

  // initialize the Leja products
  vector<double> lejaProducts;
  for (int x = 0; x < plyPoints.size(); x++)
  {
    double distance = (plyPoints[x] - currentLeja).magnitude();
    lejaProducts.push_back(distance);
  }

  // add the next Leja Points
  const int totalLejaPoints = 100;
  for (int x = 0; x < totalLejaPoints; x++)
  {
    // find the smallest product
    int largestIndex = -1;
    double largestDistance = 0;
    for (int y = 0; y < lejaProducts.size(); y++)
    {
      if (picked.find(y) != picked.end())
        continue;

      if (lejaProducts[y] > largestDistance)
      {
        largestIndex = y;
        largestDistance = lejaProducts[y];
      }
    }
    assert(largestIndex != -1);

    // add it to the list
    topSpheres.push_back(SPHERE(plyPoints[largestIndex], 0.02));
    topPoints.push_back(plyPoints[largestIndex]);
    picked[largestIndex] = true;
    cout << " Add point " << x << "\t" << plyPoints[largestIndex] << "\t with distance " << largestDistance << endl;

    // update the rest of the Leja products
    const VEC3F& currentLeja = topSpheres.back().center();
    for (int y = 0; y < lejaProducts.size(); y++)
    {
      double distance = (plyPoints[y] - currentLeja).magnitude();
      lejaProducts[y] *= (distance);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void buildPoissonSourceSimulation()
{
  // do it in a way with a known implementation
  //srand(123456);
  //random_shuffle(plyPoints.begin(), plyPoints.end());
  cout << " Shuffling ... " << flush;
  plyPoints = shuffle(plyPoints);
  cout << endl;

  TIMER functionTimer(__FUNCTION__);
  // build the residual vector
  VECTOR targetPotential(shellIndices.size());

  // lift the original potential from the field
  for (unsigned int x = 0; x < shellIndices.size(); x++)
    targetPotential[x] = shellPotential[shellIndices[x]];

  VECTOR residual = targetPotential;
  vector<VECTOR> bestColumns;

  VECTOR solution;
  Real relative = 1;

  while (relative > relativeErrorThreshold)
  {
    // build the candidate set
    int candidateSetSize = 100;
    vector<VEC3F> candidateSet;
    for (int x = 0; x < candidateSetSize; x++)
    {
      candidateSet.push_back(plyPoints.back());
      plyPoints.pop_back();
    }

    cout << " Source " << bestColumns.size() << ", " << flush;

    // build the column for each candidate site
    cout << " Building columns ... " << flush;
    vector<VECTOR> candidateColumns(candidateSetSize);
#pragma omp parallel
#pragma omp for  schedule(dynamic)
    for (int x = 0; x < candidateSetSize; x++)
    {
      // compute the potential here
      VECTOR potential(targetPotential.size());
      
      // where is the center of the pole
      VEC3F pole = candidateSet[x];

      for (int y = 0; y < potential.size(); y++)
      {
        // get the radius
        Real radius = (shellPositions[y] - pole).magnitude();

        // assume a unit charge
#if USING_SQUARED
        potential[y] = 1.0 / (radius * radius);
#else
        potential[y] = 1.0 / radius;
#endif
      }
      candidateColumns[x] = potential;
    }
    cout << " done. " << flush;

    // find the site with the best projection onto the residual
    TIMER dotTimer("Dot Products");

    VECTOR dots(candidateSetSize);
#pragma omp parallel
#pragma omp for  schedule(dynamic)
    for (int x = 0; x < candidateSetSize; x++)
    {
      // get the projection onto the residual
      Real dot = fabs(candidateColumns[x].dot(residual));
      dots[x] = dot;
    }

    Real bestSeen = 0;
    int bestSeenIndex = -1;
    for (int x = 0; x < dots.size(); x++)
    {
      if (dots[x] > bestSeen)
      {
        bestSeen = dots[x];
        bestSeenIndex = x;
      }
    }

    dotTimer.stop();
    //cout << " Best candidate seen: " << bestSeen << "\t index: " << bestSeenIndex << endl;
    assert(bestSeenIndex >= 0);

    rootPositions.push_back(candidateSet[bestSeenIndex]);
    bestColumns.push_back(candidateColumns[bestSeenIndex]);

    const int rows = bestColumns[0].size();
    const int cols = bestColumns.size();
    MATRIX A(rows, cols);

    for (int y = 0; y < cols; y++)
      for (int x = 0; x < rows; x++)
        A(x,y) = bestColumns[y][x];

    MATRIX newA = A;

    // do it the Eigen way
    VECTOR b(targetPotential);
    //bool success = A.solveLeastSquares(b);
    //assert(success);
    //solution = VECTOR(bestColumns.size());
    //for (int x = 0; x < bestColumns.size(); x++)
    //  solution[x] = b[x];
    solution = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

    VECTOR fit = newA * solution;
    residual = targetPotential - fit;
    relative = residual.norm() / targetPotential.norm();
    cout << " Weight: " << solution[solution.size() - 1] << " \t Residual after fit: " << residual.norm() << " relative: " << relative << endl;

    cout << " New root magnitude: " << candidateSet[bestSeenIndex].magnitude() << endl;
  }

  FIELD_3D residualField(shellPotential);
  residualField = 0;
  for (unsigned int x = 0; x < shellIndices.size(); x++)
    residualField[shellIndices[x]] = residual[x];

  int positives = 0;
  Real topSum = 0;
  Real bottomSum = 0;
  topSpheres.clear();
  topWeights.clear();
  bottomSpheres.clear();
  bottomWeights.clear();
  for (int x = 0; x < solution.size(); x++)
  {
    if (solution[x] > 0)
    {
      positives++;
      topSpheres.push_back(SPHERE(rootPositions[x], 0.02));
      topWeights.push_back(solution[x]);
      topSum += solution[x];
    }
    else
    {
      bottomSpheres.push_back(SPHERE(rootPositions[x], 0.02));
      bottomWeights.push_back(solution[x]);
      bottomSum += solution[x];
    }
  }

  cout << " positive roots: " << topWeights.size() << endl;
  cout << " negative roots: " << bottomWeights.size() << endl;

  cout << " top sum:    " << topSum << endl;
  cout << " bottom sum: " << bottomSum << endl;

  TIMER::printTimings();
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  TIMER functionTimer(__FUNCTION__);
  if (argc != 2)
  {
    cout << " USAGE: " << argv[0] << " <cfg file> " << endl;
    exit(0);
  }

  SIMPLE_PARSER parser(argv[1]);
  string path = parser.getString("path", "./temp/");

  string distanceFilename = path + parser.getString("distance field", "dummy.field3D");
  string potentialFilename = path + parser.getString("potential field", "dummy.field3D");
  string pointsFilename = path + parser.getString("poisson output", "dummy.ply");

  relativeErrorThreshold = parser.getFloat("relative error threshold", relativeErrorThreshold);
  cout << " Using relative error threshold: " << relativeErrorThreshold << endl;

  distanceField = FIELD_3D(distanceFilename.c_str());
  //colorField = FIELD_3D(potentialFilename.c_str());

  // get the points
  readPointsPLY(pointsFilename.c_str(), plyPoints);

  /*
  FIELD_3D fieldCopy(distanceField);
  fieldCopy -= marchingSurface;
  triangleMesh = new TRIANGLE_MESH(fieldCopy);

  int potentialOffset = parser.getInt("potential offset", 3);

  cout << " Using a potential offset surface of: " << potentialOffset << endl;

  // create an offset distance field
  FIELD_3D offsetField(distanceField);
  offsetField -= potentialOffset * offsetField.dx();

  // get its closest point field
  VECTOR3_FIELD_3D closestPoints = VECTOR3_FIELD_3D::computeClosestPoints(offsetField);
  //closestPoints.write("./temp/closest.vector.field3d");

  // do an extension of the color field using the offset distance
  offsetColors = VECTOR3_FIELD_3D::setToClosestPointValues(colorField, closestPoints);
  //offsetColors.write("./temp/offset.field3d");

  FIELD_3D originalExtensionField = offsetColors;

  // create a shell potential
  FIELD_3D scaledDistance(offsetField);
  shellPotential = FIELD_3D(offsetColors);
  shellPotential.setZeroBorder();

  const int xRes = shellPotential.xRes();
  const int slabSize = shellPotential.slabSize();
  for (int z = 1; z < shellPotential.zRes() - 1; z++)
    for (int y = 1; y < shellPotential.yRes() - 1; y++)
      for (int x = 1; x < shellPotential.xRes() - 1; x++)
      {
        int index = x + y * xRes + z * slabSize;
        if (scaledDistance[index] < 0) 
        {
          shellPotential[index] = 0;
          continue;
        }

        if (scaledDistance[index + 1] < 0 ||
            scaledDistance[index - 1] < 0 ||
            scaledDistance[index + xRes] < 0 ||
            scaledDistance[index - xRes] < 0 ||
            scaledDistance[index + slabSize] < 0 ||
            scaledDistance[index - slabSize] < 0)
        {
          // do nothing, just store the location
          shellIndices.push_back(index);
          shellPositions.push_back(scaledDistance.cellCenter(index));
        }
        else
          shellPotential[index] = 0;
      }
  surfaceShell = distanceField.outerShell();

  for (int x = 0; x < surfaceShell.totalCells(); x++)
    if (surfaceShell[x] > 0.5)
    {
      surfaceIndices.push_back(x);
      surfacePositions.push_back(surfaceShell.cellCenter(x));
    }
    */
 
  //buildPoissonSourceSimulation();
  buildLejaPoints();

  writePlyPoints("./temp/leja.ply", topPoints);

  // write the roots out
  OPTIMIZE_3D optimize3D;
  //optimize3D.setScoreFields(distanceField, offsetColors);
  optimize3D.setScoreFields(distanceField, distanceField);
  vector<QUATERNION> topRoots;
  for (unsigned int x = 0; x < topSpheres.size(); x++)
  {
    VEC3F center = topSpheres[x].translation();
    topRoots.push_back(QUATERNION(center[0], center[1], center[2], 0));
  }
  POLYNOMIAL_4D topPolynomial(topRoots);
  for (unsigned int x = 0; x < topSpheres.size(); x++)
    topPolynomial.powersMutable()[x] = 1.0;
  optimize3D.top() = topPolynomial;

  /*
  vector<QUATERNION> bottomRoots;
  for (unsigned int x = 0; x < bottomSpheres.size(); x++)
  {
    VEC3F center = bottomSpheres[x].translation();
    bottomRoots.push_back(QUATERNION(center[0], center[1], center[2], 0));
  }
  POLYNOMIAL_4D bottomPolynomial(bottomRoots);
  for (unsigned int x = 0; x < bottomWeights.size(); x++)
    bottomPolynomial.powersMutable()[x] = fabs(bottomWeights[x]);
  optimize3D.bottom() = bottomPolynomial;
  */

  // add a zero root to the front
  cout << " ADDING A ZERO ROOT " << endl;
  QUATERNION zeroRoot(0,0,0,0);
  optimize3D.top().addFrontRoot(zeroRoot);
  /*
  optimize3D.bottom().addFrontRoot(zeroRoot);
  */

  string outputFilename = path + parser.getString("source simulation filename", "temp.o3d");
  cout << " Outputting to file " << outputFilename.c_str() << endl;
  optimize3D.write(outputFilename.c_str());

  return 0;
}
