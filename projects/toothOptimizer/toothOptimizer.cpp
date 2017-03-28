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

#include "SETTINGS.h"
#include "TAO_SOLVER.h"
#include "BLMVM.h"
#include "MORE_THUENTE.h"

#include "SIMPLE_PARSER.h"
#include "OPTIMIZE_3D.h"
#include "TRIANGLE_MESH.h"
#include "TIMER.h"

string path;
string cfgFilename;
int res = 25;
int maxTaoIterations = 100;
OPTIMIZE_3D optimize3D;
Real fractalRadius = 4.0;

Real maxScore = -1;
vector<Real> scores;

TRIANGLE_MESH* triangleMesh = NULL;
TRIANGLE_MESH* distanceMesh = NULL;
TRIANGLE_MESH* fractalMesh = NULL;

////////////////////////////////////////////////////////////////////////////////////////
// monitor function that dumps out the current optimization state
////////////////////////////////////////////////////////////////////////////////////////
void monitor()
{
  static int counter = 0;

  if (counter % 10 == 0)
  {
    char buffer[256];
    //sprintf(buffer, "%soutput.%04i.o3d", path.c_str(), counter);
    sprintf(buffer, "%soutput.o3d", path.c_str());
    optimize3D.write(buffer);
  }

  counter++;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void FormScaledFunctionGradient(const VECTOR& state, Real& function, VECTOR& gradient)
{
  int scalingStart = 1 + optimize3D.top().totalRoots() + optimize3D.bottom().totalRoots();
  int solutionSize = 1 + optimize3D.top().totalRoots() + optimize3D.bottom().totalRoots() + 2;

  // copy values back from the solver
  optimize3D.expScaling() = exp(state[0]);
  int xIndex = 1;
  for (int i = 0; i < optimize3D.bottom().totalRoots(); i++, xIndex++)
    optimize3D.bottom().powersMutable()[i] = state[xIndex];
  for (int i = 0; i < optimize3D.top().totalRoots(); i++, xIndex++)
    optimize3D.top().powersMutable()[i] = state[xIndex];
  optimize3D.top().powerScalar() = state[scalingStart];
  optimize3D.bottom().powerScalar() = state[scalingStart + 1];

  // compute the new score
  function = optimize3D.computeScoreAP();
  //cout.precision(16);
  cout << "." << flush;
  //cout << function << endl;

  gradient.resize(solutionSize);
  gradient[0] = optimize3D.computeExpGradientAP();

  TIMER bottomTimer("Bottom Gradient");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int i = 0; i < optimize3D.bottom().totalRoots(); i++)
    gradient[i + 1] = optimize3D.computeBottomGradientAP(i);
  bottomTimer.stop();

  TIMER topTimer("Top Gradient");
  const int offset = 1 + optimize3D.bottom().totalRoots();
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int i = 0; i < optimize3D.top().totalRoots(); i++)
  {
    const int xIndex = i + offset;
    gradient[xIndex] = optimize3D.computeTopGradientAP(i);
  }
  topTimer.stop();

  gradient[scalingStart] = optimize3D.computeTopBulkScalarGradientAP();
  gradient[scalingStart + 1] = optimize3D.computeBottomBulkScalarGradientAP();
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void FormScaledOnlyFunctionGradient(const VECTOR& state, Real& function, VECTOR& gradient)
{
  // copy values back from the solver
  optimize3D.expScaling() = exp(state[0]);
  optimize3D.top().powerScalar() = state[1];
  optimize3D.bottom().powerScalar() = state[2];

  // compute the new score
  function = optimize3D.computeScoreAP();
  //cout.precision(16);
  cout << "." << flush;
  //cout << function  << endl;

  gradient.resize(3);
  gradient[0] = optimize3D.computeExpGradientAP();
  gradient[1] = optimize3D.computeTopBulkScalarGradientAP();
  gradient[2] = optimize3D.computeBottomBulkScalarGradientAP();
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int callConstrainedScaledOnlyTao(OPTIMIZE_3D& optimize3D)
{
  TIMER functionTimer(__FUNCTION__);
  int solutionSize = 3;

  vector<Real> values(3);
  values[0] = log(optimize3D.expScaling());
  values[1] = optimize3D.top().powerScalar();
  values[2] = optimize3D.bottom().powerScalar();

  VECTOR solution(solutionSize);
  for (int x = 0; x < solutionSize; x++)
    solution[x] = values[x];
  VECTOR upperBound = VECTOR::Constant(solution.size(), FLT_MAX);
  VECTOR lowerBound = VECTOR::Zero(solution.size());
  lowerBound[0] = -FLT_MAX;

  //cout << " New Tolerances: " << fatol << " " << frtol << " " << gatol << " " << grtol << " " << gttol << endl;
  Real fatol = 1e-8;
  Real frtol = 1e-8;
  Real gatol = 1e-8;
  Real grtol = 1e-8;
  Real gttol = 0;

  TAO_SOLVER taoSolver;
  taoSolver.SetTolerances(fatol, frtol, gatol, grtol, gttol);
  taoSolver.maxScore() = maxScore;
  taoSolver.optimizeBLMVM(solution, upperBound, lowerBound, FormScaledOnlyFunctionGradient, monitor);

  optimize3D.expScaling() = exp(solution[0]);
  optimize3D.top().powerScalar() = solution[1];
  optimize3D.bottom().powerScalar() = solution[2];

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int callConstrainedScaledTao(OPTIMIZE_3D& optimize3D)
{
  TIMER functionTimer(__FUNCTION__);
  int scalingStart = 1 + optimize3D.top().totalRoots() + optimize3D.bottom().totalRoots();
  int solutionSize = 1 + optimize3D.top().totalRoots() + optimize3D.bottom().totalRoots() + 2;
  Real values[solutionSize];

  values[0] = log(optimize3D.expScaling());
  int xIndex = 1;
  for (int i = 0; i < optimize3D.bottom().totalRoots(); i++, xIndex++)
    values[xIndex] = optimize3D.bottom().powersMutable()[i];
  for (int i = 0; i < optimize3D.top().totalRoots(); i++, xIndex++)
    values[xIndex] = optimize3D.top().powersMutable()[i];

  values[scalingStart] = optimize3D.top().powerScalar();
  values[scalingStart + 1] = optimize3D.bottom().powerScalar();

  VECTOR solution(solutionSize);
  for (int x = 0; x < solutionSize; x++)
    solution[x] = values[x];
  VECTOR upperBound = VECTOR::Constant(solution.size(), FLT_MAX);
  VECTOR lowerBound = VECTOR::Zero(solution.size());
  lowerBound[0] = -FLT_MAX;

  Real fatol = 1e-8;
  Real frtol = 1e-8;
  Real gatol = 1e-8;
  Real grtol = 1e-8;
  Real gttol = 0;

  TAO_SOLVER taoSolver;
  taoSolver.SetTolerances(fatol, frtol, gatol, grtol, gttol);
  taoSolver.maxScore() = maxScore;
  taoSolver.optimizeBLMVM(solution, upperBound, lowerBound, FormScaledFunctionGradient, monitor);

  // copy values back from the solver
  optimize3D.expScaling() = exp(solution[0]);
  xIndex = 1;
  for (int i = 0; i < optimize3D.bottom().totalRoots(); i++, xIndex++)
    optimize3D.bottom().powersMutable()[i] = solution[xIndex];
  for (int i = 0; i < optimize3D.top().totalRoots(); i++, xIndex++)
    optimize3D.top().powersMutable()[i] = solution[xIndex];
  optimize3D.top().powerScalar() = solution[scalingStart];
  optimize3D.bottom().powerScalar() = solution[scalingStart + 1];

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int solveConstrainedScaledOnlyTao(int argc, char** argv, OPTIMIZE_3D& optimize3D)
{
  TIMER functionTimer(__FUNCTION__);
  callConstrainedScaledOnlyTao(optimize3D);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int solveConstrainedScaledTao(int argc, char** argv, OPTIMIZE_3D& optimize3D)
{
  TIMER functionTimer(__FUNCTION__);
  callConstrainedScaledTao(optimize3D);
  return 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void computeFractal(int res)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing fractal with res " << res << " ... " << flush;

  optimize3D.maxIterations() = 1;
  optimize3D.computeLogScaledPowerRationalMap(true);

  cout << " and max iterations: " << optimize3D.maxIterations() << "..." << flush;
 
  FIELD_3D fieldCopy(optimize3D.fractal());
  fieldCopy.clampInfsToMean();
  fieldCopy.setBorders(fieldCopy.totalCells());

  triangleMesh = new TRIANGLE_MESH(fieldCopy);
  cout << " Marching cubes found " << triangleMesh->triangles().size() << " triangles " << endl;

  if (triangleMesh->triangles().size() == 0)
  {
    FIELDVIEW3D(fieldCopy);
  }
  VEC3F minBox, maxBox;
  triangleMesh->boundingBox(minBox, maxBox);
  cout << " final fractal bounding box: " << minBox << " " << maxBox << endl;

  if (optimize3D.distanceField().totalCells() > 0)
  {
    FIELD_3D distanceField = optimize3D.distanceField();

    if (distanceMesh) delete distanceMesh;
    distanceMesh = new TRIANGLE_MESH(distanceField);
  }
  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void buildReweightedInsideField()
{
  cout << " ========================================================================================= " << endl;
  cout << "  BUILDING REWEIGHTED INSIDE ONLY SCORING FUNCTION " << endl;
  cout << " ========================================================================================= " << endl;
  
  // build the triangle mesh so we can build an SDF
  computeFractal(97);
  triangleMesh->writeOBJ("temp.obj");

  VEC3F minBox, maxBox;
  triangleMesh->boundingBox(minBox, maxBox);

  cout << " Current fractal's triangle mesh bounding box, min: " << minBox << " max: " << maxBox << endl; 

  char buffer[256];
  sprintf(buffer, "./bin/sdfGen %s 1", cfgFilename.c_str());
  cout << " Calling: " << buffer << endl;
  system(buffer);
  FIELD_3D fractalDistance("temp.distance.field3d");

  cout << " Using fractal radius: " << fractalRadius << endl;
  fractalDistance.center() += VEC3F(-0.5, -0.5, -0.5);
  fractalDistance.scaleLengths(fractalRadius);
  fractalDistance.center() *= fractalRadius;

  fractalMesh = new TRIANGLE_MESH(fractalDistance);

  const FIELD_3D& meshDistance = optimize3D.distanceField();
  FIELD_3D poorlyCaptured(meshDistance);
  for (int x = 0; x < poorlyCaptured.totalCells(); x++)
    poorlyCaptured[x] = (fractalDistance[x] > 0 && meshDistance[x] < 0) ? 1 : 0;

  FIELD_3D wellCaptured(meshDistance);
  for (int x = 0; x < wellCaptured.totalCells(); x++)
    wellCaptured[x] = (fractalDistance[x] < 0 && meshDistance[x] < 0) ? 1 : 0;

  FIELD_3D reweighted(optimize3D.weightField());

  // See what the regions sum to
  Real wellSum = 0;
  Real poorSum = 0;
  for (int x = 0; x < reweighted.totalCells(); x++)
    if (meshDistance[x] < 0)
    {
      if (poorlyCaptured[x] > 0.5)
        poorSum += reweighted[x];
      if (wellCaptured[x] > 0.5)
        wellSum += reweighted[x];
    }
  
  cout << " well-captured sum: " << wellSum << endl;
  cout << " poorly-captured sum: " << poorSum << endl;
  wellSum = fabs(wellSum);
  poorSum = fabs(poorSum);

  // how many outside cells are there?
  int totalOutside = 0;
  for (int x = 0; x < fractalDistance.totalCells(); x++)
    if (meshDistance[x] > 0)
      totalOutside++;

  // how much is the entire side worth
  Real totalWorth = 0.5;

  // how much do you want the outside to be worth?
  Real outsideWorth = 0.5;
  Real outsidePerCell = ((outsideWorth * totalWorth) / (1.0 - outsideWorth)) / totalOutside;
  cout << " Setting ALL outside cells to: " << outsidePerCell << endl;

  // what fraction of total worth do you want to allocate to the poorly fit regions
  Real poorFraction = 0.5;
  Real wellFraction = 1.0 - poorFraction;

  Real poorReweight = totalWorth * poorFraction / poorSum;
  Real wellReweight = totalWorth * wellFraction / wellSum;

  // reweight based on region sums
  //
  // The 0.25 factors show up because the inside is worth
  // 50% of the score, so now we're redistributing that score
  // to half of the well-captured regions, and half to the poorly-captured regions.
  for (int x = 0; x < reweighted.totalCells(); x++)
    if (meshDistance[x] < 0)
    {
      if (poorlyCaptured[x] > 0.5)
        reweighted[x] *= poorReweight;
      if (wellCaptured[x] > 0.5)
        reweighted[x] *= wellReweight;
    }
    else
    { 
      // set the outside to some flat values. There still needs to be some kind of penalty here so
      // that the fractal doesn't expand to fill all of space, but concentrating repulsor penalty right by the
      // boundary is counter-productive.
      reweighted[x] = outsidePerCell;
    }

  cout << " Original weight dims: " << optimize3D.weightField().center() << " " << optimize3D.weightField().lengths() << endl;
  cout << " Reweight dims:        " << reweighted.center() << " " << reweighted.lengths() << endl;

  // commit the new weight field
  optimize3D.weightField() = reweighted;
}

//////////////////////////////////////////////////////////////////////////////
// run a solve for a specific resolution
//////////////////////////////////////////////////////////////////////////////
void optimizeSpecificResolution(int argc, char* argv[], string path, SIMPLE_PARSER& parser, int res)
{
  cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;
  cout << " OPTIMIZING FOR RESOLUTION: " << res << endl;
  cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;

  // clear any existing score runs
  scores.clear();

  // print out the initial root powers
  optimize3D.printRootPowers();

  VEC3F center = optimize3D.fractal().center();
  VEC3F lengths = optimize3D.fractal().lengths();
  optimize3D.fractal() = FIELD_3D(res,res,res, center, lengths);
  optimize3D.auxiliary() = FIELD_3D(res,res,res, center, lengths);
 
  optimize3D.maxIterations() = 1;

  if (res == 25)
  {
    optimize3D.estimateConformalRadius();

    computeFractal(res);
    TIMER::printTimings();
  }
  Real score = optimize3D.computeLogScaledPowerRationalScore();
  cout << " Initial score is: " << score << endl;

  // get what the max score is so we can normalize to -1
  maxScore = optimize3D.computeMaximumScaledPowerRationalScore();
  cout << " Maximum possible score is: " << maxScore << endl;

  // for the very first solve, do the scalar-only solve
  solveConstrainedScaledOnlyTao(argc, argv, optimize3D);

  // optimize over everybody
  solveConstrainedScaledTao(argc, argv, optimize3D);
 
  // now do a reweighted run
  FIELD_3D originalWeights = optimize3D.weightField();
  buildReweightedInsideField();
  
  maxScore = optimize3D.computeMaximumScaledPowerRationalScore();
  cout << " Maximum possible score is now: " << maxScore << endl;
  cout << " ========================================================================================= " << endl;
  cout << "  RUNNING REWEIGHTED OPTIMIZATION " << endl;
  cout << " ========================================================================================= " << endl;
  solveConstrainedScaledTao(argc, argv, optimize3D);

  // do another optimization over everybody for good measure
  optimize3D.weightField() = originalWeights;
  maxScore = optimize3D.computeMaximumScaledPowerRationalScore();
  cout << " Maximum possible score is now: " << maxScore << endl;
  cout << " ========================================================================================= " << endl;
  cout << "  RE-RUNNING UNWEIGHTED OPTIMIZATION " << endl;
  cout << " ========================================================================================= " << endl;
  solveConstrainedScaledTao(argc, argv, optimize3D);

  // print out the root powers, just to see if it all checks out
  optimize3D.printRootPowers();

  // write out an optimized version
  string optimizedFilename = path + parser.getString("optimization output filename", "temp.o3d");

  // append the current level to the filename
  char buffer[256];
  sprintf(buffer, "_res_%i.o3d", res);
  optimizedFilename = optimizedFilename + string(buffer);

  cout << " Outputting final result to " << optimizedFilename.c_str() << endl;
  optimize3D.write(optimizedFilename);
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  TIMER functionTimer(__FUNCTION__);
  if (argc < 2)
  {
    cout << " USAGE: " << argv[0] << " <cfg file> <restart file (optional)>" << endl;
    return 0;
  }

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

  cfgFilename = string(argv[1]);
  SIMPLE_PARSER parser(argv[1]);
  path = parser.getString("path", "./temp/");

  string inputFilename = path + parser.getString("source simulation filename", "dummy.o3d");

  res = parser.getInt("optimization res", res);
  cout << " Using optimization resolution " << res << "^3" << endl;

  maxTaoIterations = parser.getInt("max TAO iterations", maxTaoIterations);
  cout << " Using maximum TAO iterations " << maxTaoIterations << endl;
  cout << " Reading in source simulation result " << inputFilename.c_str() << endl;
  optimize3D.read(inputFilename.c_str());

  // get the radius that we want the fractal to be
  fractalRadius = parser.getFloat("fractal radius", 4.0);
  assert(fractalRadius > 0);

  if (optimize3D.lengths()[0] < 0.95 * fractalRadius)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " NOT TRANSLATING FIRST ROOT" << endl;
    optimize3D.translateEverythingExceptFirst(VEC3F(-0.5, -0.5, -0.5));
    optimize3D.scaleEverything(fractalRadius);
  }

  Real amp = -1;
  amp = parser.getFloat("initial bulk scalar amplitude", amp);
  cout << " Using initial bulk scalar amplitude " << amp << endl;

  POLYNOMIAL_4D& top = optimize3D.top();
  POLYNOMIAL_4D& bottom = optimize3D.bottom();

  for (int x = 1; x < top.totalRoots(); x++)
    top.powersMutable()[x] *= amp;
  for (int x = 1; x < bottom.totalRoots(); x++)
    bottom.powersMutable()[x] *= amp;
  
  top.powerScalar() = 1;
  bottom.powerScalar() = 1;

  optimizeSpecificResolution(argc, argv, path, parser, 25);
  optimizeSpecificResolution(argc, argv, path, parser, 50);
  optimizeSpecificResolution(argc, argv, path, parser, 97);

  // run headless
  return 0;
}
