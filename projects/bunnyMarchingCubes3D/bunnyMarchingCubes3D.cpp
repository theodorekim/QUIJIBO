// OpenGL Graphics includes
#include <GLUT/glut.h>
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
//#include "glvu.h"
#include "SIMPLE_PARSER.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Meshes and integrators
//////////////////////////////////////////////////////////////////////////////
TRIANGLE_MESH* triangleMesh = NULL;
Real marchingSurface = 0;

OPTIMIZE_3D optimize3D;
string outputPrefix("./temp/temp");

float minValue;
//float transformValue = 1.41;
float transformValue = 0;
float maxValue;
int steps;
int maxIterations = 1;
int res = 97;

Real isosurface = 0;

#define USING_LOWMEMORY 1

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void computeFractal(int res)
{
  TIMER functionTimer(__FUNCTION__);
  optimize3D.maxIterations() = maxIterations;
#if !USING_LOWMEMORY
  cout << " Computing fractal with res " << res << " ... " << flush;

  optimize3D.computeLogScaledPowerRationalMap(true);
  TIMER::printTimings();
  system("purge");

  cout << " and max iterations: " << optimize3D.maxIterations() << "..." << flush;

  cout << " overflows: " << (int)optimize3D.overflowed().sum() << endl;
  cout << " infs: " << optimize3D.fractal().totalInfs() << endl;
  cout << " NaNs: " << optimize3D.fractal().totalNans() << endl;
#endif
 
#if 0
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Current score: " << optimize3D.computeLogScaledPowerRationalScore() << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
#endif

  //string outputField3D = outputPrefix+ string(".field3d");
  //optimize3D.fractal().write(outputField3D);

  if (triangleMesh) delete triangleMesh;
#if USING_LOWMEMORY
  triangleMesh = new TRIANGLE_MESH(optimize3D.fractal().center(),
                                   optimize3D.fractal().lengths(),
                                   VEC3I(res,res,res),
                                   optimize3D.top(), optimize3D.bottom(), optimize3D.expScaling(), optimize3D.maxIterations(), optimize3D.slice(), isosurface);
#else
  triangleMesh = new TRIANGLE_MESH(optimize3D.fractal(), optimize3D.top(), optimize3D.bottom(), optimize3D.expScaling(), optimize3D.maxIterations(), optimize3D.slice(), isosurface);
#endif
  cout << " Marching cubes found " << triangleMesh->triangles().size() << " triangles " << endl;

  //VEC3F minBox, maxBox;
  //triangleMesh->boundingBox(minBox, maxBox);
  //cout << " final fractal bounding box: " << minBox << " " << maxBox << endl;

  TIMER::printTimings();
}

//////////////////////////////////////////////////////////////////////////////
// TOOTH example:
//
// interpolating between 30 and 0 gives interesting results. A "cocoon" shape
// emerges, where there is both an inner and outer shell. Maybe interesting
// to render as glass?
//
// BUNNY example:
//
// 7000 steps:
// at 2 iterations, 0 to -2400 looks okay
// at 3 iterations, function goes from 0 to -4000 in space of one internal
//   grid cell
// 
// 5000 steps:
// at 3 iterations, 5000 steps, -400 to 0 might be interesting
// at 3 iterations, 5000 steps, 0 to 10 is interesting, AND fits inside
//   the existing grid bounds. Does the same cocoon thing as the tooth.
///////////////////////////////////////// /////////////////////////////////////
void marchIsosurfaces(const string& descriptor)
{
  float delta = (maxValue - minValue) / (steps - 1);
  
  // for the TOOTH example: extend the size to capture the entire matroyshka shape
  VEC3F lengths = optimize3D.fractal().lengths();
  lengths *= 1.15;
  optimize3D.fractal().setLengths(lengths);

  // go ahead and compute the fractal
  computeFractal(res);
  string outputField3D = outputPrefix + string(".isosurfaces.field3d");
  optimize3D.fractal().write(outputField3D);

  for (int x = 0; x < steps; x++)
  {
    isosurface = minValue + x * delta;
    cout << " Extracing isosurface: " << isosurface << endl;

    FIELD_3D fractal = optimize3D.fractal();
    fractal -= isosurface;

    if (triangleMesh) delete triangleMesh;
    triangleMesh = new TRIANGLE_MESH(fractal, optimize3D.top(), optimize3D.bottom(), optimize3D.expScaling(), optimize3D.maxIterations(), optimize3D.slice(), isosurface);
    cout << " Marching cubes found " << triangleMesh->triangles().size() << " triangles " << endl;

    // write out the OBJ
    char buffer[256];
    sprintf(buffer, ".res.%i.iterations.%i_%s.%04i.obj", res, maxIterations, descriptor.c_str(), x);
    string outputFilename = outputPrefix + string(buffer);
    //cout << " Writing OBJ " << outputFilename.c_str() << endl;
    triangleMesh->writeOBJ(outputFilename.c_str());
  }
}

/////////////////////////////////////////////////////////////////////////////
// BUNNY: -0.65, at res 50, something shows up.
//        -0.7, at res 100, nothing shows up.
//        Perhaps best to do -0.65 to 0 or 0 to 0.65
//
//        0 to 0.65 looks fine at res = 50
//////////////////////////////////////////////////////////////////////////////
void marchSlices()
{
  float delta = (maxValue - minValue) / (steps - 1);

  for (int x = 0; x < steps; x++)
  {
    // set the slice value
    Real slice = minValue + x * delta;
    cout << " Using quaternion slice value: " << slice << endl;
    optimize3D.slice() = slice;
    
    // go ahead and compute the fractal
    computeFractal(res);

    // write out the field
    char buffer[256];
    sprintf(buffer, ".res.%i.iterations.%i_slice.%04i.field3d", res, maxIterations, x);
    string outputFieldname = outputPrefix + string(buffer);
    //cout << " Writing to field file " << outputFieldname.c_str() << endl;
    //optimize3D.fractal().write(outputFieldname);

    // write out the OBJ
    sprintf(buffer, ".res.%i.iterations.%i_slice.%04i.obj", res, maxIterations, x);
    string outputFilename = outputPrefix + string(buffer);
    //cout << " Writing OBJ " << outputFilename.c_str() << endl;
    triangleMesh->writeOBJ(outputFilename.c_str());
  }
}

//////////////////////////////////////////////////////////////////////////////
// TOOTH example:
//
// for x, from 0.5 to 1.5, went from just a featureless tube to just a few
// isolated spheres
//
// for xyz from 0.5 to 1.5, 
//    at 0.5, the shape is entirely just a sphere
//    after 1.25, zero triangle found
//    by 1.0625, not a lot of triangles left
//    at 0.5, perhaps more iterations give more features?
//
//    between 1 and 1.0625 seems like the way to go
//
//    jumping from 2 to 4 iterations didn't make a big difference
//
// BUNNY example
//
// for x, 0.5 to 1.2 is interesting
//
// for xyz 0.15 to 1.1 is interesting, need to expand box size for 0.25
//////////////////////////////////////////////////////////////////////////////
void marchScaling(const VEC3F& axis, const string& descriptor)
{
  float delta = (maxValue - minValue) / (steps - 1);
  cout << " Using the scaling axis: " << axis << " delta: " << delta << endl;

  // back up the original root positions
  const POLYNOMIAL_4D topOriginal = optimize3D.top();
  const POLYNOMIAL_4D bottomOriginal = optimize3D.bottom();

  for (int x = 0; x < steps; x++)
  {
    Real scaling = minValue + x * delta;

    MATRIX3 scale = MATRIX3::I();

    if (axis[0] > 0)
      scale(0,0) = scaling;
    if (axis[1] > 0)
      scale(1,1) = scaling;
    if (axis[2] > 0)
      scale(2,2) = scaling;

    cout << " Scaling matrix: " << scale << endl;

    for (int i = 0; i < topOriginal.roots().size(); i++)
    {
      QUATERNION q0 = topOriginal.roots()[i];
      VEC3F v0(q0[0], q0[1], q0[2]);
      v0 = scale * v0;
      optimize3D.top().rootsMutable()[i] = QUATERNION(v0[0], v0[1], v0[2], 0);
    }
    for (int i = 0; i < bottomOriginal.roots().size(); i++)
    {
      QUATERNION q0 = bottomOriginal.roots()[i];
      VEC3F v0(q0[0], q0[1], q0[2]);
      v0 = scale * v0;
      optimize3D.bottom().rootsMutable()[i] = QUATERNION(v0[0], v0[1], v0[2], 0);
    }
    
    // go ahead and compute the fractal
    computeFractal(res);

    // write out the field
    char buffer[256];
    sprintf(buffer, ".res.%i.iterations.%i_%s.%04i.field3d", res, maxIterations, descriptor.c_str(), x);
    string outputFieldname = outputPrefix + string(buffer);
    //cout << " Writing to field file " << outputFieldname.c_str() << endl;
    //optimize3D.fractal().write(outputFieldname);

    // write out the OBJ
    sprintf(buffer, ".res.%i.iterations.%i_%s.%04i.obj", res, maxIterations, descriptor.c_str(), x);
    string outputFilename = outputPrefix + string(buffer);
    //cout << " Writing OBJ " << outputFilename.c_str() << endl;
    triangleMesh->writeOBJ(outputFilename.c_str());
  }
}

//////////////////////////////////////////////////////////////////////////////
// rotations by any axes don't seem to have a major influence.
// they may warp a bit, but not really enough to be that noticeable.
//////////////////////////////////////////////////////////////////////////////
void marchRotations(const VEC3F& axis, const string& descriptor)
{
  float delta = (maxValue - minValue) / (steps - 1);
  cout << " Using the rotation axis: " << axis << " delta: " << delta << endl;

  // back up the original root positions
  const POLYNOMIAL_4D topOriginal = optimize3D.top();
  const POLYNOMIAL_4D bottomOriginal = optimize3D.bottom();

  for (int x = 0; x < steps; x++)
  {
    Real angle = x * delta;

    // convert angle to radians
    Real theta = (angle / 360.0) * 2.0 * M_PI;

    MATRIX3 rotation = MATRIX3::rotation(axis, theta);
    for (int i = 0; i < topOriginal.roots().size(); i++)
    {
      QUATERNION q0 = topOriginal.roots()[i];
      VEC3F v0(q0[0], q0[1], q0[2]);
      v0 = rotation * v0;
      optimize3D.top().rootsMutable()[i] = QUATERNION(v0[0], v0[1], v0[2], 0);
    }
    for (int i = 0; i < bottomOriginal.roots().size(); i++)
    {
      QUATERNION q0 = bottomOriginal.roots()[i];
      VEC3F v0(q0[0], q0[1], q0[2]);
      v0 = rotation * v0;
      optimize3D.bottom().rootsMutable()[i] = QUATERNION(v0[0], v0[1], v0[2], 0);
    }
    
    // go ahead and compute the fractal
    computeFractal(res);

    // write out the field
    char buffer[256];
    sprintf(buffer, ".res.%i.iterations.%i_%s.%04i.field3d", res, maxIterations, descriptor.c_str(), x);
    string outputFieldname = outputPrefix + string(buffer);
    //cout << " Writing to field file " << outputFieldname.c_str() << endl;
    //optimize3D.fractal().write(outputFieldname);

    // write out the OBJ
    sprintf(buffer, ".res.%i.iterations.%i_%s.%04i.obj", res, maxIterations, descriptor.c_str(), x);
    string outputFilename = outputPrefix + string(buffer);
    //cout << " Writing OBJ " << outputFilename.c_str() << endl;
    triangleMesh->writeOBJ(outputFilename.c_str());
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void marchConformal(const string& descriptor)
{
  float delta = (maxValue - minValue) / (steps - 1);

  cout << " Initial conformal radius: " << optimize3D.expScaling() << endl;
  cout << " Computing initial range ... " << endl;

  optimize3D.expScaling() = 0.1;

  for (int x = 0; x < steps; x++)
  {
    optimize3D.expScaling() *= 10;
    cout << " New scaling: " << optimize3D.expScaling() << endl;

    // go ahead and compute the fractal
    computeFractal(res);

    // write out the field
    char buffer[256];
    sprintf(buffer, "res.%i.iterations.%i_%s.%04i.field3d", res, maxIterations, descriptor.c_str(), x);
    string outputFieldname = outputPrefix + string(buffer);
    //cout << " Writing to field file " << outputFieldname.c_str() << endl;
    //optimize3D.fractal().write(outputFieldname);

    // write out the OBJ
    sprintf(buffer, "res.%i.iterations.%i_%s.%04i.obj", res, maxIterations, descriptor.c_str(), x);
    string outputFilename = outputPrefix + string(buffer);
    //cout << " Writing OBJ " << outputFilename.c_str() << endl;
    triangleMesh->writeOBJ(outputFilename.c_str());
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void marchZeroRoot(const VEC3F& direction, const string& descriptor)
{
  float delta = (maxValue - minValue) / (steps - 1);
  VEC3F      vDelta = direction * (Real)delta;

  //QUATERNION qDelta(0, 0, delta, 0);
  QUATERNION qDelta(vDelta[0], vDelta[1], vDelta[2], 0);
  cout << " Using the zero-only translation delta: " << delta <<  " quaternion version: " << qDelta << endl;

  const VEC3F centerOriginal = optimize3D.fractal().center();

  // back up the original root positions
  const POLYNOMIAL_4D topOriginal = optimize3D.top();
  const POLYNOMIAL_4D bottomOriginal = optimize3D.bottom();

  for (int x = 0; x < steps; x++)
  {
    optimize3D.top().rootsMutable()[0] = topOriginal.roots()[0] + x * qDelta;
    optimize3D.bottom().rootsMutable()[0] = bottomOriginal.roots()[0] + x * qDelta;
    
    //optimize3D.fractal().center() = centerOriginal + (Real)x * vDelta;

    cout << " New delta: " << x * qDelta << endl;

    // go ahead and compute the fractal
    computeFractal(res);

    // write out the field
    char buffer[256];
    sprintf(buffer, ".res.%i.iterations.%i_%s.%04i.field3d", res, maxIterations, descriptor.c_str(), x);
    string outputFieldname = outputPrefix + string(buffer);
    //cout << " Writing to field file " << outputFieldname.c_str() << endl;
    //optimize3D.fractal().write(outputFieldname);

    // write out the OBJ
    sprintf(buffer, ".res.%i.iterations.%i_%s.%04i.obj", res, maxIterations, descriptor.c_str(), x);
    string outputFilename = outputPrefix + string(buffer);
    //cout << " Writing OBJ " << outputFilename.c_str() << endl;
    triangleMesh->writeOBJ(outputFilename.c_str());
  }
}

//////////////////////////////////////////////////////////////////////////////
// For the tooth: 
//    0 to 1.0 along (0,0,1,0) gives interesting results
//
//    0 to -1.0 along (0,0,1,0) also gives interesting results,
//      not the same as 0 to 1.0
//
//    0 to 2.0 along (0,1,0,0) makes it disappear eventually, though
//      it doesn't seem as interesting 
//
//    0 to -0.75 along (0,1,0,0) also makes it disappear, and it gets interesting
//
//    0 to 1.0 along (1,0,0,0) makes it disappear
//
// BUNNY:
//    0 to 2.0 along (1,0,0,0) makes it disappear
//    0 to 2.0 along (0,1,0,0) makes it disappear, even more quickly than x
//    0 to 2.0 along (0,0,1,0) makes it disappear, also very quickly
//////////////////////////////////////////////////////////////////////////////
void marchTranslations(const VEC3F& translation, const string& descriptor)
{
  float delta = transformValue;
  //VEC3F      vDelta(0, 0, delta);
  //VEC3F      vDelta(0, delta, 0);
  //VEC3F      vDelta(delta, 0, 0);
  //VEC3F      vDelta = direction * (Real)delta;
  VEC3F      vDelta = translation;

  //QUATERNION qDelta(0, 0, delta, 0);
  QUATERNION qDelta(vDelta[0], vDelta[1], vDelta[2], 0);
  cout << " Using the translation delta: " << delta <<  " quaternion version: " << qDelta << endl;

  const VEC3F centerOriginal = optimize3D.fractal().center();

  // back up the original root positions
  const POLYNOMIAL_4D topOriginal = optimize3D.top();
  const POLYNOMIAL_4D bottomOriginal = optimize3D.bottom();

  //for (int x = 0; x < steps; x++)
  {
    for (int i = 0; i < topOriginal.roots().size(); i++)
      optimize3D.top().rootsMutable()[i] = topOriginal.roots()[i] + qDelta;
    for (int i = 0; i < bottomOriginal.roots().size(); i++)
      optimize3D.bottom().rootsMutable()[i] = bottomOriginal.roots()[i] + qDelta;
    
    optimize3D.fractal().center() = centerOriginal + vDelta;

    cout << " New delta: " << qDelta << endl;

    // go ahead and compute the fractal
    computeFractal(res);

    // write out the field
    char buffer[256];
    sprintf(buffer, ".res.%i.iterations.%i_%s.field3d", res, maxIterations, descriptor.c_str());
    string outputFieldname = outputPrefix + string(buffer);
    //cout << " Writing to field file " << outputFieldname.c_str() << endl;
    //optimize3D.fractal().write(outputFieldname);

    // write out the OBJ
    sprintf(buffer, ".res.%i.iterations.%i_%s.obj", res, maxIterations, descriptor.c_str());
    string outputFilename = outputPrefix + string(buffer);
    //cout << " Writing OBJ " << outputFilename.c_str() << endl;
    triangleMesh->writeOBJ(outputFilename.c_str());
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int mainOld(int argc, char* argv[])
{
  string cfg(argv[1]);
  SIMPLE_PARSER parser(cfg);

  VEC3 parsed = parser.getVector3("root translation", VEC3(0,0,0));
  exit(0);

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
  if (argc != 5)
  {
    cout << " USAGE: " << argv[0] << " <optimize 3D file> <res> <max iterations> <output prefix>" << endl;
    cout << " Transform types are: slices translation" << endl;
    return 0;
  }

  /*
  // which transformation are we doing?
  string transform(argv[5]);
  transformValue = atof(argv[6]); 
  cout << " Using transform type: " << transform.c_str() << endl;
  //cout << " Min: " << minValue << " Max: " << maxValue << " Steps: " << steps << endl;
  cout << " Transform value: " << transformValue << endl;
  */

  optimize3D.read(argv[1]);

  //Real slice = 0.25;
  //cout << " Using slice value: " << slice << endl;
  //optimize3D.slice() = slice;

  if (argc >= 3)
    res = atoi(argv[2]);
  cout << " Using resolution: " << res << endl;

  if (argc >= 4)
    maxIterations = atoi(argv[3]);
  cout << " Using max iterations: " << maxIterations << endl;
 
  if (argc >= 5)
    outputPrefix = string(argv[4]);
  cout << " Using output prefix: " << outputPrefix.c_str() << endl;

  VEC3F& center = optimize3D.fractal().center();
  VEC3F& lengths = optimize3D.fractal().lengths();

  cout << " Original center: " << center << endl;
  cout << " Original lengths: " << lengths << endl;

  // first zoom
  //center = VEC3F(-0.953853, 1.096036, 2.767263 - transformValue);
  //lengths *= 0.5;

  // second zoom
  //center = VEC3F(-0.923137, 1.171546, 2.695665 - transformValue);
  //lengths *= 0.1;

  // third zoom
  //center = VEC3F(-0.855658, 1.117847, 2.689409 - transformValue);
  //lengths *= 0.1 * 0.1;

  // third zoom, second attempt
  //center = VEC3F(-0.855658, 1.117847, 2.689409- transformValue);
  //lengths *= 0.1 * 0.25;

  // debugging 10 iterations on the centered bunny
  //center  = VEC3F(-1.21029,-0.763689,0);
  //lengths = VEC3F(0.593117,0.593117,0.593117);
  center  = VEC3F(-0.066911,-0.867573,-0.160586); 
  lengths = VEC3F(0.247994,0.247994,0.247994);

  //center = lengths;
  //center *= 0.25;
  //center[0] *= -1;
  //lengths *= 0.5;
  cout << " center:" << center << " lengths: " << lengths << endl;

  cout << " Allocating initial fields ... " << flush;
  optimize3D.fractal() = FIELD_3D(res,res,res, center, lengths);
  //optimize3D.auxiliary() = FIELD_3D(res,res,res, center, lengths);
  cout << " done. " << endl;

  /*
  if (transform.compare("slices") == 0)
  {
    cout << " Rendering quaternion slices " << endl;
    marchSlices();
  }
  else if (transform.compare("translationx") == 0)
  {
    cout << " Rendering root X translations " << endl;
    marchTranslations(VEC3F(1,0,0), "translationx");
  }
  else if (transform.compare("translationy") == 0)
  {
    cout << " Rendering root Y translations " << endl;
    marchTranslations(VEC3F(0,1,0), "translationy");
  }
  else
  if (transform.compare("translationz") == 0)
 */ 
  {
    cout << " Rendering root Z translations " << endl;
    marchTranslations(VEC3F(0,0,1), "translationz");
  }
  /*
  else if (transform.compare("rotationx") == 0)
  {
    cout << " Rendering root x rotations" << endl;
    marchRotations(VEC3F(1,0,0), "rotationx");
  }
  else if (transform.compare("rotationy") == 0)
  {
    cout << " Rendering root y rotations" << endl;
    marchRotations(VEC3F(0,1,0), "rotationy");
  }
  else if (transform.compare("rotationz") == 0)
  {
    cout << " Rendering root z rotations" << endl;
    marchRotations(VEC3F(0,0,1), "rotationz");
  }
  else if (transform.compare("scalex") == 0)
  {
    cout << " Rendering root x scalings" << endl;
    marchScaling(VEC3F(1,0,0), "scalex");
  }
  else if (transform.compare("scaley") == 0)
  {
    cout << " Rendering root y scalings" << endl;
    marchScaling(VEC3F(0,1,0), "scaley");
  }
  else if (transform.compare("scalez") == 0)
  {
    cout << " Rendering root z scalings" << endl;
    marchScaling(VEC3F(0,0,1), "scalez");
  }
  else if (transform.compare("scalexyz") == 0)
  {
    cout << " Rendering root xyz scalings" << endl;
    marchScaling(VEC3F(1,1,1), "scalexyz");
  }
  else if (transform.compare("translateZeroX") == 0)
  {
    cout << " Translating zero root along x" << endl;
    marchZeroRoot(VEC3F(1,0,0), "translateZeroX");
  }
  else if (transform.compare("conformal") == 0)
  {
    cout << " Trying different conformal radii " << endl;
    marchConformal("conformal");
  }
  else if (transform.compare("isosurface") == 0)
  {
    cout << " Rendering different isosurfaces" << endl;
    marchIsosurfaces("isosurface");
  }
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " UNKNOWN TRANSFORMATION " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  }
  */

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    cout << " USAGE: " << argv[0] << " <config file>" << endl;
    return 0;
  }

  string cfg(argv[1]);
  SIMPLE_PARSER parser(cfg);

  /*
  // which transformation are we doing?
  string transform(argv[5]);
  transformValue = atof(argv[6]); 
  cout << " Using transform type: " << transform.c_str() << endl;
  //cout << " Min: " << minValue << " Max: " << maxValue << " Steps: " << steps << endl;
  cout << " Transform value: " << transformValue << endl;
  */

  string filename       = parser.getString("optimization file", string(""), true);
  VEC3F rootTranslation = parser.getVector3("root translation", VEC3(0,0,0));
  outputPrefix     = parser.getString("output prefix", string("./temp/temp"));
  res              = parser.getInt("field res", 50);
  maxIterations    = parser.getInt("max iterations", 3);

  cout << " Using input file:       " << filename.c_str() << endl;
  cout << " Using root translation: " << rootTranslation << endl;
  cout << " Using output prefix:    " << outputPrefix.c_str() << endl;
  cout << " Using resolution:       " << res << endl;
  cout << " Using max iterations:   " << maxIterations << endl;

  optimize3D.read(filename.c_str());
  VEC3F& center = optimize3D.fractal().center();
  VEC3F& lengths = optimize3D.fractal().lengths();

  VEC3 centerNew  = parser.getVector3("field center", center);
  VEC3 lengthsNew = parser.getVector3("field lengths", lengths);

  cout << " Original center: " << center << endl;
  cout << " Original lengths: " << lengths << endl;
  cout << " New center: " << centerNew << endl;
  cout << " New lengths: " << lengthsNew << endl;

  // need to translate the center by the root translation
  //centerNew -= rootTranslation;

  //center = centerNew;
  lengths = lengthsNew;

#if !USING_LOWMEMORY
  cout << " Allocating initial fields ... " << flush;
  optimize3D.fractal() = FIELD_3D(res,res,res, center, lengths);
  cout << " done. " << endl;
#endif

  cout << " Rendering root Z translations " << endl;
  marchTranslations(rootTranslation, "translationz");

  return 0;
}
