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
 
  if (triangleMesh) delete triangleMesh;
  string cachePath;
  char buffer[256];
  sprintf(buffer, ".res.%i.iterations.%i.cache", res, maxIterations);
  string cacheFilename = outputPrefix + string(buffer);
  triangleMesh = new TRIANGLE_MESH(optimize3D.fractal().center(),
                                   optimize3D.fractal().lengths(),
                                   VEC3I(res,res,res),
                                   optimize3D.top(), 
                                   optimize3D.bottom(), 
                                   optimize3D.expScaling(), 
                                   optimize3D.maxIterations(), 
                                   optimize3D.slice(), 
                                   isosurface,
#if 1
                                   optimize3D.fractal().rotation());
#else
                                   optimize3D.fractal().rotation(), 
                                   cacheFilename);
#endif
  cout << " Marching cubes found " << triangleMesh->triangles().size() << " triangles " << endl;

  TIMER::printTimings();
}

void marchTranslations(const VEC3F& translation, const string& descriptor)
{
  float delta = transformValue;
  VEC3F      vDelta = translation;

  QUATERNION qDelta(vDelta[0], vDelta[1], vDelta[2], 0);
  cout << " Using the translation delta: " << delta <<  " quaternion version: " << qDelta << endl;

  const VEC3F centerOriginal = optimize3D.fractal().center();

  // back up the original root positions
  const POLYNOMIAL_4D topOriginal = optimize3D.top();
  const POLYNOMIAL_4D bottomOriginal = optimize3D.bottom();

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

  // write out the OBJ
  sprintf(buffer, ".res.%i.iterations.%i_%s.obj", res, maxIterations, descriptor.c_str());
  string outputFilename = outputPrefix + string(buffer);
  triangleMesh->writeOBJ(outputFilename.c_str());
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
  QUATERNION& rotation = optimize3D.fractal().rotation();

  VEC3 centerNew  = parser.getVector3("field center", center);
  VEC3 lengthsNew = parser.getVector3("field lengths", lengths);
  QUATERNION rotationNew = parser.getQuaternion("rotation", rotation);

  cout << " Original center: " << center << endl;
  cout << " Original lengths: " << lengths << endl;
  cout << " Original rotation: " << rotation << endl;
  cout << " New center: " << centerNew << endl;
  cout << " New lengths: " << lengthsNew << endl;
  cout << " New rotation: " << rotationNew << endl;

  bool recenter = parser.getBool("recenter", false);

  // need to translate the center by the root translation
  if (recenter)
  {
    cout << " Recentering " << centerNew << " to ";
    centerNew -= rootTranslation;
    cout << centerNew << endl;
  }
  else
    cout << " Not recentering" << endl;

  center = centerNew;
  lengths = lengthsNew;
  rotation = rotationNew;

  cout << " Allocating initial fields ... " << flush;
  optimize3D.fractal() = FIELD_3D(100,100,100, center, lengths);
  optimize3D.fractal().rotation() = rotationNew;
  cout << " done. " << endl;

  cout << " Rendering root Z translations " << endl;
  marchTranslations(rootTranslation, "translationz");

  return 0;
}
