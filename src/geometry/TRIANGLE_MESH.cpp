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
#include <algorithm>
#include <fstream>
//#include <omp.h>
#include "TIMER.h"

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH() :
  _glTextureHandle(0)
{
}

//////////////////////////////////////////////////////////////////////
// OBJ mesh constructor
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH(const string& filename) :
  _glTextureHandle(0)
{
  bool success = readOBJ(filename);

  if (!success)
    return;
}

//////////////////////////////////////////////////////////////////////
// Marching cubes constructor
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH(const FIELD_3D& field) :
  _glTextureHandle(0),
  _res(field.res()),
  _lengths(field.lengths()), 
  _center(field.center()), 
  _dxs(field.dx(), field.dy(), field.dz())
{
  computeStagedMarchingCubes(field);
}

//////////////////////////////////////////////////////////////////////
// Marching cubes constructor
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH(const FIELD_3D& field, const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const Real expScaling, const int maxIterations, const Real slice, const Real isosurface) :
  _glTextureHandle(0),
  _res(field.res()),
  _lengths(field.lengths()), 
  _center(field.center()), 
  _dxs(field.dx(), field.dy(), field.dz())
{
  _top = top;
  _bottom = bottom;
  _expScaling = expScaling;
  _maxIterations = maxIterations;
  _quaternionSlice = slice;
  _isosurface = isosurface;
  _escapeRadius = 2000.0;

  computeNonlinearMarchingCubes(field);
}

//////////////////////////////////////////////////////////////////////
// do a non-linear marching cubes on just two slabs at a shot
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH(const VEC3F& center, const VEC3F& lengths, const VEC3I& res, const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const Real expScaling, const int maxIterations, const Real slice, const Real isosurface, const string& cacheFilename) :
  _res(res),
  _lengths(lengths), 
  _center(center), 
  _dxs(lengths[0] / res[0], lengths[1] / res[1], lengths[2] / res[2])
{
  _top = top;
  _bottom = bottom;
  _expScaling = expScaling;
  _maxIterations = maxIterations;
  _quaternionSlice = slice;
  _isosurface = isosurface;
  _escapeRadius = 2000.0;
  _cacheFilename = cacheFilename;

  //computeNonlinearMarchingCubesLowMemory();
  computeNonlinearMarchingCubesLowMemoryHuge();
}

//////////////////////////////////////////////////////////////////////
// Marching cubes constructor, but just for quadratic Julia sets
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH(const FIELD_3D& field, const QUATERNION& qConst, const int maxIterations) :
  _glTextureHandle(0),
  _quadraticConst(qConst),
  _res(field.res()),
  _lengths(field.lengths()), 
  _center(field.center()), 
  _dxs(field.dx(), field.dy(), field.dz())
{
  _maxIterations = maxIterations;

  computeQuadraticMarchingCubes(field);
}

//////////////////////////////////////////////////////////////////////
// destructor
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::~TRIANGLE_MESH()
{
}

//////////////////////////////////////////////////////////////////////
// perform marching cubes
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearMarchingCubes(const FIELD_3D& field, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
    cout << " Marching cubes ...";flush(cout);

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexHash.clear();

  _fieldDeltas[0] = field.dx();
  _fieldDeltas[1] = field.dy();
  _fieldDeltas[2] = field.dz();

  // set "outside" to something a lot bigger than the known grid bounds
  _outside = field.center()[0] + field.lengths()[0] * 10000;
 
  _xRes = field.xRes();
  _yRes = field.yRes();
  _zRes = field.zRes();
  _slabSize = _xRes * _yRes;

  _toMarch = &field;

  // build all the vertex pairs 
  _vertexPairs.clear();
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();

        CUBE cube;
        cube.NNN = field(x,y,z);
        cube.NNP = field(x,y,z + 1);
        cube.NPN = field(x,y + 1,z);
        cube.NPP = field(x,y + 1,z + 1);
        cube.PNN = field(x + 1,y,z);
        cube.PNP = field(x + 1,y,z+1);
        cube.PPN = field(x + 1,y + 1,z);
        cube.PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));
		
        // three vertices are added to _vertexPairs here  
        switch (flag)
#include "MARCHING_CUBES_VERTICES.include" 
      }
    //if (verbose)
    //  cout << z << "/" << field.zRes() << " "; flush(cout);
  }

  // compute the interpolations along the marching cubes edges
  computeNonlinearEdgeInterpolations();
 
  // build the actual triangles
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();

        CUBE cube;
        cube.NNN = field(x,y,z);
        cube.NNP = field(x,y,z + 1);
        cube.NPN = field(x,y + 1,z);
        cube.NPP = field(x,y + 1,z + 1);
        cube.PNN = field(x + 1,y,z);
        cube.PNP = field(x + 1,y,z+1);
        cube.PPN = field(x + 1,y + 1,z);
        cube.PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));
		  
        switch (flag)
#include "MARCHING_CUBES_TRIANGLES.include" 
      }
    //if (verbose)
    //  cout << z << "/" << field.zRes() << " "; flush(cout);
  }

  // create the final triangles based on the vertex indices -- this
  // couldn't be done in the inner loop because the vector keeps
  // resizing and changing the vertex addresses
  assert(_triangleVertices.size() % 3 == 0);
	
  if (verbose) cout << "computed triangles: " << _triangleVertices.size() << endl;
	
  for (unsigned int x = 0; x < _triangleVertices.size() / 3; x++)
  {
    VEC3F* v0 = &_vertices[_triangleVertices[3 * x]];
    VEC3F* v1 = &_vertices[_triangleVertices[3 * x + 1]];
    VEC3F* v2 = &_vertices[_triangleVertices[3 * x + 2]];

    // ignore any degenerate triangles
    Real dist0 = norm((*v0) - (*v1));
    Real dist1 = norm((*v0) - (*v2));
    Real dist2 = norm((*v1) - (*v2));
    Real eps = 1e-7;
    if (dist0 < eps || dist1 < eps || dist2 < eps)
      continue;
    _triangles.push_back(TRIANGLE(v0, v1, v2));
  }

  // all done -- throw away the indices
  _triangleVertices.clear();

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  if (verbose)
    cout << "done. " << endl;
}


//////////////////////////////////////////////////////////////////////
// support function for computeNonlinearMarchingCubesLowMemory(), 
// which computes a single slice of the potential function
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearSlice(const int z, FIELD_2D& field)
{
  int xRes = _res[0];
  int yRes = _res[1];

  if (field.xRes() != xRes || field.yRes() != yRes)
    field.resizeAndWipe(xRes, yRes, _center, _lengths);

  //Real escape = 20.0;
  Real escape = _escapeRadius;
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      VEC3F point = cellCenter(x,y,z);
      Real xReal = point[0];
      Real yReal = point[1];
      Real zReal = point[2];

      QUATERNION iterate(xReal, yReal, zReal, _quaternionSlice);
   
      Real magnitude = iterate.magnitude();
      int totalIterations = 0;
      //Real magnitudeTrial = 0;
      while (magnitude < escape && totalIterations < _maxIterations) 
      {
        QUATERNION topEval = _top.evaluateScaledPowerFactored(iterate);
        QUATERNION bottomEval;
        
        if (_bottom.totalRoots() > 0)
        {
          bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
          iterate = (topEval / bottomEval);
        }
        else
          iterate = topEval;

        //magnitudeTrial = iterate.magnitude();

        iterate *= _expScaling;
        magnitude = iterate.magnitude();

        totalIterations++;

        // see if it fell into the black hole at the origin
        if (magnitude < 10.0 * REAL_MIN)
          totalIterations = _maxIterations;
      }
      field(x,y) = log(magnitude);
    }
}

//////////////////////////////////////////////////////////////////////
// try to read in an edge cache
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readEdgeCache(vector<pair<int, int> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No cache file found: " << _cacheFilename.c_str() << " " << endl;
    return false;
  }
  cout << " Cache file found: " << _cacheFilename.c_str() << endl;
 
  // read in: map<pair<int, int>, bool> _vertexPairs;
  int totalPairs;
  fread((void*)&(totalPairs), sizeof(int), 1, file);

  _vertexPairs.clear();
  for (int x = 0; x < totalPairs; x++)
  {
    pair<int,int> vertexPair;
    fread((void*)&(vertexPair.first), sizeof(int), 1, file);
    fread((void*)&(vertexPair.second), sizeof(int), 1, file);
    _vertexPairs[vertexPair] = true;
  }

  // read in: vector<pair<int, int> > flags
  int totalFlags;
  fread((void*)&(totalFlags), sizeof(int), 1, file);

  flags.clear();
  for (int x = 0; x < totalFlags; x++)
  {
    pair<int, int> flagPair;
    fread((void*)&(flagPair.first), sizeof(int), 1, file);
    fread((void*)&(flagPair.second), sizeof(int), 1, file);
    flags.push_back(flagPair);
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// try to read in an edge cache
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readEdgeCacheHuge(vector<pair<int, VEC3I> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No cache file found: " << _cacheFilename.c_str() << " " << endl;
    return false;
  }
  cout << " Cache file found: " << _cacheFilename.c_str() << endl;
 
  // read in: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets;
  fread((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " Reading in " << totalTriplets << " triplets " << endl;

  _vertexTriplets.clear();
  for (int x = 0; x < totalTriplets; x++)
  {
    pair<VEC3I, VEC3I> vertexTriplet;
    VEC3I first;
    VEC3I second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(first[y]), sizeof(int), 1, file);
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);

    vertexTriplet.first = first;
    vertexTriplet.second = second;
    _vertexTriplets[vertexTriplet] = true;
  }
  cout << " Stored " << _vertexTriplets.size() << " triplets" << endl;

  // read in: vector<pair<int, int> > flags
  int totalFlags;
  fread((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " Reading in " << totalFlags << " flags " << endl;

  flags.clear();
  for (int x = 0; x < totalFlags; x++)
  {
    pair<int, VEC3I> flagPair;
    fread((void*)&(flagPair.first), sizeof(int), 1, file);
    VEC3I& second = flagPair.second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);
    flags.push_back(flagPair);
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// write out an edge cache
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeEdgeCache(const vector<pair<int, int> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FAILED TO CACHE OUT EDGES: " << _cacheFilename.c_str() << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  cout << " Writing out cache file " << _cacheFilename.c_str() << " ... " << flush;
 
  // write out: map<pair<int, int>, bool> _vertexPairs;
  int totalPairs = _vertexPairs.size(); 
  fwrite((void*)&(totalPairs), sizeof(int), 1, file);

  map<pair<int, int>, bool>::iterator iter;
  for (iter = _vertexPairs.begin(); iter != _vertexPairs.end(); iter++)
  {
    pair<int,int> vertexPair = iter->first;
    fwrite((void*)&(vertexPair.first), sizeof(int), 1, file);
    fwrite((void*)&(vertexPair.second), sizeof(int), 1, file);
  }
  
  // write out: vector<pair<int, int> > flags
  int totalFlags = flags.size(); 
  fwrite((void*)&(totalFlags), sizeof(int), 1, file);
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    pair<int, int> flagPair = flags[x];
    fwrite((void*)&(flagPair.first), sizeof(int), 1, file);
    fwrite((void*)&(flagPair.second), sizeof(int), 1, file);
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// write out an edge cache
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeEdgeCacheHuge(const vector<pair<int, VEC3I> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FAILED TO CACHE OUT EDGES: " << _cacheFilename.c_str() << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  cout << " Writing out cache file " << _cacheFilename.c_str() << " ... " << flush;
 
  // write out: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets = _vertexTriplets.size(); 
  fwrite((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " total triplets: " << totalTriplets << " ";

  map<pair<VEC3I, VEC3I>, bool>::iterator iter;
  for (iter = _vertexTriplets.begin(); iter != _vertexTriplets.end(); iter++)
  {
    pair<VEC3I,VEC3I> vertexTriplet = iter->first;

    VEC3I& first  = vertexTriplet.first;
    VEC3I& second = vertexTriplet.second;
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(first[x]), sizeof(int), 1, file);
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(second[x]), sizeof(int), 1, file);
  }
  
  // write out: vector<pair<int, int> > flags
  int totalFlags = flags.size(); 
  fwrite((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " total flags: " << totalFlags << " ";
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    pair<int, VEC3I> flagPair = flags[x];
    fwrite((void*)&(flagPair.first), sizeof(int), 1, file);
    VEC3I second = flagPair.second;
    for (int y = 0; y < 3; y++)
      fwrite((void*)&(second[y]), sizeof(int), 1, file);
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// perform marching cubes
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearMarchingCubesLowMemory()
{
  bool verbose = true;
  TIMER functionTimer(__FUNCTION__);

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexHash.clear();

  // set "outside" to something a lot bigger than the known grid bounds
  //_outside = field.center()[0] + field.lengths()[0] * 10000;
  _outside = _center[0] + _lengths[0] * 10000;
 
  _xRes = _res[0];
  _yRes = _res[1];
  _zRes = _res[2];
  _slabSize = _xRes * _yRes;

  // store all the meaningful flags <flag, index>
  vector<pair<int, int> > flags;

  // see if a previous run was cached before computing all the slices
  if (readEdgeCache(flags) == false)
  {
    if (verbose)
      cout << " Low-memory marching cubes ..." << flush;
    computeAllLowMemorySlices(flags);
    //writeEdgeCache(flags);
  }
  TIMER::printTimings();
  cout << " Found " << flags.size() << " flags. " << endl;
  cout << " Found " << _vertexPairs.size() << " vertex pairs. " << endl;

  // compute the interpolations along the marching cubes edges
  computeNonlinearEdgeInterpolations();

  // go back over the vertex pairs and emit the triangles
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    int flag = flags[x].first;
    int index = flags[x].second;
  
    switch (flag)
#include "MARCHING_CUBES_TRIANGLES.include"
  }

  // create the final triangles based on the vertex indices -- this
  // couldn't be done in the inner loop because the vector keeps
  // resizing and changing the vertex addresses
  assert(_triangleVertices.size() % 3 == 0);

  if (verbose) cout << "computed triangles: " << _triangleVertices.size() << endl;
	
  int totalDegenerate = 0;  
  for (unsigned int x = 0; x < _triangleVertices.size() / 3; x++)
  {
    VEC3F* v0 = &_vertices[_triangleVertices[3 * x]];
    VEC3F* v1 = &_vertices[_triangleVertices[3 * x + 1]];
    VEC3F* v2 = &_vertices[_triangleVertices[3 * x + 2]];

    // ignore any degenerate triangles
    Real dist0 = norm((*v0) - (*v1));
    Real dist1 = norm((*v0) - (*v2));
    Real dist2 = norm((*v1) - (*v2));
    Real eps = 1e-7;
    if (dist0 < eps || dist1 < eps || dist2 < eps)
    {
      totalDegenerate++;
      continue;
    }
    _triangles.push_back(TRIANGLE(v0, v1, v2));
  }
  cout << " Found " << totalDegenerate << " degenerate triangles " << endl;

  // all done -- throw away the indices
  _triangleVertices.clear();

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  if (verbose)
    cout << "done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// perform marching cubes
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearMarchingCubesLowMemoryHuge()
{
  bool verbose = true;
  TIMER functionTimer(__FUNCTION__);

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexTripletHash.clear();

  // set "outside" to something a lot bigger than the known grid bounds
  //_outside = field.center()[0] + field.lengths()[0] * 10000;
  _outside = _center[0] + _lengths[0] * 10000;
 
  _xRes = _res[0];
  _yRes = _res[1];
  _zRes = _res[2];
  _slabSize = _xRes * _yRes;

  // store all the meaningful flags <flag, index>
  vector<pair<int, VEC3I> > flags;

  // see if a previous run was cached before computing all the slices
  if (readEdgeCacheHuge(flags) == false)
  {
    if (verbose)
      cout << " Low-memory marching cubes ..." << flush;
    computeAllLowMemorySlicesHuge(flags);
    writeEdgeCacheHuge(flags);
  }
  //computeAllLowMemorySlicesHuge(flags);
  TIMER::printTimings();
  cout << " Found " << flags.size() << " flags. " << endl;
  cout << " Found " << _vertexTriplets.size() << " vertex triplets. " << endl;

  // compute the interpolations along the marching cubes edges
  computeNonlinearEdgeInterpolationsHuge();

  // go back over the vertex pairs and emit the triangles
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    int flag = flags[x].first;
    VEC3I index = flags[x].second;
  
    switch (flag)
//#include "MARCHING_CUBES_TRIANGLES.include"
#include "MARCHING_CUBES_TRIANGLES.include.huge"
  }

  // create the final triangles based on the vertex indices -- this
  // couldn't be done in the inner loop because the vector keeps
  // resizing and changing the vertex addresses
  assert(_triangleVertices.size() % 3 == 0);

  if (verbose) cout << "computed triangles: " << _triangleVertices.size() << endl;

  int totalDegenerate = 0;  
  for (unsigned int x = 0; x < _triangleVertices.size() / 3; x++)
  {
    VEC3F* v0 = &_vertices[_triangleVertices[3 * x]];
    VEC3F* v1 = &_vertices[_triangleVertices[3 * x + 1]];
    VEC3F* v2 = &_vertices[_triangleVertices[3 * x + 2]];

    // ignore any degenerate triangles
    Real dist0 = norm((*v0) - (*v1));
    Real dist1 = norm((*v0) - (*v2));
    Real dist2 = norm((*v1) - (*v2));
    Real eps = 1e-7;
    if (dist0 < eps || dist1 < eps || dist2 < eps)
    {
      totalDegenerate++;
      continue;
    }
    _triangles.push_back(TRIANGLE(v0, v1, v2));
  }
  cout << " Found " << totalDegenerate << " degenerate triangles " << endl;

  // all done -- throw away the indices
  _triangleVertices.clear();

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  if (verbose)
    cout << "done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// perform marching cubes
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeQuadraticMarchingCubes(const FIELD_3D& field, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
    cout << " Marching cubes ...";flush(cout);

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexHash.clear();

  _fieldDeltas[0] = field.dx();
  _fieldDeltas[1] = field.dy();
  _fieldDeltas[2] = field.dz();

  // set "outside" to something a lot bigger than the known grid bounds
  _outside = field.center()[0] + field.lengths()[0] * 10000;
 
  _xRes = field.xRes();
  _yRes = field.yRes();
  _zRes = field.zRes();
  _slabSize = _xRes * _yRes;

  _toMarch = &field;

  // build all the vertex pairs 
  _vertexPairs.clear();
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();

        //_cellCenter = field.cellCenter(x,y,z);
        //const VEC3F center = field.cellCenter(x,y,z);

        CUBE cube;
        cube.NNN = field(x,y,z);
        cube.NNP = field(x,y,z + 1);
        cube.NPN = field(x,y + 1,z);
        cube.NPP = field(x,y + 1,z + 1);
        cube.PNN = field(x + 1,y,z);
        cube.PNP = field(x + 1,y,z+1);
        cube.PPN = field(x + 1,y + 1,z);
        cube.PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));
		  
        switch (flag)
#include "MARCHING_CUBES_VERTICES.include" 
      }
    //if (verbose)
    //  cout << z << "/" << field.zRes() << " "; flush(cout);
  }

  // compute the interpolations along the marching cubes edges
  computeQuadraticEdgeInterpolations();
 
  // build the actual triangles
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();

        CUBE cube;
        cube.NNN = field(x,y,z);
        cube.NNP = field(x,y,z + 1);
        cube.NPN = field(x,y + 1,z);
        cube.NPP = field(x,y + 1,z + 1);
        cube.PNN = field(x + 1,y,z);
        cube.PNP = field(x + 1,y,z+1);
        cube.PPN = field(x + 1,y + 1,z);
        cube.PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));
		  
        switch (flag)
#include "MARCHING_CUBES_TRIANGLES.include" 
      }
    //if (verbose)
    //  cout << z << "/" << field.zRes() << " "; flush(cout);
  }

  // create the final triangles based on the vertex indices -- this
  // couldn't be done in the inner loop because the vector keeps
  // resizing and changing the vertex addresses
  assert(_triangleVertices.size() % 3 == 0);
	
  if (verbose) cout << "computed triangles: " << _triangleVertices.size() << endl;
	
  for (unsigned int x = 0; x < _triangleVertices.size() / 3; x++)
  {
    VEC3F* v0 = &_vertices[_triangleVertices[3 * x]];
    VEC3F* v1 = &_vertices[_triangleVertices[3 * x + 1]];
    VEC3F* v2 = &_vertices[_triangleVertices[3 * x + 2]];

    // ignore any degenerate triangles
    Real dist0 = norm((*v0) - (*v1));
    Real dist1 = norm((*v0) - (*v2));
    Real dist2 = norm((*v1) - (*v2));
    Real eps = 1e-7;
    if (dist0 < eps || dist1 < eps || dist2 < eps)
      continue;
    _triangles.push_back(TRIANGLE(v0, v1, v2));
  }

  // all done -- throw away the indices
  _triangleVertices.clear();

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  if (verbose)
    cout << "done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// perform marching cubes
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeStagedMarchingCubes(const FIELD_3D& field, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
    cout << " Marching cubes ...";flush(cout);

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexHash.clear();

  _fieldDeltas[0] = field.dx();
  _fieldDeltas[1] = field.dy();
  _fieldDeltas[2] = field.dz();

  // set "outside" to something a lot bigger than the known grid bounds
  _outside = field.center()[0] + field.lengths()[0] * 10000;
 
  _xRes = field.xRes();
  _yRes = field.yRes();
  _zRes = field.zRes();
  _slabSize = _xRes * _yRes;

  _toMarch = &field;

  // build all the vertex pairs 
  _vertexPairs.clear();
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();

        //_cellCenter = field.cellCenter(x,y,z);
        //const VEC3F center = field.cellCenter(x,y,z);

        CUBE cube;
        cube.NNN = field(x,y,z);
        cube.NNP = field(x,y,z + 1);
        cube.NPN = field(x,y + 1,z);
        cube.NPP = field(x,y + 1,z + 1);
        cube.PNN = field(x + 1,y,z);
        cube.PNP = field(x + 1,y,z+1);
        cube.PPN = field(x + 1,y + 1,z);
        cube.PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));
		  
        switch (flag)
#include "MARCHING_CUBES_VERTICES.include" 
      }
    //if (verbose)
    //  cout << z << "/" << field.zRes() << " "; flush(cout);
  }

  // compute the interpolations along the marching cubes edges
  computeEdgeInterpolations();
 
  // build the actual triangles
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();

        CUBE cube;
        cube.NNN = field(x,y,z);
        cube.NNP = field(x,y,z + 1);
        cube.NPN = field(x,y + 1,z);
        cube.NPP = field(x,y + 1,z + 1);
        cube.PNN = field(x + 1,y,z);
        cube.PNP = field(x + 1,y,z+1);
        cube.PPN = field(x + 1,y + 1,z);
        cube.PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));
		  
        switch (flag)
#include "MARCHING_CUBES_TRIANGLES.include" 
      }
    //if (verbose)
    //  cout << z << "/" << field.zRes() << " "; flush(cout);
  }

  // create the final triangles based on the vertex indices -- this
  // couldn't be done in the inner loop because the vector keeps
  // resizing and changing the vertex addresses
  assert(_triangleVertices.size() % 3 == 0);
	
  if (verbose) cout << "computed triangles: " << _triangleVertices.size() << endl;
	
  for (unsigned int x = 0; x < _triangleVertices.size() / 3; x++)
  {
    VEC3F* v0 = &_vertices[_triangleVertices[3 * x]];
    VEC3F* v1 = &_vertices[_triangleVertices[3 * x + 1]];
    VEC3F* v2 = &_vertices[_triangleVertices[3 * x + 2]];

    // ignore any degenerate triangles
    Real dist0 = norm((*v0) - (*v1));
    Real dist1 = norm((*v0) - (*v2));
    Real dist2 = norm((*v1) - (*v2));
    Real eps = 1e-7;
    if (dist0 < eps || dist1 < eps || dist2 < eps)
      continue;
    _triangles.push_back(TRIANGLE(v0, v1, v2));
  }

  // all done -- throw away the indices
  _triangleVertices.clear();

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  if (verbose)
    cout << "done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// perform marching cubes
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeMarchingCubes(const FIELD_3D& field, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
    cout << " Marching cubes ...";flush(cout);

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexHash.clear();

  _fieldDeltas[0] = field.dx();
  _fieldDeltas[1] = field.dy();
  _fieldDeltas[2] = field.dz();

  // set "outside" to something a lot bigger than the known grid bounds
  _outside = field.center()[0] + field.lengths()[0] * 10000;
 
  _xRes = field.xRes();
  _yRes = field.yRes();
  _zRes = field.zRes();
  _slabSize = _xRes * _yRes;
  
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();
        _cellCenter = field.cellCenter(x,y,z);

        _NNN = field(x,y,z);
        _NNP = field(x,y,z + 1);
        _NPN = field(x,y + 1,z);
        _NPP = field(x,y + 1,z + 1);
        _PNN = field(x + 1,y,z);
        _PNP = field(x + 1,y,z+1);
        _PPN = field(x + 1,y + 1,z);
        _PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((_NNN > 0) + 2 *   (_NNP > 0) + 4  * (_NPN > 0) +
                   8 * (_NPP > 0) + 16 *  (_PNN > 0) + 32 * (_PNP > 0) +
                   64 *(_PPN > 0) + 128 * (_PPP > 0));
		  
//		if (verbose && flag != 0 && flag != 255) {
//			cout << "x: " << x << ", y: " << y << ", z: " << z;
//			cout << ", marching flag:" << flag;
//			cout << " _NNN: " <<  _NNN <<  "_NNP: " <<  _NNP;
//			cout << endl;
//		}
		  
        switch (flag) 
        {
        case 0:  case 255: break;
        case 1: addTriangle(2,1,10,index); break;
        case 2: addTriangle(2,11,3,index); break;
        case 3: addTriangle(10,11,1,index); addTriangle(11,3,1,index); break;
        case 4: addTriangle(1,4,9,index); break;
        case 5: addTriangle(10,2,4,index); addTriangle(4,9,10,index); break;
        case 6: addTriangle(1,4,9,index); addTriangle(2,11,3,index); break;
        case 7: addTriangle(10,4,9,index); addTriangle(10,3,4,index); addTriangle(10,11,3,index); break;
        case 8: addTriangle(3,12,4,index); break;
        case 9: addTriangle(3,12,4,index); addTriangle(2,1,10,index); break;
        case 10: addTriangle(11,12,4,index); addTriangle(4,2,11,index); break;
        case 11: addTriangle(11,1,10,index); addTriangle(11,4,1,index); addTriangle(11,12,4,index); break;
        case 12: addTriangle(9,1,3,index); addTriangle(3,12,9,index); break;
        case 13: addTriangle(9,3,12,index); addTriangle(9,2,3,index); addTriangle(9,10,2,index); break;
        case 14: addTriangle(12,9,1,index); addTriangle(12,1,2,index); addTriangle(12,2,11,index); break;
        case 15: addTriangle(11,12,9,index); addTriangle(9,10,11,index); break;
        case 16: addTriangle(10,5,6,index); break; 
        case 17: addTriangle(5,6,2,index); addTriangle(2,1,5,index); break;
        case 18: addTriangle(10,5,6,index); addTriangle(2,11,3,index); break;
        case 19: addTriangle(1,5,6,index); addTriangle(1,6,11,index); addTriangle(1,11,3,index); break;
        case 20: addTriangle(10,5,6,index); addTriangle(1,4,9,index); break;
        case 21: addTriangle(2,4,9,index); addTriangle(2,9,5,index); addTriangle(2,5,6,index); break;
        case 22: addTriangle(2,11,3,index); addTriangle(1,4,9,index); addTriangle(10,5,6,index); break;
        case 23: addTriangle(5,6,11,index); addTriangle(9,5,11,index); addTriangle(11,3,9,index); addTriangle(3,4,9,index); break;
        case 24: addTriangle(10,5,6,index); addTriangle(3,12,4,index); break;
        case 25: addTriangle(5,6,2,index); addTriangle(2,1,5,index); addTriangle(3,12,4,index); break;
        case 26: addTriangle(11,12,4,index); addTriangle(4,2,11,index); addTriangle(10,5,6,index); break;
        case 27: addTriangle(5,6,11,index); addTriangle(11,1,5,index); addTriangle(1,11,12,index); addTriangle(12,4,1,index); break;
        case 28: addTriangle(9,1,3,index); addTriangle(3,12,9,index); addTriangle(10,5,6,index); break;
        case 29: addTriangle(3,12,9,index); addTriangle(9,2,3,index); addTriangle(2,9,5,index); addTriangle(5,6,2,index); break;
        case 30: addTriangle(12,9,1,index); addTriangle(12,1,2,index); addTriangle(12,2,11,index); addTriangle(10,5,6,index); break;
        case 31: addTriangle(12,9,5,index); addTriangle(12,5,6,index); addTriangle(12,6,11,index); break;
        case 32: addTriangle(6,7,11,index); break;
        case 33: addTriangle(6,7,11,index); addTriangle(2,1,10,index); break;
        case 34: addTriangle(6,7,3,index); addTriangle(3,2,6,index); break;
        case 35: addTriangle(3,1,10,index); addTriangle(3,10,6,index); addTriangle(3,6,7,index); break;
        case 36: addTriangle(6,7,11,index); addTriangle(1,4,9,index); break;
        case 37: addTriangle(10,2,4,index); addTriangle(4,9,10,index); addTriangle(6,7,11,index); break;
        case 38: addTriangle(6,7,3,index); addTriangle(3,2,6,index); addTriangle(1,4,9,index); break;
        case 39: addTriangle(4,9,3,index); addTriangle(3,9,10,index); addTriangle(10,7,3,index); addTriangle(10,6,7,index); break;
        case 40: addTriangle(6,7,11,index); addTriangle(3,12,4,index); break;
        case 41: addTriangle(2,1,10,index); addTriangle(3,12,4,index); addTriangle(6,7,11,index); break;
        case 42: addTriangle(2,6,7,index); addTriangle(2,7,12,index); addTriangle(2,12,4,index); break;
        case 43: addTriangle(1,12,4,index); addTriangle(1,7,12,index); addTriangle(1,10,7,index); addTriangle(10,6,7,index); break;
        case 44: addTriangle(9,1,3,index); addTriangle(3,12,9,index); addTriangle(6,7,11,index); break;
        case 45: addTriangle(9,3,12,index); addTriangle(9,2,3,index); addTriangle(9,10,2,index); addTriangle(6,7,11,index); break;
        case 46: addTriangle(9,6,12,index); addTriangle(12,6,7,index); addTriangle(9,2,6,index); addTriangle(9,1,2,index); break;
        case 47: addTriangle(9,10,6,index); addTriangle(9,6,7,index); addTriangle(9,7,12,index); break;
        case 48: addTriangle(11,10,5,index); addTriangle(5,7,11,index); break; 
        case 49: addTriangle(5,7,11,index); addTriangle(5,11,2,index); addTriangle(5,2,1,index); break; 
        case 50: addTriangle(7,3,2,index); addTriangle(7,2,10,index); addTriangle(7,10,5,index); break; 
        case 51: addTriangle(5,7,3,index); addTriangle(3,1,5,index); break; 
        case 52: addTriangle(11,10,5,index); addTriangle(5,7,11,index); addTriangle(1,4,9,index); break; 
        case 53: addTriangle(11,2,5,index); addTriangle(11,5,7,index); addTriangle(2,4,9,index); addTriangle(9,5,2,index); break; 
        case 54: addTriangle(7,3,2,index); addTriangle(7,2,10,index); addTriangle(7,10,5,index); addTriangle(1,4,9,index); break; 
        case 55: addTriangle(7,3,4,index); addTriangle(7,4,9,index); addTriangle(7,9,5,index); break; 
        case 56: addTriangle(11,10,5,index); addTriangle(5,7,11,index); addTriangle(3,12,4,index); break; 
        case 57: addTriangle(5,7,11,index); addTriangle(5,11,2,index); addTriangle(5,2,1,index); addTriangle(3,12,4,index); break; 
        case 58: addTriangle(2,10,5,index); addTriangle(2,5,7,index); addTriangle(7,4,2,index); addTriangle(7,12,4,index); break; 
        case 59: addTriangle(5,7,12,index); addTriangle(5,12,4,index); addTriangle(5,4,1,index); break; 
        case 60: addTriangle(9,1,3,index); addTriangle(3,12,9,index); addTriangle(11,10,5,index); addTriangle(5,7,11,index); break; 
        case 61: addTriangle(2,9,5,index); addTriangle(9,2,3,index); addTriangle(3,12,9,index); addTriangle(5,7,11,index); addTriangle(11,2,5,index); break; 
        case 62: addTriangle(2,7,12,index); addTriangle(7,2,10,index); addTriangle(10,5,7,index); addTriangle(12,9,1,index); addTriangle(1,2,12,index); break; 
        case 63: addTriangle(12,9,5,index); addTriangle(5,7,12,index); break; 
        case 64: addTriangle(9,8,5,index); break; 
        case 65: addTriangle(9,8,5,index); addTriangle(2,1,10,index); break; 
        case 66: addTriangle(9,8,5,index); addTriangle(2,11,3,index); break; 
        case 67: addTriangle(10,11,1,index); addTriangle(11,3,1,index); addTriangle(9,8,5,index); break; 
        case 68: addTriangle(1,4,8,index); addTriangle(8,5,1,index); break; 
        case 69: addTriangle(4,8,5,index); addTriangle(4,5,10,index); addTriangle(4,10,2,index); break; 
        case 70: addTriangle(1,4,8,index); addTriangle(8,5,1,index); addTriangle(2,11,3,index); break; 
        case 71: addTriangle(3,4,8,index); addTriangle(8,11,3,index); addTriangle(5,10,11,index); addTriangle(11,8,5,index); break; 
        case 72: addTriangle(9,8,5,index); addTriangle(3,12,4,index); break; 
        case 73: addTriangle(2,1,10,index); addTriangle(3,12,4,index); addTriangle(9,8,5,index); break; 
        case 74: addTriangle(11,12,4,index); addTriangle(4,2,11,index); addTriangle(9,8,5,index); break; 
        case 75: addTriangle(11,1,10,index); addTriangle(11,4,1,index); addTriangle(11,12,4,index); addTriangle(9,8,5,index); break; 
        case 76: addTriangle(1,3,12,index); addTriangle(1,12,8,index); addTriangle(1,8,5,index); break; 
        case 77: addTriangle(5,10,2,index); addTriangle(5,2,3,index); addTriangle(5,3,8,index); addTriangle(8,3,12,index); break; 
        case 78: addTriangle(12,1,2,index); addTriangle(2,11,12,index); addTriangle(1,12,8,index); addTriangle(8,5,1,index); break; 
        case 79: addTriangle(11,12,8,index); addTriangle(11,8,5,index); addTriangle(11,5,10,index); break; 
        case 80: addTriangle(6,10,9,index); addTriangle(9,8,6,index); break; 
        case 81: addTriangle(6,2,1,index); addTriangle(6,1,9,index); addTriangle(6,9,8,index); break; 
        case 82: addTriangle(6,10,9,index); addTriangle(9,8,6,index); addTriangle(2,11,3,index); break; 
        case 83: addTriangle(8,6,11,index); addTriangle(11,3,8,index); addTriangle(3,1,9,index); addTriangle(9,8,3,index); break; 
        case 84: addTriangle(8,6,10,index); addTriangle(8,10,1,index); addTriangle(8,1,4,index); break; 
        case 85: addTriangle(2,4,8,index); addTriangle(8,6,2,index); break; 
        case 86: addTriangle(8,6,10,index); addTriangle(8,10,1,index); addTriangle(8,1,4,index); addTriangle(2,11,3,index); break; 
        case 87: addTriangle(8,6,11,index); addTriangle(8,11,3,index); addTriangle(8,3,4,index); break; 
        case 88: addTriangle(6,10,9,index); addTriangle(9,8,6,index); addTriangle(3,12,4,index); break; 
        case 89: addTriangle(6,2,1,index); addTriangle(6,1,9,index); addTriangle(6,9,8,index); addTriangle(3,12,4,index); break; 
        case 90: addTriangle(11,12,4,index); addTriangle(4,2,11,index); addTriangle(6,10,9,index); addTriangle(9,8,6,index); break; 
        case 91: addTriangle(1,6,11,index); addTriangle(11,12,4,index); addTriangle(4,1,11,index); addTriangle(6,1,9,index); addTriangle(9,8,6,index); break; 
        case 92: addTriangle(3,6,1,index); addTriangle(6,3,8,index); addTriangle(3,12,8,index); addTriangle(6,10,1,index); break; 
        case 93: addTriangle(6,2,3,index); addTriangle(6,3,12,index); addTriangle(6,12,8,index); break; 
        case 94: addTriangle(2,11,12,index); addTriangle(2,12,1,index); addTriangle(1,12,8,index); addTriangle(10,1,8,index); addTriangle(10,8,6,index); break; 
        case 95: addTriangle(8,6,11,index); addTriangle(11,12,8,index); break; 
        case 96: addTriangle(9,8,5,index); addTriangle(6,7,11,index); break; 
        case 97: addTriangle(2,1,10,index); addTriangle(6,7,11,index); addTriangle(9,8,5,index); break; 
        case 98: addTriangle(6,7,3,index); addTriangle(3,2,6,index); addTriangle(9,8,5,index); break; 
        case 99: addTriangle(3,1,10,index); addTriangle(3,10,6,index); addTriangle(3,6,7,index); addTriangle(9,8,5,index); break; 
        case 100: addTriangle(1,4,8,index); addTriangle(8,5,1,index); addTriangle(6,7,11,index); break; 
        case 101: addTriangle(4,8,5,index); addTriangle(4,5,10,index); addTriangle(4,10,2,index); addTriangle(6,7,11,index); break; 
        case 102: addTriangle(6,7,3,index); addTriangle(3,2,6,index); addTriangle(1,4,8,index); addTriangle(8,5,1,index); break; 
        case 103: addTriangle(6,7,3,index); addTriangle(10,6,3,index); addTriangle(5,4,8,index); addTriangle(10,4,5,index); addTriangle(10,3,4,index); break; 
        case 104: addTriangle(9,8,5,index); addTriangle(6,7,11,index); addTriangle(3,12,4,index); break; 
        case 105: addTriangle(2,1,10,index); addTriangle(3,12,4,index); addTriangle(6,7,11,index); addTriangle(9,8,5,index); break; 
        case 106: addTriangle(2,6,7,index); addTriangle(2,7,12,index); addTriangle(2,12,4,index); addTriangle(9,8,5,index); break; 
        case 107: addTriangle(9,8,5,index); addTriangle(1,12,4,index); addTriangle(1,7,12,index); addTriangle(1,10,7,index); addTriangle(10,6,7,index); break; 
        case 108: addTriangle(1,3,12,index); addTriangle(1,12,8,index); addTriangle(1,8,5,index); addTriangle(11,6,7,index); break; 
        case 109: addTriangle(2,3,8,index); addTriangle(8,5,2,index); addTriangle(3,12,8,index); addTriangle(5,10,2,index); addTriangle(6,7,11,index); break; 
        case 110: addTriangle(2,6,7,index); addTriangle(2,7,12,index); addTriangle(1,12,8,index); addTriangle(1,8,5,index); addTriangle(1,2,12,index); break; 
        case 111: addTriangle(10,6,7,index); addTriangle(10,7,12,index); addTriangle(12,8,5,index); addTriangle(12,5,10,index); break; 
        case 112: addTriangle(10,9,8,index); addTriangle(10,8,7,index); addTriangle(10,7,11,index); break; 
        case 113: addTriangle(8,7,11,index); addTriangle(9,8,2,index); addTriangle(8,11,2,index); addTriangle(1,9,2,index); break; 
        case 114: addTriangle(3,9,8,index); addTriangle(3,8,7,index); addTriangle(2,9,3,index); addTriangle(2,10,9,index); break; 
        case 115: addTriangle(3,1,9,index); addTriangle(3,9,8,index); addTriangle(3,8,7,index); break; 
        case 116: addTriangle(8,7,11,index); addTriangle(8,11,4,index); addTriangle(10,4,11,index); addTriangle(4,10,1,index); break; 
        case 117: addTriangle(4,8,7,index); addTriangle(4,7,11,index); addTriangle(4,11,2,index); break; 
        case 118: addTriangle(3,2,7,index); addTriangle(2,10,7,index); addTriangle(1,4,8,index); addTriangle(10,8,7,index); addTriangle(10,1,8,index); break; 
        case 119: addTriangle(8,7,3,index); addTriangle(3,4,8,index); break; case 120: addTriangle(10,9,8,index); addTriangle(10,8,7,index); addTriangle(10,7,11,index); addTriangle(3,12,4,index); break;
        case 121: addTriangle(8,7,11,index); addTriangle(9,8,2,index); addTriangle(8,11,2,index); addTriangle(1,9,2,index); addTriangle(3,12,4,index); break; 
        case 122: addTriangle(2,12,4,index); addTriangle(2,7,12,index); addTriangle(10,9,8,index); addTriangle(10,8,7,index); addTriangle(2,10,7,index); break; 
        case 123: addTriangle(1,7,12,index); addTriangle(1,12,4,index); addTriangle(7,1,9,index); addTriangle(7,9,8,index); break; 
        case 124: addTriangle(10,1,8,index); addTriangle(10,8,7,index); addTriangle(7,11,10,index); addTriangle(1,3,12,index); addTriangle(12,8,1,index); break; 
        case 125: addTriangle(2,8,7,index); addTriangle(7,11,2,index); addTriangle(3,12,8,index); addTriangle(8,2,3,index); break; 
        case 126: addTriangle(12,8,7,index); addTriangle(10,1,2,index); break; 
        case 127: addTriangle(12,8,7,index); break; 
        case 128: addTriangle(7,8,12,index); break; 
        case 129: addTriangle(7,8,12,index); addTriangle(2,1,10,index); break; 
        case 130: addTriangle(7,8,12,index); addTriangle(2,11,3,index); break; 
        case 131: addTriangle(10,11,1,index); addTriangle(11,3,1,index); addTriangle(7,8,12,index); break; 
        case 132: addTriangle(1,4,9,index); addTriangle(7,8,12,index); break; 
        case 133: addTriangle(10,2,4,index); addTriangle(4,9,10,index); addTriangle(7,8,12,index); break; 
        case 134: addTriangle(2,11,3,index); addTriangle(1,4,9,index); addTriangle(7,8,12,index); break; 
        case 135: addTriangle(10,4,9,index); addTriangle(10,3,4,index); addTriangle(10,11,3,index); addTriangle(7,8,12,index); break; 
        case 136: addTriangle(3,7,8,index); addTriangle(8,4,3,index); break; 
        case 137: addTriangle(3,7,8,index); addTriangle(8,4,3,index); addTriangle(2,1,10,index); break; 
        case 138: addTriangle(4,2,11,index); addTriangle(4,11,7,index); addTriangle(4,7,8,index); break; 
        case 139: addTriangle(11,7,8,index); addTriangle(8,4,10,index); addTriangle(10,11,8,index); addTriangle(4,1,10,index); break; 
        case 140: addTriangle(3,7,8,index); addTriangle(3,8,9,index); addTriangle(3,9,1,index); break; 
        case 141: addTriangle(3,7,8,index); addTriangle(8,9,3,index); addTriangle(2,3,9,index); addTriangle(9,10,2,index); break; 
        case 142: addTriangle(11,7,8,index); addTriangle(2,11,8,index); addTriangle(2,8,9,index); addTriangle(2,9,1,index); break; 
        case 143: addTriangle(10,11,7,index); addTriangle(10,7,8,index); addTriangle(10,8,9,index); break; 
        case 144: addTriangle(7,8,12,index); addTriangle(10,5,6,index); break; 
        case 145: addTriangle(5,6,2,index); addTriangle(2,1,5,index); addTriangle(7,8,12,index); break; 
        case 146: addTriangle(2,11,3,index); addTriangle(10,5,6,index); addTriangle(7,8,12,index); break; 
        case 147: addTriangle(1,5,6,index); addTriangle(1,6,11,index); addTriangle(1,11,3,index); addTriangle(7,8,12,index); break; 
        case 148: addTriangle(1,4,9,index); addTriangle(10,5,6,index); addTriangle(7,8,12,index); break; 
        case 149: addTriangle(2,4,9,index); addTriangle(2,9,5,index); addTriangle(2,5,6,index); addTriangle(7,8,12,index); break; 
        case 150: addTriangle(7,8,12,index); addTriangle(2,11,3,index); addTriangle(1,4,9,index); addTriangle(10,5,6,index); break; 
        case 151: addTriangle(5,6,11,index); addTriangle(9,5,11,index); addTriangle(11,3,9,index); addTriangle(3,4,9,index); addTriangle(7,8,12,index); break; 
        case 152: addTriangle(3,7,8,index); addTriangle(8,4,3,index); addTriangle(10,5,6,index); break; 
        case 153: addTriangle(5,6,2,index); addTriangle(2,1,5,index); addTriangle(3,7,8,index); addTriangle(8,4,3,index); break;
        case 154: addTriangle(4,2,11,index); addTriangle(4,11,7,index); addTriangle(4,7,8,index); addTriangle(10,5,6,index); break; 
        case 155: addTriangle(11,4,1,index); addTriangle(7,8,4,index); addTriangle(4,11,7,index); addTriangle(1,5,6,index); addTriangle(6,11,1,index); break; 
        case 156: addTriangle(3,7,8,index); addTriangle(3,8,9,index); addTriangle(3,9,1,index); addTriangle(10,5,6,index); break; 
        case 157: addTriangle(9,2,3,index); addTriangle(3,7,8,index); addTriangle(8,9,3,index); addTriangle(5,6,2,index); addTriangle(2,9,5,index); break; 
        case 158: addTriangle(11,7,8,index); addTriangle(2,11,8,index); addTriangle(2,8,9,index); addTriangle(2,9,1,index); break; 
        case 159: addTriangle(8,9,11,index); addTriangle(11,7,8,index); addTriangle(9,6,11,index); addTriangle(9,5,6,index); break; 
        case 160: addTriangle(11,6,8,index); addTriangle(8,12,11,index); break; 
        case 161: addTriangle(11,6,8,index); addTriangle(8,12,11,index); addTriangle(2,1,10,index); break; 
        case 162: addTriangle(6,8,12,index); addTriangle(6,12,3,index); addTriangle(6,3,2,index); break; 
        case 163: addTriangle(8,12,3,index); addTriangle(3,1,6,index); addTriangle(6,8,3,index); addTriangle(1,10,6,index); break; 
        case 164: addTriangle(11,6,8,index); addTriangle(8,12,11,index); addTriangle(1,4,9,index); break; 
        case 165: addTriangle(10,2,4,index); addTriangle(4,9,10,index); addTriangle(11,6,8,index); addTriangle(8,12,11,index); break; 
        case 166: addTriangle(6,8,12,index); addTriangle(6,12,3,index); addTriangle(6,3,2,index); addTriangle(1,4,9,index); break; 
        case 167: addTriangle(6,3,10,index); addTriangle(6,8,12,index); addTriangle(12,3,6,index); addTriangle(4,9,10,index); addTriangle(10,3,4,index); break; 
        case 168: addTriangle(8,4,3,index); addTriangle(8,3,11,index); addTriangle(8,11,6,index); break; 
        case 169: addTriangle(8,4,3,index); addTriangle(8,3,11,index); addTriangle(8,11,6,index); addTriangle(2,1,10,index); break; 
        case 170: addTriangle(2,6,8,index); addTriangle(8,4,2,index); break; 
        case 171: addTriangle(8,4,1,index); addTriangle(8,1,10,index); addTriangle(8,10,6,index); break; 
        case 172: addTriangle(3,11,6,index); addTriangle(6,8,3,index); addTriangle(1,3,8,index); addTriangle(8,9,1,index); break; 
        case 173: addTriangle(9,3,8,index); addTriangle(8,3,11,index); addTriangle(11,6,8,index); addTriangle(2,3,9,index); addTriangle(9,10,2,index); break; 
        case 174: addTriangle(6,8,9,index); addTriangle(6,9,1,index); addTriangle(6,1,2,index); break; 
        case 175: addTriangle(9,10,6,index); addTriangle(6,8,9,index); break; 
        case 176: addTriangle(11,10,5,index); addTriangle(11,5,8,index); addTriangle(11,8,12,index); break; 
        case 177: addTriangle(2,1,5,index); addTriangle(5,11,2,index); addTriangle(11,5,8,index); addTriangle(8,12,11,index); break; 
        case 178: addTriangle(2,10,5,index); addTriangle(2,5,8,index); addTriangle(8,3,2,index); addTriangle(8,12,3,index); break; 
        case 179: addTriangle(1,5,8,index); addTriangle(1,8,12,index); addTriangle(1,12,3,index); break; 
        case 180: addTriangle(11,10,5,index); addTriangle(11,5,8,index); addTriangle(11,8,12,index); addTriangle(1,4,9,index); break; 
        case 181: addTriangle(5,11,2,index); addTriangle(11,5,8,index); addTriangle(8,12,11,index); addTriangle(2,4,9,index); addTriangle(9,5,2,index); break; 
        case 182: addTriangle(2,10,5,index); addTriangle(2,5,8,index); addTriangle(8,3,2,index); addTriangle(8,12,3,index); addTriangle(1,4,9,index); break; 
        case 183: addTriangle(3,4,9,index); addTriangle(9,5,3,index); addTriangle(5,8,12,index); addTriangle(12,3,5,index); break; 
        case 184: addTriangle(11,10,5,index); addTriangle(5,8,11,index); addTriangle(8,4,3,index); addTriangle(3,11,8,index); break; 
        case 185: addTriangle(8,11,5,index); addTriangle(5,11,2,index); addTriangle(2,1,5,index); addTriangle(8,4,3,index); addTriangle(3,11,8,index); break; 
        case 186: addTriangle(4,2,10,index); addTriangle(4,10,5,index); addTriangle(4,5,8,index); break; 
        case 187: addTriangle(8,4,1,index); addTriangle(1,5,8,index); break; 
        case 188: addTriangle(8,3,11,index); addTriangle(8,9,3,index); addTriangle(9,1,3,index); addTriangle(5,8,11,index); addTriangle(11,10,5,index); break; 
        case 189: addTriangle(5,8,9,index); addTriangle(3,11,2,index); break; 
        case 190: addTriangle(2,10,5,index); addTriangle(5,8,2,index); addTriangle(8,9,1,index); addTriangle(1,2,8,index); break; 
        case 191: addTriangle(5,8,9,index); break; 
        case 192: addTriangle(5,9,12,index); addTriangle(12,7,5,index); break; 
        case 193: addTriangle(5,9,12,index); addTriangle(12,7,5,index); addTriangle(2,1,10,index); break; 
        case 194: addTriangle(5,9,12,index); addTriangle(12,7,5,index); addTriangle(2,11,3,index); break; 
        case 195: addTriangle(10,11,1,index); addTriangle(11,3,1,index); addTriangle(5,9,12,index); addTriangle(12,7,5,index); break; 
        case 196: addTriangle(5,1,4,index); addTriangle(5,4,12,index); addTriangle(5,12,7,index); break; 
        case 197: addTriangle(4,12,7,index); addTriangle(7,5,2,index); addTriangle(5,10,2,index); addTriangle(2,4,7,index); break; 
        case 198: addTriangle(5,1,4,index); addTriangle(5,4,12,index); addTriangle(5,12,7,index); addTriangle(2,11,3,index); break; 
        case 199: addTriangle(10,4,5,index); addTriangle(10,11,3,index); addTriangle(3,4,10,index); addTriangle(5,4,12,index); addTriangle(12,7,5,index); break; 
        case 200: addTriangle(7,5,9,index); addTriangle(7,9,4,index); addTriangle(7,4,3,index); break; 
        case 201: addTriangle(7,5,9,index); addTriangle(7,9,4,index); addTriangle(7,4,3,index); addTriangle(2,1,10,index); break; 
        case 202: addTriangle(5,9,4,index); addTriangle(4,7,5,index); addTriangle(4,2,7,index); addTriangle(7,2,11,index); break; 
        case 203: addTriangle(4,11,7,index); addTriangle(11,4,1,index); addTriangle(1,10,11,index); addTriangle(7,5,9,index); addTriangle(9,4,7,index); break; 
        case 204: addTriangle(1,3,7,index); addTriangle(7,5,1,index); break; 
        case 205: addTriangle(7,5,10,index); addTriangle(7,10,2,index); addTriangle(7,2,3,index); break; 
        case 206: addTriangle(5,1,2,index); addTriangle(5,2,11,index); addTriangle(5,11,7,index); break; 
        case 207: addTriangle(11,7,5,index); addTriangle(5,10,11,index); break; 
        case 208: addTriangle(9,12,7,index); addTriangle(9,7,6,index); addTriangle(9,6,10,index); break; 
        case 209: addTriangle(2,1,9,index); addTriangle(9,12,2,index); addTriangle(2,12,6,index); addTriangle(12,7,6,index); break; 
        case 210: addTriangle(9,12,7,index); addTriangle(9,7,6,index); addTriangle(9,6,10,index); addTriangle(2,11,3,index); break; 
        case 211: addTriangle(1,9,6,index); addTriangle(9,12,7,index); addTriangle(7,6,9,index); addTriangle(1,6,11,index); addTriangle(11,3,1,index); break; 
        case 212: addTriangle(4,12,7,index); addTriangle(1,4,7,index); addTriangle(7,6,1,index); addTriangle(6,10,1,index); break; 
        case 213: addTriangle(2,4,12,index); addTriangle(2,12,7,index); addTriangle(2,7,6,index); break; 
        case 214: addTriangle(4,12,7,index); addTriangle(1,4,7,index); addTriangle(7,6,1,index); addTriangle(6,10,1,index); addTriangle(2,11,3,index); break; 
        case 215: addTriangle(3,4,6,index); addTriangle(6,11,3,index); addTriangle(4,12,7,index); addTriangle(7,6,4,index); break; 
        case 216: addTriangle(7,6,10,index); addTriangle(3,7,10,index); addTriangle(10,9,3,index); addTriangle(9,4,3,index); break;
        case 217: addTriangle(9,7,6,index); addTriangle(1,9,6,index); addTriangle(6,2,1,index); addTriangle(4,3,7,index); addTriangle(7,9,4,index); break; 
        case 218: addTriangle(9,4,7,index); addTriangle(6,10,9,index); addTriangle(9,7,6,index); addTriangle(4,2,11,index); addTriangle(11,7,4,index); break; 
        case 219: addTriangle(11,7,6,index); addTriangle(9,4,1,index); break; 
        case 220: addTriangle(3,7,6,index); addTriangle(3,6,10,index); addTriangle(3,10,1,index); break; 
        case 221: addTriangle(6,2,3,index); addTriangle(3,7,6,index); break; 
        case 222: addTriangle(11,7,1,index); addTriangle(1,2,11,index); addTriangle(7,6,10,index); addTriangle(10,1,7,index); break; 
        case 223: addTriangle(11,7,6,index); break;
        case 224: addTriangle(12,11,6,index); addTriangle(12,6,5,index); addTriangle(12,5,9,index); break; 
        case 225: addTriangle(12,11,6,index); addTriangle(12,6,5,index); addTriangle(12,5,9,index); addTriangle(2,1,10,index); break; 
        case 226: addTriangle(6,12,3,index); addTriangle(3,2,6,index); addTriangle(12,6,5,index); addTriangle(5,9,12,index); break; 
        case 227: addTriangle(6,12,3,index); addTriangle(5,9,12,index); addTriangle(12,6,5,index); addTriangle(10,6,3,index); addTriangle(3,1,10,index); break; 
        case 228: addTriangle(1,4,12,index); addTriangle(1,12,11,index); addTriangle(1,11,5,index); addTriangle(11,6,5,index); break; 
        case 229: addTriangle(4,12,5,index); addTriangle(12,11,6,index); addTriangle(6,5,12,index); addTriangle(4,5,10,index); addTriangle(10,2,4,index); break; 
        case 230: addTriangle(12,6,5,index); addTriangle(6,12,3,index); addTriangle(3,2,6,index); addTriangle(5,1,4,index); addTriangle(4,12,5,index); break; 
        case 231: addTriangle(6,5,10,index); addTriangle(4,12,3,index); break; 
        case 232: addTriangle(11,6,5,index); addTriangle(11,5,9,index); addTriangle(9,3,11,index); addTriangle(9,4,3,index); break; 
        case 233: addTriangle(11,6,5,index); addTriangle(11,5,9,index); addTriangle(9,3,11,index); addTriangle(9,4,3,index); addTriangle(2,1,10,index); break; 
        case 234: addTriangle(2,6,5,index); addTriangle(2,5,9,index); addTriangle(2,9,4,index); break; 
        case 235: addTriangle(10,6,4,index); addTriangle(4,1,10,index); addTriangle(9,4,6,index); addTriangle(6,5,9,index); break; 
        case 236: addTriangle(1,3,11,index); addTriangle(1,11,6,index); addTriangle(1,6,5,index); break; 
        case 237: addTriangle(6,5,3,index); addTriangle(3,11,6,index); addTriangle(2,3,5,index); addTriangle(5,10,2,index); break; 
        case 238: addTriangle(2,6,5,index); addTriangle(5,1,2,index); break; 
        case 239: addTriangle(6,5,10,index); break; 
        case 240: addTriangle(11,10,9,index); addTriangle(9,12,11,index); break; 
        case 241: addTriangle(12,11,2,index); addTriangle(12,2,1,index); addTriangle(12,1,9,index); break; 
        case 242: addTriangle(9,12,3,index); addTriangle(9,3,2,index); addTriangle(9,2,10,index); break; 
        case 243: addTriangle(3,1,9,index); addTriangle(9,12,3,index); break; 
        case 244: addTriangle(11,10,1,index); addTriangle(11,1,4,index); addTriangle(11,4,12,index); break; 
        case 245: addTriangle(4,12,11,index); addTriangle(11,2,4,index); break; 
        case 246: addTriangle(2,10,12,index); addTriangle(12,3,2,index); addTriangle(4,12,10,index); addTriangle(10,1,4,index); break; 
        case 247: addTriangle(4,12,3,index); break; 
        case 248: addTriangle(10,9,4,index); addTriangle(10,4,3,index); addTriangle(10,3,11,index); break; 
        case 249: addTriangle(11,9,4,index); addTriangle(4,3,11,index); addTriangle(2,1,9,index); addTriangle(9,11,2,index); break; 
        case 250: addTriangle(4,2,10,index); addTriangle(10,9,4,index); break; 
        case 251: addTriangle(9,4,1,index); break; 
        case 252: addTriangle(1,11,10,index); addTriangle(1,3,11,index); break; 
        case 253: addTriangle(3,11,2,index); break; 
        case 254: addTriangle(10,1,2,index); break;
        }
      }
    //if (verbose)
    //  cout << z << "/" << field.zRes() << " "; flush(cout);
  }

  // create the final triangles based on the vertex indices -- this
  // couldn't be done in the inner loop because the vector keeps
  // resizing and changing the vertex addresses
  assert(_triangleVertices.size() % 3 == 0);
	
  if (verbose) cout << "computed triangles: " << _triangleVertices.size() << endl;
	
  for (unsigned int x = 0; x < _triangleVertices.size() / 3; x++)
  {
    VEC3F* v0 = &_vertices[_triangleVertices[3 * x]];
    VEC3F* v1 = &_vertices[_triangleVertices[3 * x + 1]];
    VEC3F* v2 = &_vertices[_triangleVertices[3 * x + 2]];

    // ignore any degenerate triangles
    Real dist0 = norm((*v0) - (*v1));
    Real dist1 = norm((*v0) - (*v2));
    Real dist2 = norm((*v1) - (*v2));
    Real eps = 1e-7;
    if (dist0 < eps || dist1 < eps || dist2 < eps)
      continue;
    _triangles.push_back(TRIANGLE(v0, v1, v2));
  }

  // all done -- throw away the indices
  _triangleVertices.clear();

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  if (verbose)
    cout << "done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// add vertex triplets to be interpolated
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexTriplets(int i, int j, int k, VEC3I index)
{
  pair<VEC3I, VEC3I> toAdd;
  toAdd = getVertexTriplets(i, index);
  _vertexTriplets[toAdd] = true;

  toAdd = getVertexTriplets(j, index);
  _vertexTriplets[toAdd] = true;

  toAdd = getVertexTriplets(k, index);
  _vertexTriplets[toAdd] = true;
}

//////////////////////////////////////////////////////////////////////
// add vertex pairs to be interpolated
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexPairs(int i, int j, int k, int index)
{
  pair<int, int> toAdd;
  toAdd = getVertexPair(i, index);
  _vertexPairs[toAdd] = true;

  toAdd = getVertexPair(j, index);
  _vertexPairs[toAdd] = true;

  toAdd = getVertexPair(k, index);
  _vertexPairs[toAdd] = true;
}

//////////////////////////////////////////////////////////////////////
// get the indices for the first and second vertices in the
// interpolation
//////////////////////////////////////////////////////////////////////
pair<int, int> TRIANGLE_MESH::getVertexPair(int i, int index)
{
  pair<int, int> toAdd;
  toAdd.first = index;
  toAdd.second = index;

  // in general: P = (coord+1), N = (coord)
	switch (i) {
	case 1: // (x,y,z) to (x,y+1,z)
    toAdd.first += 0;
    toAdd.second += _xRes;
		break;
	case 2: // (x,y,z) to (x,y,z+1)
    toAdd.first += 0;
    toAdd.second += _slabSize;
		break;
	case 3: // (x,y,z+1) to (x,y+1,z+1)
    toAdd.first += _slabSize;
    toAdd.second += _xRes + _slabSize;
		break;
	case 4: // (x,y+1,z) to (x,y+1,z+1)
    toAdd.first += _xRes;
    toAdd.second += _xRes + _slabSize;
		break;
	case 5: // (x+1,y,z) to (x+1,y+1,z)
    toAdd.first += 1;
    toAdd.second += 1 + _xRes;
		break;
	case 6: // (x+1,y,z) to (x+1,y,z+1)
    toAdd.first += 1;
    toAdd.second += 1 + _slabSize;
		break;
	case 7: // (x+1,y,z+1) to (x+1,y+1,z+1)
    toAdd.first += 1 + _slabSize;
    toAdd.second += 1 + _xRes + _slabSize;
		break;
	case 8: // (x+1,y+1,z) to (x+1,y+1,z+1)
    toAdd.first += 1 + _xRes;
    toAdd.second += 1 + _xRes + _slabSize;
		break;
	case 9: // (x,y+1,z) to (x+1,y+1,z)
    toAdd.first += _xRes;
    toAdd.second += 1 + _xRes;
		break;
	case 10: // (x,y,z) to (x+1,y,z)
    toAdd.first += 0;
    toAdd.second += 1;
		break;
	case 11: // (x,y,z+1) to (x+1,y,z+1)
    toAdd.first += _slabSize;
    toAdd.second += 1 + _slabSize;
		break;
	case 12: // (x,y+1,z+1) to (x+1,y+1,z+1)
    toAdd.first += _xRes + _slabSize;
    toAdd.second += 1 + _xRes + _slabSize;
		break;
	}

  return toAdd;  
}

//////////////////////////////////////////////////////////////////////
// get the indices for the first and second vertices in the
// interpolation
//////////////////////////////////////////////////////////////////////
pair<VEC3I, VEC3I> TRIANGLE_MESH::getVertexTriplets(int i, VEC3I v)
{
  //pair<vector<int>, vector<int> > toAdd;
  pair<VEC3I, VEC3I> toAdd;
  toAdd.first = v;
  toAdd.second = v;

  // in general: P = (coord+1), N = (coord)
	switch (i) {
	case 1: // (x,y,z) to (x,y+1,z)
    //toAdd.first += 0;
    //toAdd.second += _xRes;
    toAdd.second[1] += 1;
		break;
	case 2: // (x,y,z) to (x,y,z+1)
    //toAdd.first += 0;
    //toAdd.second += _slabSize;
    toAdd.second[2] += 1;
		break;
	case 3: // (x,y,z+1) to (x,y+1,z+1)
    //toAdd.first += _slabSize;
    //toAdd.second += _xRes + _slabSize;
    toAdd.first[2] += 1;
    toAdd.second[1] += 1;
    toAdd.second[2] += 1;
		break;
	case 4: // (x,y+1,z) to (x,y+1,z+1)
    //toAdd.first += _xRes;
    //toAdd.second += _xRes + _slabSize;
    toAdd.first[1] += 1;
    toAdd.second[1] += 1;
    toAdd.second[2] += 1;
		break;
	case 5: // (x+1,y,z) to (x+1,y+1,z)
    //toAdd.first += 1;
    //toAdd.second += 1 + _xRes;
    toAdd.first[0] += 1;
    toAdd.second[0] += 1;
    toAdd.second[1] += 1;
		break;
	case 6: // (x+1,y,z) to (x+1,y,z+1)
    //toAdd.first += 1;
    //toAdd.second += 1 + _slabSize;
    toAdd.first[0] += 1;
    toAdd.second[0] += 1;
    toAdd.second[2] += 1;
		break;
	case 7: // (x+1,y,z+1) to (x+1,y+1,z+1)
    //toAdd.first += 1 + _slabSize;
    //toAdd.second += 1 + _xRes + _slabSize;
    toAdd.first[0] += 1;
    toAdd.first[2] += 1;
    toAdd.second[0] += 1;
    toAdd.second[1] += 1;
    toAdd.second[2] += 1;
		break;
	case 8: // (x+1,y+1,z) to (x+1,y+1,z+1)
    //toAdd.first += 1 + _xRes;
    //toAdd.second += 1 + _xRes + _slabSize;
    toAdd.first[0] += 1;
    toAdd.first[1] += 1;
    toAdd.second[0] += 1;
    toAdd.second[1] += 1;
    toAdd.second[2] += 1;
		break;
	case 9: // (x,y+1,z) to (x+1,y+1,z)
    //toAdd.first += _xRes;
    //toAdd.second += 1 + _xRes;
    toAdd.first[1] += 1;
    toAdd.second[0] += 1;
    toAdd.second[1] += 1;
		break;
	case 10: // (x,y,z) to (x+1,y,z)
    //toAdd.first += 0;
    //toAdd.second += 1;
    toAdd.second[0] += 1;
		break;
	case 11: // (x,y,z+1) to (x+1,y,z+1)
    //toAdd.first += _slabSize;
    //toAdd.second += 1 + _slabSize;
    toAdd.first[2] += 1;
    toAdd.second[0] += 1;
    toAdd.second[2] += 1;
		break;
	case 12: // (x,y+1,z+1) to (x+1,y+1,z+1)
    //toAdd.first += _xRes + _slabSize;
    //toAdd.second += 1 + _xRes + _slabSize;
    toAdd.first[1] += 1;
    toAdd.first[2] += 1;
    toAdd.second[0] += 1;
    toAdd.second[1] += 1;
    toAdd.second[2] += 1;
		break;
	}

  return toAdd;  
}

//////////////////////////////////////////////////////////////////////
// get the (x,y,z) of an index
//////////////////////////////////////////////////////////////////////
VEC3I TRIANGLE_MESH::getXYZ(const int index) const
{
  VEC3I final;
  final[2] = index / _slabSize;
  final[1] = (index % _slabSize) / _xRes;
  final[0] = (index % _slabSize) % _xRes;

  assert(final[0] + final[1] * _xRes + final[2] * _slabSize == index);
  return final;
}

//////////////////////////////////////////////////////////////////////
// get the nonlinear function value here
//////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::nonlinearValue(const VEC3F& position, const bool debug)
{
  QUATERNION iterate(position[0], position[1], position[2], _quaternionSlice);
 
  // TODO: fractal should really be passing this in ...
  const int maxIterations = _maxIterations;
  //Real escape = 20.0;
  const Real escape = _escapeRadius;

  Real magnitude = iterate.magnitude();
  int totalIterations = 0;
  while (magnitude < escape && totalIterations < maxIterations)
  {
    if (debug)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " iteration " << totalIterations << ": " << iterate << endl;
    }
    QUATERNION topEval = _top.evaluateScaledPowerFactored(iterate);
    QUATERNION bottomEval;
    
    if (_bottom.totalRoots() > 0)
    {
      bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
      iterate = (topEval / bottomEval);
    }
    else
      iterate = topEval;

    iterate *= _expScaling;
    magnitude = iterate.magnitude();
    totalIterations++;

    // see if it fell into the black hole at the origin
    if (magnitude < 10.0 * REAL_MIN)
      totalIterations = maxIterations;
  }
  return log(magnitude) - _isosurface;
}

//////////////////////////////////////////////////////////////////////
// get the quadratic function value here
//////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::quadraticValue(const VEC3F& position, const bool debug)
{
  QUATERNION iterate(position[0], position[1], position[2], 0);
 
  // TODO: fractal should really be passing this in ...
  int maxIterations = _maxIterations;
  //Real escape = 2.0;
  Real escape = _escapeRadius;

  Real magnitude = iterate.magnitude();
  int totalIterations = 0;
  while (magnitude < escape && totalIterations < maxIterations)
  {
    iterate = iterate * iterate + _quadraticConst;
    magnitude = iterate.magnitude();

    totalIterations++;
  }
  return magnitude < escape ? -1 : 1;
}

//////////////////////////////////////////////////////////////////////
// compute the linear interpolations for matching cubes edges
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearEdgeInterpolations()
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing non-linear edge interpolations ... " << flush;
  map<pair<int, int>, bool>::iterator iter;
 
  // flatten the map out to an array now that collision are resolved
  vector<pair<int, int> > pairs;
  for (iter = _vertexPairs.begin(); iter != _vertexPairs.end(); iter++)
  {
    pair<int,int> vertexPair = iter->first;
    pairs.push_back(vertexPair);
  }

  // compute the actual vertices
  _vertices.clear();
  _vertices.resize(pairs.size());
  const int size = pairs.size();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int x = 0; x < size; x++)
  {
    pair<int,int> vertexPair = pairs[x];

    int firstIndex = vertexPair.first;
    VEC3I firstXYZ = getXYZ(firstIndex);
    VEC3F firstVertex = cellCenter(firstXYZ[0], firstXYZ[1], firstXYZ[2]);
    Real firstValue = nonlinearValue(firstVertex, false);
    
    int secondIndex = vertexPair.second;
    VEC3I secondXYZ = getXYZ(secondIndex);
    VEC3F secondVertex = cellCenter(secondXYZ[0], secondXYZ[1], secondXYZ[2]);
    Real secondValue = nonlinearValue(secondVertex, false);

    VEC3F positiveVertex = firstVertex;
    Real positiveValue = firstValue;
    VEC3F negativeVertex = secondVertex;
    Real negativeValue = secondValue;

    if (positiveValue * negativeValue >= 0.0)
    {
    #pragma omp critical
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " field dims: " << endl;
        cout << " res: " << _res << endl;
        cout << " lengths: " << _lengths << endl;
        cout << " center:  " << _center << endl;
        cout << " dxs:     " << _dxs << endl;
        cout << " p vertex: " << positiveVertex << " n vertex: " << negativeVertex << endl;
        cout << " first XYZ: " << firstXYZ << " second XYZ:" << secondXYZ << endl;
        cout << " firstIndex: " << firstIndex << endl;
        cout << " secondIndex: " << secondIndex << endl;

        cout << " positive:           " << positiveValue                                << " negative:           " << negativeValue << endl;

        exit(0);
      }
    }

    //assert(positiveValue * negativeValue < 0.0);

    if (firstValue < 0)
    {
      positiveVertex = secondVertex;
      positiveValue = secondValue;
      negativeVertex = firstVertex;
      negativeValue = firstValue;
    }

    // this turns the midpoint search on and off. If you want to compare to just traditional
    // marching cubes with linear interpolation, set this to 0.
#if 1
    VEC3F finalVertex = midpointSearch(positiveVertex, positiveValue, negativeVertex, negativeValue);
#else
    VEC3F offset = secondVertex - firstVertex;
   
    // get the distances at the two cell centers 
    //Real dist[] = {field(firstXYZ[0], firstXYZ[1], firstXYZ[2]),
    //               field(secondXYZ[0], secondXYZ[1], secondXYZ[2])};
    Real dist[] = {nonlinearValue(firstVertex), nonlinearValue(secondVertex)};
    Real interpolation = dist[0] / (dist[0] - dist[1]);

    // store the final vertex for this edge
    VEC3F finalVertex = firstVertex + interpolation * offset;
#endif 
    _vertices[x] = finalVertex;
    if (x % (int)(size / 10) == 0)
      cout << 100 * ((Real)x / size) << "% " << flush;
  }

  // hash where each vertex is so the triangle construction looking for it later can find it
  _vertexPairHash.clear();
  for (unsigned int x = 0; x < pairs.size(); x++)
  {
    pair<int,int> vertexPair = pairs[x];
    int firstIndex = vertexPair.first;
    int secondIndex = vertexPair.second;
    _vertexPairHash[pair<int,int>(firstIndex, secondIndex)] = x;
    _vertexPairHash[pair<int,int>(secondIndex, firstIndex)] = x;
  }
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// compute the linear interpolations for matching cubes edges
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearEdgeInterpolationsHuge()
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing huge non-linear edge interpolations ... " << flush;
  map<pair<VEC3I, VEC3I>, bool>::iterator iter;
 
  // flatten the map out to an array now that collision are resolved
  vector<pair<VEC3I, VEC3I> > pairs;
  for (iter = _vertexTriplets.begin(); iter != _vertexTriplets.end(); iter++)
  {
    pair<VEC3I,VEC3I> vertexPair = iter->first;
    pairs.push_back(vertexPair);
  }

  // compute the actual vertices
  _vertices.clear();
  _vertices.resize(pairs.size());
  const int size = pairs.size();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int x = 0; x < size; x++)
  {
    pair<VEC3I,VEC3I>& vertexPair = pairs[x];

    VEC3I& firstXYZ = vertexPair.first;
    VEC3F firstVertex = cellCenter(firstXYZ[0], firstXYZ[1], firstXYZ[2]);
    Real firstValue = nonlinearValue(firstVertex, false);
    
    VEC3I& secondXYZ = vertexPair.second;
    VEC3F secondVertex = cellCenter(secondXYZ[0], secondXYZ[1], secondXYZ[2]);
    Real secondValue = nonlinearValue(secondVertex, false);

    VEC3F positiveVertex = firstVertex;
    Real positiveValue = firstValue;
    VEC3F negativeVertex = secondVertex;
    Real negativeValue = secondValue;

    if (positiveValue * negativeValue >= 0.0)
    {
    #pragma omp critical
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " field dims: " << endl;
        cout << " res: " << _res << endl;
        cout << " lengths: " << _lengths << endl;
        cout << " center:  " << _center << endl;
        cout << " dxs:     " << _dxs << endl;
        cout << " p vertex: " << positiveVertex << " n vertex: " << negativeVertex << endl;
        cout << " first XYZ: " << firstXYZ << " second XYZ:" << secondXYZ << endl;

        cout << " positive:           " << positiveValue                                << " negative:           " << negativeValue << endl;

        exit(0);
      }
    }

    //assert(positiveValue * negativeValue < 0.0);

    if (firstValue < 0)
    {
      positiveVertex = secondVertex;
      positiveValue = secondValue;
      negativeVertex = firstVertex;
      negativeValue = firstValue;
    }

    // this turns the midpoint search on and off. If you want to compare to just traditional
    // marching cubes with linear interpolation, set this to 0.
    VEC3F finalVertex = midpointSearch(positiveVertex, positiveValue, negativeVertex, negativeValue);
    
    _vertices[x] = finalVertex;
    if (x % (int)(size / 10) == 0)
      cout << 100 * ((Real)x / size) << "% " << flush;
  }

  // hash where each vertex is so the triangle construction looking for it later can find it
  _vertexTripletHash.clear();
  for (unsigned int x = 0; x < pairs.size(); x++)
  {
    pair<VEC3I,VEC3I> vertexTriplet = pairs[x];
    VEC3I firstIndex = vertexTriplet.first;
    VEC3I secondIndex = vertexTriplet.second;
    _vertexTripletHash[pair<VEC3I,VEC3I>(firstIndex, secondIndex)] = x;
    _vertexTripletHash[pair<VEC3I,VEC3I>(secondIndex, firstIndex)] = x;
  }
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// compute the linear interpolations for matching cubes edges
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeQuadraticEdgeInterpolations()
{
  cout << " Computing quadratic edge interpolations ... " << flush;
  const FIELD_3D& field = *_toMarch;
  map<pair<int, int>, bool>::iterator iter;
 
  // flatten the map out to an array now that collision are resolved
  vector<pair<int, int> > pairs;
  for (iter = _vertexPairs.begin(); iter != _vertexPairs.end(); iter++)
  {
    pair<int,int> vertexPair = iter->first;
    pairs.push_back(vertexPair);
  }

  // compute the actual vertices
  _vertices.clear();
  _vertices.resize(pairs.size());
  const int size = pairs.size();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int x = 0; x < size; x++)
  {
    pair<int,int> vertexPair = pairs[x];

    int firstIndex = vertexPair.first;
    VEC3I firstXYZ = getXYZ(firstIndex);
    VEC3F firstVertex = field.cellCenter(firstXYZ[0], firstXYZ[1], firstXYZ[2]);
    Real firstValue = quadraticValue(firstVertex, false);
    
    int secondIndex = vertexPair.second;
    VEC3I secondXYZ = getXYZ(secondIndex);
    VEC3F secondVertex = field.cellCenter(secondXYZ[0], secondXYZ[1], secondXYZ[2]);
    Real secondValue = quadraticValue(secondVertex, false);

    VEC3F positiveVertex = firstVertex;
    Real positiveValue = firstValue;
    VEC3F negativeVertex = secondVertex;
    Real negativeValue = secondValue;

    if (positiveValue * negativeValue >= 0.0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " field dims: " << endl;
      cout << " res: " << field.xRes() << " " << field.yRes() << " " << field.zRes() << endl;
      cout << " lengths: " << field.lengths() << endl;
      cout << " center:  " << field.center() << endl;
      cout << " dxs:     " << field.dxs() << endl;
      cout << " p vertex: " << positiveVertex << " n vertex: " << negativeVertex << endl;
      cout << " first XYZ: " << firstXYZ << " second XYZ:" << secondXYZ << endl;

      cout << " positive:           " << positiveValue                                << " negative:           " << negativeValue << endl;
      cout << " positive should be: " << field(firstXYZ[0], firstXYZ[1], firstXYZ[2]) << " negative should be: " << field(secondXYZ[0], secondXYZ[1], secondXYZ[2]) << endl;

      exit(0);
    }

    //assert(positiveValue * negativeValue < 0.0);

    if (firstValue < 0)
    {
      positiveVertex = secondVertex;
      positiveValue = secondValue;
      negativeVertex = firstVertex;
      negativeValue = firstValue;
    }

    // this turns the midpoint search on and off. If you want to compare to just traditional
    // marching cubes with linear interpolation, set this to 0.
#if 1
    VEC3F finalVertex = quadraticMidpointSearch(positiveVertex, positiveValue, negativeVertex, negativeValue);
#else
    VEC3F offset = secondVertex - firstVertex;
   
    // get the distances at the two cell centers 
    //Real dist[] = {field(firstXYZ[0], firstXYZ[1], firstXYZ[2]),
    //               field(secondXYZ[0], secondXYZ[1], secondXYZ[2])};
    Real dist[] = {quadraticValue(firstVertex), quadraticValue(secondVertex)};
    Real interpolation = dist[0] / (dist[0] - dist[1]);

    // store the final vertex for this edge
    VEC3F finalVertex = firstVertex + interpolation * offset;
#endif 
    _vertices[x] = finalVertex;
    if (x % (int)(size / 10) == 0)
      cout << 100 * ((Real)x / size) << "% " << flush;
  }

  // hash where each vertex is so the triangle construction looking for it later can find it
  _vertexPairHash.clear();
  for (unsigned int x = 0; x < pairs.size(); x++)
  {
    pair<int,int> vertexPair = pairs[x];
    int firstIndex = vertexPair.first;
    int secondIndex = vertexPair.second;
    _vertexPairHash[pair<int,int>(firstIndex, secondIndex)] = x;
    _vertexPairHash[pair<int,int>(secondIndex, firstIndex)] = x;
  }
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// do a midpoint search
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::midpointSearch(const VEC3F& positiveVertex, const Real& positiveValue, 
                                    const VEC3F& negativeVertex, const Real& negativeValue, const int recursion)
{
  VEC3F midpointVertex = positiveVertex + negativeVertex;
  midpointVertex *= 0.5;
  if (positiveValue * negativeValue >= 0.0)
  {
    cout << " positive: " << positiveValue << " negative: " << negativeValue << endl;
    cout << " p vertex: " << positiveVertex << " n vertex: " << negativeVertex << endl;
    cout << " recursion: " << recursion << endl;

    return midpointVertex;
  }

  // this resolves roughly millimeter level details if a grid edge is a meter
  //if (recursion >= 8)
  if (recursion >= 6)
  {
    return midpointVertex;
  }

  Real midpointValue = nonlinearValue(midpointVertex);
  if (fabs(midpointValue) < 1e-8)
  {
    return midpointVertex;
  }

  if (midpointValue < 0)
    return midpointSearch(positiveVertex, positiveValue,
                          midpointVertex, midpointValue, recursion + 1);

  return midpointSearch(midpointVertex, midpointValue,
                        negativeVertex, negativeValue, recursion + 1);
}

//////////////////////////////////////////////////////////////////////
// do a midpoint search
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::quadraticMidpointSearch(const VEC3F& positiveVertex, const Real& positiveValue, 
                                             const VEC3F& negativeVertex, const Real& negativeValue, const int recursion)
{
  if (positiveValue * negativeValue >= 0.0)
  {
    cout << " positive: " << positiveValue << " negative: " << negativeValue << endl;
    cout << " p vertex: " << positiveVertex << " n vertex: " << negativeVertex << endl;
    cout << " recursion: " << recursion << endl;
    exit(0);
  }

  VEC3F midpointVertex = positiveVertex + negativeVertex;
  midpointVertex *= 0.5;

  if (recursion >= 100)
  {
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //exit(0);
    return midpointVertex;
  }

  Real midpointValue = quadraticValue(midpointVertex);
  if (fabs(midpointValue) < 1e-8)
  {
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout << " recursions: " << recursion << endl;
    //exit(0);
    return midpointVertex;
  }

  //cout << endl;
  //cout << " positive: " << positiveVertex << " value: " << positiveValue << endl;
  //cout << " negative: " << negativeVertex << " value: " << negativeValue << endl;
  //cout << " midpoint: " << midpointVertex << " value: " << midpointValue << endl;

  if (midpointValue < 0)
    return quadraticMidpointSearch(positiveVertex, positiveValue,
                          midpointVertex, midpointValue, recursion + 1);

  return quadraticMidpointSearch(midpointVertex, midpointValue,
                        negativeVertex, negativeValue, recursion + 1);
}

//////////////////////////////////////////////////////////////////////
// compute the linear interpolations for matching cubes edges
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdgeInterpolations()
{
  const FIELD_3D& field = *_toMarch;
  map<pair<int, int>, bool>::iterator iter;
 
  // flatten the map out to an array now that collision are resolved
  vector<pair<int, int> > pairs;
  for (iter = _vertexPairs.begin(); iter != _vertexPairs.end(); iter++)
  {
    pair<int,int> vertexPair = iter->first;
    pairs.push_back(vertexPair);
  }

  // compute the actual vertices
  _vertices.clear();
  _vertices.resize(pairs.size());
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int x = 0; x < (int)pairs.size(); x++)
  {
    pair<int,int> vertexPair = pairs[x];

    int firstIndex = vertexPair.first;
    VEC3I firstXYZ = getXYZ(firstIndex);
    VEC3F firstVertex = field.cellCenter(firstXYZ[0], firstXYZ[1], firstXYZ[2]);
    
    int secondIndex = vertexPair.second;
    VEC3I secondXYZ = getXYZ(secondIndex);
    VEC3F secondVertex = field.cellCenter(secondXYZ[0], secondXYZ[1], secondXYZ[2]);

    VEC3F offset = secondVertex - firstVertex;
   
    // get the distances at the two cell centers 
    Real dist[] = {field(firstXYZ[0], firstXYZ[1], firstXYZ[2]),
                   field(secondXYZ[0], secondXYZ[1], secondXYZ[2])};
    Real interpolation = dist[0] / (dist[0] - dist[1]);

    // store the final vertex for this edge
    VEC3F finalVertex = firstVertex + interpolation * offset;
    _vertices[x] = finalVertex;
  }

  // hash where each vertex is so the triangle construction looking for it later can find it
  _vertexPairHash.clear();
  for (unsigned int x = 0; x < pairs.size(); x++)
  {
    pair<int,int> vertexPair = pairs[x];
    int firstIndex = vertexPair.first;
    int secondIndex = vertexPair.second;
    _vertexPairHash[pair<int,int>(firstIndex, secondIndex)] = x;
    _vertexPairHash[pair<int,int>(secondIndex, firstIndex)] = x;
  }
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addStagedTriangle(int i, int j, int k, int index, const CUBE& cube, const VEC3F& center)
{
  VEC3F p0 = computeStagedVertex(i, index, cube, center);
  VEC3F p1 = computeStagedVertex(j, index, cube, center);
  VEC3F p2 = computeStagedVertex(k, index, cube, center);
  if (p0[0] == _outside || p1[0] == _outside || p2[0] == _outside) return;
 
  // if the vertex has been computed before, don't duplicate it 
  int v0 = storeVertex(p0, index);
  int v1 = storeVertex(p1, index);
  int v2 = storeVertex(p2, index);

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addTriangle(int i, int j, int k, int index)
{
  VEC3F p0 = computeVertex(i, index);
  VEC3F p1 = computeVertex(j, index);
  VEC3F p2 = computeVertex(k, index);
  if (p0[0] == _outside || p1[0] == _outside || p2[0] == _outside) return;
 
  // if the vertex has been computed before, don't duplicate it 
  int v0 = storeVertex(p0, index);
  int v1 = storeVertex(p1, index);
  int v2 = storeVertex(p2, index);

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list whose vertices were all precomputed
// by computeEdgeInterpolations()
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexTripletTriangle(int i, int j, int k, VEC3I index)
{
  pair<VEC3I,VEC3I> vp0 = getVertexTriplets(i, index);
  pair<VEC3I,VEC3I> vp1 = getVertexTriplets(j, index);
  pair<VEC3I,VEC3I> vp2 = getVertexTriplets(k, index);

  int v0 = _vertexTripletHash[vp0];
  int v1 = _vertexTripletHash[vp1];
  int v2 = _vertexTripletHash[vp2];

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list whose vertices were all precomputed
// by computeEdgeInterpolations()
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexPairTriangle(int i, int j, int k, int index)
{
  pair<int,int> vp0 = getVertexPair(i, index);
  pair<int,int> vp1 = getVertexPair(j, index);
  pair<int,int> vp2 = getVertexPair(k, index);

  int v0 = _vertexPairHash[vp0];
  int v1 = _vertexPairHash[vp1];
  int v2 = _vertexPairHash[vp2];

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// get the edge point
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::computeStagedVertex(int i, int index, const CUBE& cube, const VEC3F& center) 
{
  VEC3F point(_outside);
  Real dist[2];

  // the base is the current cell center
  const VEC3F& base = center;

  // in general: P = (coord+1), N = (coord)
	switch (i) {
	case 1: // (x,y,z) to (x,y+1,z)
		point[0] = base[0];
    dist[0] = cube.NNN;
    dist[1] = cube.NPN;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2];
		break;
	case 2: // (x,y,z) to (x,y,z+1)
		point[0] = base[0];
		point[1] = base[1];
    dist[0] = cube.NNN;
    dist[1] = cube.NNP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 3: // (x,y,z+1) to (x,y+1,z+1)
		point[0] = base[0];
    dist[0] = cube.NNP;
    dist[1] = cube.NPP;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	case 4: // (x,y+1,z) to (x,y+1,z+1)
		point[0] = base[0];
		point[1] = base[1] + _fieldDeltas[1];
    dist[0] = cube.NPN;
    dist[1] = cube.NPP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 5: // (x+1,y,z) to (x+1,y+1,z)
		point[0] = base[0] + _fieldDeltas[0];
    dist[0] = cube.PNN;
    dist[1] = cube.PPN;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2];
		break;
	case 6: // (x+1,y,z) to (x+1,y,z+1)
		point[0] = base[0] + _fieldDeltas[0];
		point[1] = base[1];
    dist[0] = cube.PNN;
    dist[1] = cube.PNP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 7: // (x+1,y,z+1) to (x+1,y+1,z+1)
		point[0] = base[0] + _fieldDeltas[0];
    dist[0] = cube.PNP;
    dist[1] = cube.PPP;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	case 8: // (x+1,y+1,z) to (x+1,y+1,z+1)
		point[0] = base[0] + _fieldDeltas[0];
		point[1] = base[1] + _fieldDeltas[1];
    dist[0] = cube.PPN;
    dist[1] = cube.PPP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 9: // (x,y+1,z) to (x+1,y+1,z)
    dist[0] = cube.NPN;
    dist[1] = cube.PPN;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1] + _fieldDeltas[1];
		point[2] = base[2];
		break;
	case 10: // (x,y,z) to (x+1,y,z)
    dist[0] = cube.NNN;
    dist[1] = cube.PNN;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1];
		point[2] = base[2];
		break;
	case 11: // (x,y,z+1) to (x+1,y,z+1)
    dist[0] = cube.NNP;
    dist[1] = cube.PNP;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	case 12: // (x,y+1,z+1) to (x+1,y+1,z+1)
    dist[0] = cube.NPP;
    dist[1] = cube.PPP;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1] + _fieldDeltas[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	}
  return point;
}

//////////////////////////////////////////////////////////////////////
// get the edge point
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::computeVertex(int i, int index) 
{
  VEC3F point(_outside);
  Real dist[2];

  // the base is the current cell center
  const VEC3F& base = _cellCenter;

	switch (i) {
	case 1:
		point[0] = base[0];
    dist[0] = _NNN;
    dist[1] = _NPN;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2];
		break;
	case 2:
		point[0] = base[0];
		point[1] = base[1];
    dist[0] = _NNN;
    dist[1] = _NNP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 3:
		point[0] = base[0];
    dist[0] = _NNP;
    dist[1] = _NPP;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	case 4:
		point[0] = base[0];
		point[1] = base[1] + _fieldDeltas[1];
    dist[0] = _NPN;
    dist[1] = _NPP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 5:
		point[0] = base[0] + _fieldDeltas[0];
    dist[0] = _PNN;
    dist[1] = _PPN;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2];
		break;
	case 6:
		point[0] = base[0] + _fieldDeltas[0];
		point[1] = base[1];
    dist[0] = _PNN;
    dist[1] = _PNP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 7:
		point[0] = base[0] + _fieldDeltas[0];
    dist[0] = _PNP;
    dist[1] = _PPP;
		point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	case 8:
		point[0] = base[0] + _fieldDeltas[0];
		point[1] = base[1] + _fieldDeltas[1];
    dist[0] = _PPN;
    dist[1] = _PPP;
		point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
		break;
	case 9:
    dist[0] = _NPN;
    dist[1] = _PPN;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1] + _fieldDeltas[1];
		point[2] = base[2];
		break;
	case 10:
    dist[0] = _NNN;
    dist[1] = _PNN;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1];
		point[2] = base[2];
		break;
	case 11:
    dist[0] = _NNP;
    dist[1] = _PNP;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	case 12:
    dist[0] = _NPP;
    dist[1] = _PPP;
		point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
		point[1] = base[1] + _fieldDeltas[1];
		point[2] = base[2] + _fieldDeltas[2];
		break;
	}
  return point;
}

//////////////////////////////////////////////////////////////////////
// see if the vertex has been computed before, and if not, store it
//////////////////////////////////////////////////////////////////////
int TRIANGLE_MESH::storeVertex(VEC3F& vertex, int index)
{
#if 1
  // get vertices that have been created near this cell before
  vector<int>& nearest = _vertexHash[index];

  // check them all to see if any are near
  for (unsigned int x = 0; x < nearest.size(); x++)
  {
    int nearestIndex = nearest[x];
    VEC3F diff = _vertices[nearestIndex] - vertex;
    Real magnitude = norm(diff);

    // if they're really close to each other, just return that.
    if (magnitude < 1e-6)
      return nearestIndex;
  }

  // go ahead and add it
  _vertices.push_back(vertex);
 
  // add it to the hash table as well 
  int vectorIndex = _vertices.size() - 1;

  int fatten = 1;
  for (int z = -fatten; z <= fatten; z++)
    for (int y = -fatten; y <= fatten; y++)
      for (int x = -fatten; x <= fatten; x++)
        _vertexHash[index + x + _xRes * y + _slabSize * z].push_back(vectorIndex);

  return vectorIndex;
#else
  // this does a slow linear search, but if it becomes a problem,
  // the obvious thing to do is to put a KD-tree search in instead
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    VEC3F diff = _vertices[x] - vertex;
    Real magnitude = norm(diff);

    // if they're really close to each other, just return that.
    if (magnitude < 1e-7)
      return x;
  }

  _vertices.push_back(vertex);

  return _vertices.size() - 1;
#endif
}

//////////////////////////////////////////////////////////////////////
// read an obj file using stdlib
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeOBJ(const string& filename)
{
  FILE* file = NULL;

  file = fopen(filename.c_str(), "w");
  cout << " Writing OBJ file " << filename.c_str() << " ... " << flush;

	for (unsigned int i = 0; i < _vertices.size(); i++)
  {
    //double vertex[] = {_vertices[i][0], _vertices[i][1], _vertices[i][2]};
    double vertex[] = {(double)_vertices[i][0], 
                       (double)_vertices[i][1], 
                       (double)_vertices[i][2]};
    fprintf(file, "v %.16f %.16f %.16f\n", vertex[0], vertex[1], vertex[2]);
		//out << "v " << _vertices[i][0] << " " << _vertices[i][1] << " " << _vertices[i][2] << endl;
  }
	// normals
	for (unsigned int i = 0; i < _normals.size(); i++)
  {
    //double normal[] = {_normals[i][0], _normals[i][1], _normals[i][2]};
    double normal[] = {(double)_normals[i][0], 
                       (double)_normals[i][1], 
                       (double)_normals[i][2]};
    fprintf(file, "vn %.16f %.16f %.16f\n", normal[0], normal[1], normal[2]);
  }
	// faces
	for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    int indices[] = { _vertexIndices[_triangles[i].vertex(0)] + 1,
                  	  _vertexIndices[_triangles[i].vertex(1)] + 1,
                  	  _vertexIndices[_triangles[i].vertex(2)] + 1};

    fprintf(file, "f %i %i %i\n", indices[0], indices[1], indices[2]);
	}

  fclose(file);
  cout << "done. " << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////
// read an obj file using stdlib
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeOBJgz(const string& filename)
{
  string filenameFinal(filename);
  int size = filenameFinal.size();
  if (filenameFinal[size - 1] != 'z' || filenameFinal[size - 2] != 'g')
    filenameFinal = filenameFinal + string(".gz");

  gzFile file ;
  file = gzopen(filenameFinal.c_str(), "w");
  cout << " Writing gzipped OBJ file " << filenameFinal.c_str() << " ... " << flush;

	for (unsigned int i = 0; i < _vertices.size(); i++)
  {
    //double vertex[] = {_vertices[i][0], _vertices[i][1], _vertices[i][2]};
    double vertex[] = {(double)_vertices[i][0], 
                       (double)_vertices[i][1], 
                       (double)_vertices[i][2]};
    gzprintf(file, "v %.16f %.16f %.16f\n", vertex[0], vertex[1], vertex[2]);
		//out << "v " << _vertices[i][0] << " " << _vertices[i][1] << " " << _vertices[i][2] << endl;
  }
	// normals
	for (unsigned int i = 0; i < _normals.size(); i++)
  {
    //double normal[] = {_normals[i][0], _normals[i][1], _normals[i][2]};
    double normal[] = {(double)_normals[i][0], 
                       (double)_normals[i][1], 
                       (double)_normals[i][2]};
    gzprintf(file, "vn %.16f %.16f %.16f\n", normal[0], normal[1], normal[2]);
  }
	// faces
	for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    int indices[] = { _vertexIndices[_triangles[i].vertex(0)] + 1,
                  	  _vertexIndices[_triangles[i].vertex(1)] + 1,
                  	  _vertexIndices[_triangles[i].vertex(2)] + 1};

    gzprintf(file, "f %i %i %i\n", indices[0], indices[1], indices[2]);
	}

  gzclose(file);
  cout << "done. " << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////
// read an obj file using stdlib
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readOBJ(const string& filename)
{
	// clear anything that existed before
	_vertices.clear();
	_normals.clear();
  _triangles.clear();
  _texcoords.clear();

	// open up file
  FILE* file = fopen(filename.c_str(), "r");
	
	if (file == NULL)
	{
		cerr << "Can't read input file " << filename << endl;
		return false;
	}

  // make a type-dependent format string
  string format;
  if (sizeof(Real) == sizeof(double))
    format = string("%lf %lf %lf");
  else
    format = string("%f %f %f");

	// read through line by line
	int lineNumber = 0;
	bool faceSeen = false;
  while (true)
	{
		if (feof(file)) 
      break;

		char type[1024];
		lineNumber++;
    fscanf(file, "%s", type);

		if (feof(file) || ferror(file)) break;

    // see if it's a comment
    if (type[0] == '#')
    {
      fgets(type, 1024, file);
      continue;
    }

		// reading vertices
		if (strcmp(type, "v") == 0)
		{
      if (faceSeen)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Error! Vertex seen after face read begun!" << endl;
      }

      double d[3];
      fscanf(file, "%lf %lf %lf", &d[0], &d[1], &d[2]);
			VEC3F v(d[0], d[1], d[2]);
			_vertices.push_back(v);
		}

		// vertex normals
		if (strcmp(type, "vn") == 0)
		{
			VEC3F vn;
      fscanf(file, format.c_str(), &vn[0], &vn[1], &vn[2]);
			_normals.push_back(vn);
		}

		// vertex texcoords
		if (strcmp(type, "vt") == 0)
		{
			VEC3F vt;
      fscanf(file, format.c_str(), &vt[0], &vt[1], &vt[2]);
			_texcoords.push_back(vt);
    }

		// reading triangles
	  TRIANGLE f;
		if (type[0] == 'f')
		{
      char indices[256];
      faceSeen = true;
      
      for (int x = 0; x < 3; x++)
      {
        //int totalChars = fscanf(file, "%s", indices); 
        fscanf(file, "%s", indices); 
        if (feof(file) || ferror(file)) break;

        int vertexIndex = -1;
        //int texcoordIndex = -1;
        //int normalIndex = -1;
        char* texcoordExists = strchr(indices, '/');
        //char* normalExists   = strrchr(indices, '/');

        int vertexEnd = texcoordExists - indices;
        //int texcoordEnd = normalExists - indices;

        // extract the vertex index
        if (texcoordExists == NULL)
        {
          // if there is nothing but an index
          vertexIndex = atoi(indices);
        }
        else
        {
          // get the right substring
          char vertexString[256];
          strncpy(vertexString, indices, vertexEnd);

          // convert it to an int
          vertexIndex = atoi(vertexString);
        }

        /*
        // extract the texture index
        if (texcoordExists)
        {
          if (normalExists == NULL)
          {
            // extract to the end of the string
            char texcoordString[256];
            strncpy(texcoordString, &indices[vertexEnd + 1], totalChars - vertexEnd);
            
            // convert it to an int
            texcoordIndex = atoi(texcoordString);
          }
          else
          {
            // extract to the beginning of the normal index
            char texcoordString[256];
            strncpy(texcoordString, &indices[vertexEnd + 1], texcoordEnd - vertexEnd);
            
            // convert it to an int
            texcoordIndex = atoi(texcoordString);
          }
        }

        // extract the normal index
        if (normalExists)
          {
            // extract to the end of the string
            char normalString[256];
            strncpy(normalString, &indices[texcoordEnd + 1], totalChars - texcoordEnd);
              
            // convert it to an int
            normalIndex = atoi(normalString);
          }
          */

        // subtract one and store
        f.vertex(x) = &_vertices[vertexIndex - 1];
      }

      if (feof(file) || ferror(file)) break;

		  // store the triangle
		  _triangles.push_back(f);
		}
	}

  fclose(file);
  cout << filename.c_str() << " successfully loaded" << endl;
  cout << " Vertices: " << _vertices.size() << endl;
  cout << " Faces: " << _triangles.size() << endl;
  cout << " Normals: " << _normals.size() << endl;

  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

	return true;
}

//////////////////////////////////////////////////////////////////////
//get the range and center of the mesh
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::getBounds(VEC3F& maxVert, VEC3F& minVert, Real& maxLength){
	maxVert = minVert = _vertices[0];
	// get bounds
	for (unsigned int x = 0; x < _vertices.size(); x++)
		for (int y = 0; y < 3; y++)
	    {
	      if (_vertices[x][y] < minVert[y]) minVert[y] = _vertices[x][y];
	      if (_vertices[x][y] > maxVert[y]) maxVert[y] = _vertices[x][y];
	    }
	VEC3F lengths = maxVert - minVert;
	maxLength = lengths.maxElement();
}

//////////////////////////////////////////////////////////////////////
// Normalize mesh to 1x1x1 cube centered at (0,0,0)
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::normalize()
{
  VEC3F maxVert(_vertices[0]);
  VEC3F minVert(_vertices[0]);

  // get bounds
  for (unsigned int x = 0; x < _vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      if (_vertices[x][y] < minVert[y]) minVert[y] = _vertices[x][y];
      if (_vertices[x][y] > maxVert[y]) maxVert[y] = _vertices[x][y];
    }

  VEC3F half(0.5, 0.5, 0.5);
  VEC3F diff = maxVert - minVert;
  double maxDiff = diff.maxElement();
//cout << "mesh center " << (maxVert + minVert)/2.0 << " maxDiff " << maxDiff << endl;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x] -= minVert;
    _vertices[x] *= 1.0 / maxDiff;
    _vertices[x] -= half;
  }
}

//////////////////////////////////////////////////////////////////////
// fit the geometry inside a 1x1x1 box
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::normalize(Real padding)
{
  VEC3F currentMin(_vertices[0]);
  VEC3F currentMax(_vertices[0]);

	for (unsigned int i = 0; i < _vertices.size (); i++)
  {
		for (int j = 0; j < 3; j++)
		{
			currentMin[j] = min(currentMin[j],_vertices[i][j]);
			currentMax[j] = max(currentMax[j],_vertices[i][j]);
		}
    //cout << " vertex: " << _vertices[i] << endl;
  }
  cout << " Min: " << currentMin << endl;
  cout << " Max: " << currentMax << endl;

  VEC3F recenter = (currentMin + currentMax) * (Real)0.5f;

  // translate everything to the bounding box center
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] -= recenter;

  // find the maximum magnitude
  double maxVal = 0.0f;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    maxVal = (fabs(_vertices[x][0]) > maxVal) ? fabs(_vertices[x][0]) : maxVal;
    maxVal = (fabs(_vertices[x][1]) > maxVal) ? fabs(_vertices[x][1]) : maxVal;
    maxVal = (fabs(_vertices[x][2]) > maxVal) ? fabs(_vertices[x][2]) : maxVal;
  }
  cout << " Max value found: " << maxVal << endl;

  // add a little padding for the last grid cell
  //Real scale = 0.5 - (1.0 / 16.0);
  Real scale = 0.5 - padding;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] *= scale / maxVal;

  // translate everything to 0.5, 0.5, 0.5
  VEC3F half(0.5, 0.5, 0.5);
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] += half;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::inside(const VEC3F& point)
{
  unsigned int j;
  //float vec[3];
  VEC3F vec;

  //float pt[3];
  VEC3F pt;
  pt[0] = point[0];
  pt[1] = point[1];
  pt[2] = point[2];

  // try all cardinal directions, take majority
  vec[0] = 1.0f; vec[1] = 1.0f; vec[2] = 0.0f;
  int xPos = 0;
  for (j = 0; j < _triangles.size(); j++)
    if (intersectTriangle(pt, vec, j))
      xPos++;

  // y positive
  vec[0] = 0.0f; vec[1] = 1.0f; vec[2] = 0.0f;
  int yPos = 0;
  for (j = 0; j < _triangles.size(); j++)
    if (intersectTriangle(pt, vec, j))
      yPos++;

  // z positive
  vec[0] = 0.0f; vec[1] = 0.0f; vec[2] = 1.0f;
  int zPos = 0;
  for (j = 0; j < _triangles.size(); j++)
    if (intersectTriangle(pt, vec, j))
      zPos++;

  // x negative
  vec[0] = -1.0f; vec[1] = 0.0f; vec[2] = 0.0f;
  int xNeg = 0;
  for (j = 0; j < _triangles.size(); j++)
    if (intersectTriangle(pt, vec, j))
      xNeg++;
 
  // y negative
  vec[0] = 0.0f; vec[1] = -1.0f; vec[2] = 0.0f;
  int yNeg = 0;
  for (j = 0; j < _triangles.size(); j++)
    if (intersectTriangle(pt, vec, j))
      yNeg++;
  
  // z negative
  vec[0] = 0.0f; vec[1] = 0.0f; vec[2] = -1.0f;
  int zNeg = 0;
  for (j = 0; j < _triangles.size(); j++)
    if (intersectTriangle(pt, vec, j))
      zNeg++;
  
  float sum = (xPos % 2 + yPos % 2 + zPos % 2 + xNeg % 2 + yNeg % 2 + zNeg % 2);
 
  // take a majority vote 
  return (sum > 3);
}

//////////////////////////////////////////////////////////////////////
// Slow, but correct
//////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::distance(const VEC3F& point)
{
  if (_triangles.size() == 0)
    return 0;

  // do an exhaustive search
  Real final = pointFaceDistanceSq(_triangles[0], point);
  for (unsigned int x = 1; x < _triangles.size(); x++)
  {
    Real currentDistance = pointFaceDistanceSq(_triangles[x], point);
    if (currentDistance < final || x == 0)
      final = currentDistance;
  }

  return sqrt(final);
}

//////////////////////////////////////////////////////////////////////
// Slow, but correct
//////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::signedDistance(const VEC3F& point)
{
  Real final = distance(point);

  return (inside(point)) ? -final : final;
}

//////////////////////////////////////////////////////////////////////
// initialize a signed distance field by only setting values along
// the surface -- replace the slow default version with
// one that only looks at grid cells in a triangle's neighborhood
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::initializeSignedDistanceField(FIELD_3D& field)
{
  field = field.outside();

  Real diagonal = sqrt(2.);
  Real dMax = 1.0 / field.maxRes();

  // for each triangle
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    // get the bounding box
    VEC3F mins, maxs;
    _triangles[i].boundingBox(mins, maxs);

    // get the integer ranges to iterate over
    VEC3I iMins, iMaxs;
    field.boundingBoxIndices(mins, maxs, iMins, iMaxs);

    for (int z = iMins[2]; z < iMaxs[2]; z++)
      for (int y = iMins[1]; y < iMaxs[1]; y++)
        for (int x = iMins[0]; x < iMaxs[0]; x++)
        {
          int index = x + y * field.xRes() + z * field.slabSize();

          // compute distance at those points
          VEC3F point = field.cellCenter(x,y,z);
          Real distance = pointFaceDistanceSq(_triangles[i], point);
          distance = sqrt(distance);
          distance *= 1.0 / dMax;

          if (distance > diagonal) continue;

          // see if it's closer than what's been seen before
          if (distance < fabs(field[index]))
          {
            Real sign = (inside(point)) ? -1.0 : 1.0;
            field[index] = sign * distance;
          }
        }
  }
}

//////////////////////////////////////////////////////////////////////
// ray-triangle test from
// Tomas Möller and Ben Trumbore. Fast, minimum storage ray-triangle intersection. 
// Journal of graphics tools, 2(1):21-28, 1997
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::intersectTriangle(const VEC3F& orig, const VEC3F& dir, int faceIndex)
{
   VEC3F edge1, edge2, tvec, pvec, qvec;
   float det,inv_det;
   float t,u,v;

   const VEC3F& vert0 = *(_triangles[faceIndex].vertex(0));
   const VEC3F& vert1 = *(_triangles[faceIndex].vertex(1));
   const VEC3F& vert2 = *(_triangles[faceIndex].vertex(2));

   // find vectors for two edges sharing vert0
   edge1 = vert1 - vert0;
   edge2 = vert2 - vert0;

   // begin calculating determinant - also used to calculate U parameter
   pvec = cross(dir, edge2);

   // if determinant is near zero, ray lies in plane of triangle
   det = edge1 * pvec;

   float eps = 1e-8;
   if (det > -eps && det < eps)
     return false;
   inv_det = 1.0 / det;

   // calculate distance from vert0 to ray origin
   tvec = orig - vert0;

   // calculate U parameter and test bounds
   u = (tvec * pvec) * inv_det;
   if (u < 0.0 || u > 1.0)
     return false;

   // prepare to test V parameter
   qvec = cross(tvec, edge1);

   // calculate V parameter and test bounds
   v = (dir * qvec) * inv_det;
   if (v < 0.0 || u + v > 1.0)
     return false;

   t = (edge2 * qvec) * inv_det;

   if (t < 0.0f) return false;
   return true;
}

//////////////////////////////////////////////////////////////////////////////
// This is a modified version of the Wild Magic DistVector3Triangle3 class.
// This is what is used to generate the distance grid.
//
// The license info is as below:
//
// Wild Magic Source Code
// David Eberly
// http://www.geometrictools.com
// Copyright (c) 1998-2008
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
//
// Version: 4.0.1 (2007/05/06)
//----------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::pointFaceDistanceSq(TRIANGLE& f, const VEC3F& point)
{
  VEC3F v0 = *(f.vertex(0));
  VEC3F v1 = *(f.vertex(1));
  VEC3F v2 = *(f.vertex(2));
  VEC3F kDiff = v0 - point;
  VEC3F kEdge0 = v1 - v0;
  VEC3F kEdge1 = v2 - v0;
  Real fA00 = norm2(kEdge0);
  Real fA01 = kEdge0 * kEdge1;
  Real fA11 = norm2(kEdge1);
  Real fB0 = kDiff*kEdge0;
  Real fB1 = kDiff*kEdge1;
  Real fC = norm2(kDiff);
  Real fDet = fabs(fA00*fA11-fA01*fA01);
  Real fS = fA01*fB1-fA11*fB0;
  Real fT = fA01*fB0-fA00*fB1;
  Real fSqrDistance;

  if (fS + fT <= fDet) {
      if (fS < (Real)0.0) {
          if (fT < (Real)0.0)  // region 4
          {
              if (fB0 < (Real)0.0) {
                  fT = (Real)0.0;
                  if (-fB0 >= fA00) {
                      fS = (Real)1.0;
                      fSqrDistance = fA00+((Real)2.0)*fB0+fC;
                  }
                  else {
                      fS = -fB0/fA00;
                      fSqrDistance = fB0*fS+fC;
                  }
              }
              else {
                  fS = (Real)0.0;
                  if (fB1 >= (Real)0.0) {
                      fT = (Real)0.0;
                      fSqrDistance = fC;
                  }
                  else if (-fB1 >= fA11) {
                      fT = (Real)1.0;
                      fSqrDistance = fA11+((Real)2.0)*fB1+fC;
                  }
                  else {
                      fT = -fB1/fA11;
                      fSqrDistance = fB1*fT+fC;
                  }
              }
          }
          else  // region 3
          {
              fS = (Real)0.0;
              if (fB1 >= (Real)0.0) {
                  fT = (Real)0.0;
                  fSqrDistance = fC;
              }
              else if (-fB1 >= fA11) {
                  fT = (Real)1.0;
                  fSqrDistance = fA11+((Real)2.0)*fB1+fC;
              }
              else {
                  fT = -fB1/fA11;
                  fSqrDistance = fB1*fT+fC;
              }
          }
      }
      else if (fT < (Real)0.0)  // region 5
      {
          fT = (Real)0.0;
          if (fB0 >= (Real)0.0) {
              fS = (Real)0.0;
              fSqrDistance = fC;
          }
          else if (-fB0 >= fA00) {
              fS = (Real)1.0;
              fSqrDistance = fA00+((Real)2.0)*fB0+fC;
          }
          else {
              fS = -fB0/fA00;
              fSqrDistance = fB0*fS+fC;
          }
      }
      else  // region 0
      {
          // minimum at interior point
          Real fInvDet = ((Real)1.0)/fDet;
          fS *= fInvDet;
          fT *= fInvDet;
          fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
              fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
      }
  }
  else {
      Real fTmp0, fTmp1, fNumer, fDenom;

      if (fS < (Real)0.0)  // region 2
      {
          fTmp0 = fA01 + fB0;
          fTmp1 = fA11 + fB1;
          if (fTmp1 > fTmp0) {
              fNumer = fTmp1 - fTmp0;
              fDenom = fA00-2.0f*fA01+fA11;
              if (fNumer >= fDenom) {
                  fS = (Real)1.0;
                  fT = (Real)0.0;
                  fSqrDistance = fA00+((Real)2.0)*fB0+fC;
              }
              else {
                  fS = fNumer/fDenom;
                  fT = (Real)1.0 - fS;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                      fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
              }
          }
          else {
              fS = (Real)0.0;
              if (fTmp1 <= (Real)0.0) {
                  fT = (Real)1.0;
                  fSqrDistance = fA11+((Real)2.0)*fB1+fC;
              }
              else if (fB1 >= (Real)0.0) {
                  fT = (Real)0.0;
                  fSqrDistance = fC;
              }
              else {
                  fT = -fB1/fA11;
                  fSqrDistance = fB1*fT+fC;
              }
          }
      }
      else if (fT < (Real)0.0)  // region 6
      {
          fTmp0 = fA01 + fB1;
          fTmp1 = fA00 + fB0;
          if (fTmp1 > fTmp0) {
              fNumer = fTmp1 - fTmp0;
              fDenom = fA00-((Real)2.0)*fA01+fA11;
              if (fNumer >= fDenom) {
                  fT = (Real)1.0;
                  fS = (Real)0.0;
                  fSqrDistance = fA11+((Real)2.0)*fB1+fC;
              }
              else {
                  fT = fNumer/fDenom;
                  fS = (Real)1.0 - fT;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
                      fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
              }
          }
          else {
              fT = (Real)0.0;
              if (fTmp1 <= (Real)0.0) {
                  fS = (Real)1.0;
                  fSqrDistance = fA00+((Real)2.0)*fB0+fC;
              }
              else if (fB0 >= (Real)0.0) {
                  fS = (Real)0.0;
                  fSqrDistance = fC;
              }
              else {
                  fS = -fB0/fA00;
                  fSqrDistance = fB0*fS+fC;
              }
          }
      }
      else  // region 1
      {
          fNumer = fA11 + fB1 - fA01 - fB0;
          if (fNumer <= (Real)0.0) {
              fS = (Real)0.0;
              fT = (Real)1.0;
              fSqrDistance = fA11+((Real)2.0)*fB1+fC;
          }
          else {
              fDenom = fA00-2.0f*fA01+fA11;
              if (fNumer >= fDenom) {
                  fS = (Real)1.0;
                  fT = (Real)0.0;
                  fSqrDistance = fA00+((Real)2.0)*fB0+fC;
              }
              else {
                  fS = fNumer/fDenom;
                  fT = (Real)1.0 - fS;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
                      fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
              }
          }
      }
  }

  // account for numerical round-off error
  if (fSqrDistance < (Real)0.0) {
      fSqrDistance = (Real)0.0;
  }
  return fSqrDistance;
}

//////////////////////////////////////////////////////////////////////////////
// initialize an extension scalar to a field
//////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::initializeExtensionScalars(const FIELD_3D& distance, FIELD_3D& scalar)
{
  cout << " Initializing extension scalars ... "; flush(cout);

  int xRes = scalar.xRes();
  int yRes = scalar.yRes();
  int zRes = scalar.zRes();
  int slabSize = scalar.slabSize();
  //int maxRes = scalar.maxRes();

  // clear all old values
  scalar = 0;

  // count how many scalars were added to various cells
  // <cell index, how many scalars got added>
  map<int, Real> scalarCounts;

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    VEC3F gridCoords;
    VEC3F& point = _vertices[x];

    Real interp = point[0];
    if (_texcoords.size() == _vertices.size())
      interp = _texcoords[x][0];
  
    // get the lower corner position
    VEC3F corner = distance.center() - (Real)0.5 * distance.lengths();
    VEC3F dxs = distance.dxs();
    corner += (Real)0.5 * dxs;
  
    // recenter position
    gridCoords = point;
    gridCoords -= corner;
    gridCoords[0] *= 1.0 / dxs[0];
    gridCoords[1] *= 1.0 / dxs[1];
    gridCoords[2] *= 1.0 / dxs[2];

    // find which direction is farthest from an int, then interpolate
    // in that direction
    Real floors[3];
    Real ceils[3];
    Real farthest = 0.0f;
    int farthestInt = 0;
    for (int i = 0; i < 3; i++)
    {
      floors[i] = gridCoords[i] - floor(gridCoords[i]);
      ceils[i] = ceil(gridCoords[i]) - gridCoords[i];

      if (floors[i] > 0.5)
      {
        if (ceils[i] > farthest)
        {
          farthest = ceils[i];
          farthestInt = i;
        }
      }
      else
        if (floors[i] > farthest)
        {
          farthest = floors[i];
          farthestInt = i;
        }
    }

    int index = scalar.cellIndex(point);
    scalarCounts[index] = 0;
    scalar[index] = 0;

    // if it's right on top of a point, just set that point
    if (farthest < 1e-5)
    {
      int index = scalar.cellIndex(point);

      // DEBUG: SCALAR TO INTERPOLATE IS ACCESSED HERE
      scalar[index] = interp;
      scalarCounts[index] = 1.0f;
      continue;
    }

    // interpolate in z direction
    if (farthestInt == 2)
    {
      scalarCounts[index] += floors[2];

      // DEBUG: SCALAR TO INTERPOLATE IS ACCESSED HERE
      scalar[index] += floors[2] * interp;

      // get the z plus neighbor
      int neighbor = index + slabSize;

      // if there is no neighbor, or there is no zero crossing neighbor,
      // get the z minus neighbor
      if (!(neighbor > zRes - 1) || distance[neighbor] * distance[index] >= 0.0f)
        neighbor = index - slabSize;
      
      // make sure these is a zero crossing here
      if (distance[neighbor] * distance[index] >= 0.0f)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " z neighbor missing in extension velocities" << endl;
        cout << " index: " << index << endl;
        distance.printNeighborhood(index);
        cout << " vertex: " << x << endl;
        cout << " grid coords: " << gridCoords << endl;
        cout << " Point: " << point << endl;
        exit(0);
      }

      if (scalarCounts.find(neighbor) == scalarCounts.end())
      {
        scalarCounts[neighbor] = 0.0f;
        scalar[neighbor] = 0.0f;
      }

      scalarCounts[neighbor] += ceils[2];
      
      // DEBUG: SCALAR TO INTERPOLATE IS ACCESSED HERE
      scalar[neighbor] += ceils[2] * interp;
    }

    // interpolate in y direction
    if (farthestInt == 1)
    {
      scalarCounts[index] += floors[1];

      // DEBUG: SCALAR TO INTERPOLATE IS ACCESSED HERE
      scalar[index] += floors[1] * interp;

      // get the y plus neighbor
      int neighbor = index + xRes;

      // if there is no neighbor, or there is no zero crossing neighbor,
      // get the y minus neighbor
      if (!(neighbor > yRes - 1) || distance[neighbor] * distance[index] >= 0.0f)
        neighbor = index - xRes;
      
      // make sure these is a zero crossing here
      if (distance[neighbor] * distance[index] >= 0.0f)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << "y neighbor missing in extension velocities" << endl;
        cout << " index: " << index << endl;
        distance.printNeighborhood(index);
        cout << " vertex: " << x << endl;
        cout << " grid coords: " << gridCoords << endl;
        cout << " Point: " << point << endl;
        exit(0);
      }

      if (scalarCounts.find(neighbor) == scalarCounts.end())
      {
        scalarCounts[neighbor] = 0.0f;
        scalar[neighbor] = 0.0f;
      }

      scalarCounts[neighbor] += ceils[1];
      
      // DEBUG: SCALAR TO INTERPOLATE IS ACCESSED HERE
      scalar[neighbor] += ceils[1] * interp;
    }

    // interpolate in x direction
    if (farthestInt == 0)
    {
      scalarCounts[index] += floors[0];

      // DEBUG: SCALAR TO INTERPOLATE IS ACCESSED HERE
      scalar[index] += floors[0] * interp;

      // get the y plus neighbor
      int neighbor = index + 1;

      // if there is no neighbor, or there is no zero crossing neighbor,
      // get the y minus neighbor
      if (!(neighbor > xRes - 1) || distance[neighbor] * distance[index] >= 0.0f)
        neighbor = index - 1;
      
      // make sure these is a zero crossing here
      if (distance[neighbor] * distance[index] >= 0.0f)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << "x neighbor missing in extension velocities" << endl;
        cout << " index: " << index << endl;
        distance.printNeighborhood(index);
        cout << " vertex: " << x << endl;
        cout << " grid coords: " << gridCoords << endl;
        cout << " Point: " << point << endl;
        exit(0);
      }

      if (scalarCounts.find(neighbor) == scalarCounts.end())
      {
        scalarCounts[neighbor] = 0.0f;
        scalar[neighbor] = 0.0f;
      }

      scalarCounts[neighbor] += ceils[0];
      
      // DEBUG: SCALAR TO INTERPOLATE IS ACCESSED HERE
      scalar[neighbor] += ceils[0] * interp;
    }
  }

  map<int, Real>::iterator iter;
  for (iter = scalarCounts.begin(); iter != scalarCounts.end(); iter++)
  {
    int index = iter->first;
    Real count = iter->second;
    scalar[index] *= 1.0 / count;
  }
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// center using vertex means
///////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::vertexMean() const
{
  VEC3F sum;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    sum += _vertices[x];

  return sum * (Real)(1.0 / _vertices.size());
}

///////////////////////////////////////////////////////////////////////
// assign a texture coordinates with cylindrical projection
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setCylindricalTexcoords(float bottomCap, float topCap)
{
  float interval = topCap - bottomCap;
  //float capMean = (topCap + bottomCap) * 0.5;
  _texcoords.clear();

  // use the vertex mean as the center
  VEC3F center = vertexMean();

  VEC3F uvMins;
  VEC3F uvMaxs;

  //for (unsigned int x = 0; x < _vertices.size(); x++)
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    Real phi[3];
    Real theta[3];

    for (int y = 0; y < 3; y++)
    {
      VEC3F vertex = _vertices[vertexIndex(x,y)];
      VEC3F centered = vertex - center;

      theta[y] = (centered[1] - bottomCap) / interval;
      phi[y]   = atan2(centered[0], centered[2]);

      phi[y] = (phi[y] < 0) ? phi[y] + 2.0 * M_PI : phi[y];
    }

    // account for wraparound --
    // if one phi is big and the others small, make them all big
    int lessThan = 0;
    int greaterThan = 0;
    for (int y = 0; y < 3; y++)
    {
      phi[y] *= 1.0 / (2.0 * M_PI);

      if (phi[y] < 0.75 && phi[y] > 0.25)
        continue;

      if (phi[y] > 0.75)
        greaterThan++;

      if (phi[y] < 0.25)
        lessThan++;
    }

    if (lessThan > 0 && greaterThan > 0)
    {
      //clamp the less thans
      for (int y = 0; y < 3; y++)
      {
        if (phi[y] < 0.25)
          phi[y] += 1.0;
      }
    }

    for (int y = 0; y < 3; y++)
    {
      VEC3F uv;
      uv[0] = phi[y];
      uv[1] = theta[y];

      _texcoords.push_back(uv);

      if (uv[0] < uvMins[0]) uvMins[0] = uv[0];
      if (uv[1] < uvMins[1]) uvMins[1] = uv[1];

      if (uv[0] > uvMaxs[0]) uvMaxs[0] = uv[0];
      if (uv[1] > uvMaxs[1]) uvMaxs[1] = uv[1];
    }
  }
}

///////////////////////////////////////////////////////////////////////
// assign a texture coordinates with spherical projection
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setSphericalTexcoords()
{
  _texcoords.clear();

  // use the vertex mean as the center
  VEC3F center = vertexMean();

  VEC3F uvMins;
  VEC3F uvMaxs;

  //for (unsigned int x = 0; x < _vertices.size(); x++)
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    Real phi[3];
    Real theta[3];

    for (int y = 0; y < 3; y++)
    {
      VEC3F vertex = _vertices[vertexIndex(x,y)];
      VEC3F centered = vertex - center;

      centered.normalize();
      theta[y] = acos(centered[2]);
      phi[y]   = atan2(centered[1], centered[0]);

      phi[y] = (phi[y] < 0) ? phi[y] + 2.0 * M_PI : phi[y];
    }

    // account for wraparound --
    // if one phi is big and the others small, make them all big
    int lessThan = 0;
    int greaterThan = 0;
    for (int y = 0; y < 3; y++)
    {
      phi[y] *= 1.0 / (2.0 * M_PI);

      if (phi[y] < 0.75 && phi[y] > 0.25)
        continue;

      if (phi[y] > 0.75)
        greaterThan++;

      if (phi[y] < 0.25)
        lessThan++;
    }

    if (lessThan > 0 && greaterThan > 0)
    {
      //clamp the less thans
      for (int y = 0; y < 3; y++)
      {
        if (phi[y] < 0.25)
          phi[y] += 1.0;
      }
    }

    for (int y = 0; y < 3; y++)
    {
      VEC3F uv;
      uv[0] = phi[y];
      uv[1] = (M_PI - theta[y]) / M_PI;

      _texcoords.push_back(uv);

      if (uv[0] < uvMins[0]) uvMins[0] = uv[0];
      if (uv[1] < uvMins[1]) uvMins[1] = uv[1];

      if (uv[0] > uvMaxs[0]) uvMaxs[0] = uv[0];
      if (uv[1] > uvMaxs[1]) uvMaxs[1] = uv[1];
    }
  }
}

///////////////////////////////////////////////////////////////////////
// assign a texture coordinates with swirled spherical projection,
// just to get something more non-trivial for debugging purposes
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setSwirledTexcoords()
{
  _texcoords.clear();

  // use the vertex mean as the center
  VEC3F center = vertexMean();

  VEC3F uvMins;
  VEC3F uvMaxs;

  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    Real phi[3];
    Real theta[3];

    for (int y = 0; y < 3; y++)
    {
      VEC3F vertex = _vertices[vertexIndex(x,y)];
      VEC3F centered = vertex - center;

      centered.normalize();
      theta[y] = acos(centered[2]);
      phi[y]   = atan2(centered[1], centered[0]);

      phi[y] += 8 * theta[y];
      phi[y] = (phi[y] < 0) ? phi[y] + 2.0 * M_PI : phi[y];
    }

    // account for wraparound --
    // if one phi is big and the others small, make them all big
    int lessThan = 0;
    int greaterThan = 0;
    for (int y = 0; y < 3; y++)
    {
      phi[y] *= 1.0 / (2.0 * M_PI);

      if (phi[y] < 0.75 && phi[y] > 0.25)
        continue;

      if (phi[y] > 0.75)
        greaterThan++;

      if (phi[y] < 0.25)
        lessThan++;
    }

    if (lessThan > 0 && greaterThan > 0)
    {
      //clamp the less thans
      for (int y = 0; y < 3; y++)
      {
        if (phi[y] < 0.25)
          phi[y] += 1.0;
      }
    }

    for (int y = 0; y < 3; y++)
    {
      VEC3F uv;
      uv[0] = phi[y];
      uv[1] = (M_PI - theta[y]) / M_PI;

      _texcoords.push_back(uv);

      if (uv[0] < uvMins[0]) uvMins[0] = uv[0];
      if (uv[1] < uvMins[1]) uvMins[1] = uv[1];

      if (uv[0] > uvMaxs[0]) uvMaxs[0] = uv[0];
      if (uv[1] > uvMaxs[1]) uvMaxs[1] = uv[1];
    }
  }
}

///////////////////////////////////////////////////////////////////////
// write out a BREP file for Sean Mauch's closest point 
// transform (CPT) code
//
// from CPT code driver.cc:
//
//  The b-rep file contains a triangle mesh.  (Specifically, an indexed triangle
//  face set.)  The first two fields give the number of vertices and the number
//  of faces in the mesh.  This is followed by the coordinates of the vertices.
//  Next the faces are enumerated.  A face is specified by the indices of
//  three vertices.  (These indices are in the range [0..num_vertices-1].)
//
//  \verbatim
//  num_vertices
//  num_faces
//  vertex_0_x    vertex_0_y	vertex_0_z
//  vertex_1_x    vertex_1_y	vertex_1_z
//  ...
//  face_0_index_0	face_0_index_1		face_0_index_2
//  face_1_index_0	face_1_index_1		face_1_index_2
//  ...
//  \endverbatim
//
// Note that the CPT code will also need a geom file:
//
//  The geometry file contains information about the computational domain,
//  the size of the grids and options for the CPT.
//
//  \verbatim
//  xmin   ymin   zmin   xmax   ymax   zmax
//  num_grid_points_x   num_grid_points_y   num_grid_points_z
//  max_distance
//  oriented
//  \endverbatim
//
//  Description of fields:
//  - \e xmin, etc. describe the Cartesian domain spanned by the grid.
//  - \e num_grid_points_* are the number of grid points in each direction.
//  - The signed distance and closest point will be correctly computed up to 
//    \e max_distance from the surface.  A value of 0 for \e max_distance 
//    indicates that \e max_distance will be set equal to the diagonal 
//    length of the Cartesian domain.
//  - The distance for any points farther away than max_distance is set 
//    to \c std::numeric_limits<number_type>::max().  If the -l option is
//    given, the far away distances are set to +- max_distance.  Each
//    component of the gradient of the distance and closest point for
//    far away grid points are set to std::numeric_limits<number_type>::max().
//    The closest face for far away grid points is set to -1.
//  - \e oriented is a boolean value, (0 or 1), indicating if the surface 
//    is oriented.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeBRep(const char* filename)
{
  FILE* file = fopen(filename, "w");

  if (file == NULL)
  {
    cout << " Couldn't open filename " << filename << "!!!" << endl;
    fclose(file);
    exit(0);
  }

  cout << " Writing brep file " << filename << " ... "; flush(cout);

  // write dims
  fprintf(file, "%i\n", (int)_vertices.size());
  fprintf(file, "%i\n", (int)_triangles.size());

  // write vertex positions
  for (unsigned int x = 0; x < _vertices.size(); x++)
    fprintf(file, "%f %f %f\n", (double)_vertices[x][0], (double)_vertices[x][1], (double)_vertices[x][2]);

  // write face indices
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    int indices[3];
    indices[0] = _vertexIndices[_triangles[x].vertex(0)];
    indices[1] = _vertexIndices[_triangles[x].vertex(1)];
    indices[2] = _vertexIndices[_triangles[x].vertex(2)];
    fprintf(file, "%i %i %i\n", indices[0], indices[1], indices[2]);
  }

  fclose(file);

  cout << " done." << endl;
}

///////////////////////////////////////////////////////////////////////
// return a bounding box
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  mins = _vertices[0];
  maxs = _vertices[0];

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    if (_vertices[x][0] < mins[0]) mins[0] = _vertices[x][0];
    if (_vertices[x][1] < mins[1]) mins[1] = _vertices[x][1];
    if (_vertices[x][2] < mins[2]) mins[2] = _vertices[x][2];

    if (_vertices[x][0] > maxs[0]) maxs[0] = _vertices[x][0];
    if (_vertices[x][1] > maxs[1]) maxs[1] = _vertices[x][1];
    if (_vertices[x][2] > maxs[2]) maxs[2] = _vertices[x][2];
  }
}

///////////////////////////////////////////////////////////////////////
// normalize any soid texture that has been previously generated
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::normalizeSolidTexture()
{
  Real maxFound = _solidTextureMeans[0];
  for (unsigned int x = 1; x < _solidTextureMeans.size(); x++)
  {
    if (_solidTextureMeans[x] > maxFound)
      maxFound = _solidTextureMeans[x];
  }
  for (unsigned int x = 1; x < _solidTextureMeans.size(); x++)
    _solidTextureMeans[x] *= 1.0 / maxFound;
}

///////////////////////////////////////////////////////////////////////
// get the center coordinate of the triangle with the maximum 
// texture value
///////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::maxTextureCenter()
{
  Real maxFound = _solidTextureMeans[0];
  int whichTriangle = 0;
  for (unsigned int x = 1; x < _solidTextureMeans.size(); x++)
    if (_solidTextureMeans[x] > maxFound)
    {
      maxFound = _solidTextureMeans[x];
      whichTriangle = x;
    }

  // get the vertices
  vector<VEC3F> verts;
  verts.push_back(*(_triangles[whichTriangle].vertices()[0]));
  verts.push_back(*(_triangles[whichTriangle].vertices()[1]));
  verts.push_back(*(_triangles[whichTriangle].vertices()[2]));
  VEC3F vertMean = verts[0] + verts[1] + verts[2];
  vertMean *= 1.0 / 3.0;

  return vertMean;
}

///////////////////////////////////////////////////////////////////////
// write out LuxRender geometry
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeLuxRender(const char* filename, const char* path)
{
  assert(_solidTextures.size() > 0);

  string fullfile = string(path) + string(filename);
  FILE* file = fopen(fullfile.c_str(), "w");

  cout << " Dumping " << _triangles.size() << " JPG files ... "; flush(cout);
  int size = _triangles.size();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < size; x++)
  {
    char buffer[256];
    sprintf(buffer, "%i", x);
    string number(buffer);

    // get its texture
    FIELD_2D& texture = _solidTextures[x];

    // dump it to a jpg
    string jpgName = string(path) + string("triangle.") + number + string(".jpg");
    texture.writeJPG(jpgName);
  }
  cout << " done. " << endl;
   
  cout << " Writing LuxRender file ... "; flush(cout); 
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    char buffer[256];
    sprintf(buffer, "%i", x);
    string number(buffer);
    string jpgName = string(path) + string("triangle.") + number + string(".jpg");

    // create the texture in LuxRender:
    //
    //  Texture "MyImage" "color" "imagemap"
    //    "string filename" ["res/test.jpg"]
    //
    string textureName = string("Texture_") + number;
    fprintf(file, "Texture \"%s\" \"color\" \"imagemap\"\n", textureName.c_str());
    fprintf(file, "\"string filename\" [\"%s\"]\n", jpgName.c_str());
    fprintf(file, "\"string filtertype\" [\"bilinear\"]\n");
    //fprintf(file, "\"float uscale\" [1]\n");
    //fprintf(file, "\"float vscale\" [1]\n");
    
    // make the named material
    //
    //  MakeNamedMaterial "MyMaterial"
    //    "texture Kd" ["MyImage"]
    //
    string materialName = string("Material_") + number;
    fprintf(file, "MakeNamedMaterial \"%s\"\n", materialName.c_str());
    fprintf(file, "\"texture Kd\" [\"%s\"]\n", textureName.c_str());
    
    // draw the triangle
    //
    // AttributeBegin 
    //  NamedMaterial "MyMaterial"
    //  Shape "trianglemesh"
    //    "integer indices" [0 1 2]
    //    "point P" [-0.5 0.5 0  0.5 -0.5 0  -0.5 -0.5 0]
    //    "float uv" [ 1 1  0 0  1 0]
    // AttributeEnd
    //
    fprintf(file, "AttributeBegin\n");
    fprintf(file, " NamedMaterial \"%s\"\n", materialName.c_str());
    fprintf(file, " Shape \"trianglemesh\"\n");
    fprintf(file, " \"integer indices\" [0 1 2]\n");
    fprintf(file, " \"point P\" [");

    TRIANGLE& triangle = _triangles[x];
    for (int y = 0; y < 3; y++)
    {
      VEC3F& vertex = *(triangle.vertex(y));
      fprintf(file, " %f %f %f", (double)vertex[0], (double)vertex[1], (double)vertex[2]);
    }

    fprintf(file, "]\n");
    //fprintf(file, " \"float uv\" [ 1 1  1 0  0 0]\n");
    //fprintf(file, " \"float uv\" [ 1 1  0 0  1 0]\n");
    fprintf(file, " \"float uv\" [");
    for (int y = 0; y < 3; y++)
      fprintf(file, "%f %f ", (double)_texcoords[3 * x + y][0], (double)_texcoords[3 * x + y][1]);
    fprintf(file, "]\n");
    fprintf(file, "AttributeEnd\n");
  }
  cout << " done." << endl;

  fclose(file);
}

///////////////////////////////////////////////////////////////////////
// write out PBRT geometry
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writePBRT(const char* filename)
{
  FILE* file = fopen(filename, "w");

  // read in the vertices
  fprintf(file, " Shape \"trianglemesh\"\n");
  fprintf(file, "\"point P\" [\n");
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    VEC3F vertex = _vertices[x];
    fprintf(file, "%f %f %f\n", (double)vertex[0], (double)vertex[1], (double)vertex[2]);
  }
  fprintf(file, "]\n");

  fprintf(file, "\"integer indices\" [\n");
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    int indices[3];
    for (int y = 0; y < 3; y++)
      indices[y] = _vertexIndices[_triangles[x].vertex(y)];
    fprintf(file, "%i %i %i\n", indices[0], indices[1], indices[2]);
  }
  fprintf(file, "]\n");

  fclose(file);
}


///////////////////////////////////////////////////////////////////////
// write out bobj gz
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeBobjGz(const char* filename)
{
	gzFile gzf;
	gzf = gzopen(filename, "wb1");
	if (gzf == NULL)
	{
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
		cout << " TRIANGLE_MESH::writeBobjGz failed! " << endl;
		cout << " Could not open file " << filename << endl;
		exit(0);
	}

	cout << " Writing file " << filename << " ... "; flush(cout);

	// read in the scalars
	//gzread(file, (void*)&_xRes, sizeof(int));


	// write to file
	int numVerts;
	if(sizeof(numVerts)!=4) { assert(false); exit(1); } // sanity check
	numVerts = _vertices.size();

	int numTris = _triangles.size();

	gzwrite(gzf, &numVerts, sizeof(numVerts));
	for(int i=0; i<numVerts; i++) {
		for(int j=0; j<3; j++) {
			float vertp = _vertices[i][j];
			gzwrite(gzf, &vertp, sizeof(vertp)); }
	}

	// normals dont matter , ignore!
	gzwrite(gzf, &numVerts, sizeof(numVerts));
	for(int i=0; i<numVerts; i++) {
		for(int j=0; j<3; j++) {
			float normp = 0.f;
			gzwrite(gzf, &normp, sizeof(normp)); }
	}

	gzwrite(gzf, &numTris, sizeof(numTris));
	for(int i=0; i<numTris; i++) {
		for (int y = 0; y < 3; y++)
		{
			int triIndex = _vertexIndices[_triangles[i].vertex(y)];
			gzwrite(gzf, &triIndex, sizeof(triIndex)); 
		}
	}

	// done
	gzclose( gzf );
}


///////////////////////////////////////////////////////////////////////
// distinguish where to load physBam files from - move to utility file?
///////////////////////////////////////////////////////////////////////
bool checkDdfLoadEnvVar() 
{
	bool haveDdfFiles = false;
	if(getenv("LOADDDF")) 
	{
		if(atoi(getenv("LOADDDF"))>=1) 
		{
			haveDdfFiles = true;
		}
	}
	return haveDdfFiles;
}

////////////////////////////////////////////////////////////////////////////
// for debugging, set the triangle mesh to a single tet (used to test
// Loop subdivision)
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setToTet()
{
  _vertices.clear();
  _normals.clear();
  _texcoords.clear();
  _triangles.clear();

  VEC3F vertices[4];
  vertices[0][0] = 0.5;
  vertices[0][1] = sqrt(3.0) / 6;
  vertices[0][2] = sqrt(2.0 / 3.0);

  vertices[1][0] = 0.5;
  vertices[1][1] = sqrt(3.0) / 2.0;
  vertices[1][2] = 0;

  vertices[2][0] = 1.0;
  vertices[2][1] = 0;
  vertices[2][2] = 0;

  vertices[3][0] = 0;
  vertices[3][1] = 0;
  vertices[3][2] = 0;

  _vertices.push_back(vertices[0]);
  _vertices.push_back(vertices[1]);
  _vertices.push_back(vertices[2]);
  _vertices.push_back(vertices[3]);

  _triangles.push_back(TRIANGLE(&_vertices[0], & _vertices[1], &_vertices[3]));
  _triangles.push_back(TRIANGLE(&_vertices[0], & _vertices[2], &_vertices[1]));
  _triangles.push_back(TRIANGLE(&_vertices[0], & _vertices[3], &_vertices[2]));
  _triangles.push_back(TRIANGLE(&_vertices[1], & _vertices[2], &_vertices[3]));

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;
}

////////////////////////////////////////////////////////////////////////////
// for debugging, set the triangle mesh to an icosahedron
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setToIcosahedron()
{
  _vertices.clear();
  _normals.clear();
  _texcoords.clear();
  _triangles.clear();

  Real phi = (1.0 + sqrt(5.0)) * 0.5;

  _vertices.push_back(VEC3F(0,  1,  phi));   // 11 -> 0 
  _vertices.push_back(VEC3F(0,  1,  -phi));  // 10 -> 1 
  _vertices.push_back(VEC3F(0,  -1,  phi));  // 0  -> 2 
  _vertices.push_back(VEC3F(0,  -1,  -phi)); // 9  -> 3 

  _vertices.push_back(VEC3F(1,  phi,  0));   // 6  -> 4 
  _vertices.push_back(VEC3F(1,  -phi,  0));  // 7  -> 5 
  _vertices.push_back(VEC3F(-1,  phi,  0));  // 5  -> 6 
  _vertices.push_back(VEC3F(-1,  -phi,  0)); // 8  -> 7 

  _vertices.push_back(VEC3F(phi,  0,  1));   // 1  -> 8 
  _vertices.push_back(VEC3F(phi,  0,  -1));  // 2  -> 9 
  _vertices.push_back(VEC3F(-phi,  0,  1));  // 4  -> 10 
  _vertices.push_back(VEC3F(-phi,  0,  -1)); // 3  -> 11 

  vector<VEC3I> faces;
//faces.push_back(VEC3I(1, 2, 6));
  faces.push_back(VEC3I(8, 9, 4));

//faces.push_back(VEC3I(1, 7, 2));
  faces.push_back(VEC3I(8, 5, 9));

//faces.push_back(VEC3I(3, 4, 5));
  faces.push_back(VEC3I(11, 10, 6));

//faces.push_back(VEC3I(4, 3, 8));
  faces.push_back(VEC3I(10, 11, 7));

//faces.push_back(VEC3I(6, 5, 11));
  faces.push_back(VEC3I(4, 6, 0));

//faces.push_back(VEC3I(5, 6, 10));
  faces.push_back(VEC3I(6, 4, 1));

//faces.push_back(VEC3I(9, 10, 2));
  faces.push_back(VEC3I(3, 1, 9));

//faces.push_back(VEC3I(10, 9, 3));
  faces.push_back(VEC3I(1, 3, 11));

//faces.push_back(VEC3I(7, 8, 9));
  faces.push_back(VEC3I(5, 7, 3));

//faces.push_back(VEC3I(8, 7, 0));
  faces.push_back(VEC3I(7, 5, 2));

//faces.push_back(VEC3I(11, 0, 1));
  faces.push_back(VEC3I(0, 2, 8));

//faces.push_back(VEC3I(0, 11, 4));
  faces.push_back(VEC3I(2, 0, 10));

//faces.push_back(VEC3I(6, 2, 10));
  faces.push_back(VEC3I(4, 9, 1));

//faces.push_back(VEC3I(1, 6, 11));
  faces.push_back(VEC3I(8, 4, 0));

//faces.push_back(VEC3I(3, 5, 10));
  faces.push_back(VEC3I(11, 6, 1));

//faces.push_back(VEC3I(5, 4, 11));
  faces.push_back(VEC3I(6, 10, 0));

//faces.push_back(VEC3I(2, 7, 9));
  faces.push_back(VEC3I(9, 5, 3));

//faces.push_back(VEC3I(7, 1, 0));
  faces.push_back(VEC3I(5, 8, 2));

//faces.push_back(VEC3I(3, 9, 8));
  faces.push_back(VEC3I(11, 3, 7));

//faces.push_back(VEC3I(4, 8, 0));
  faces.push_back(VEC3I(10, 7, 2));

  // store as triangles
  for (unsigned int x = 0; x < faces.size(); x++)
    _triangles.push_back(TRIANGLE(&_vertices[faces[x][0]],
                                  &_vertices[faces[x][1]],
                                  &_vertices[faces[x][2]]));

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeOBJ(const string path, const int frame)
{
  // build the filename
  string filename = path;

  // check for the slash
  if (filename[filename.length() - 1] != '/')
    filename = filename + string("/");

  // add the frame number
  filename = filename + string("surface.");
  char buffer[256];
  sprintf(buffer, "%04i", frame);
  filename = filename + string(buffer) + string(".obj");

  return writeOBJ(filename);
}

////////////////////////////////////////////////////////////////////////////
// write the raw data to a file stream
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::write(FILE* file)
{
  // write out dims
  int totalVertices = _vertices.size();
  int totalTriangles = _triangles.size();
  fwrite((void*)&totalVertices, sizeof(int), 1, file);
  fwrite((void*)&totalTriangles, sizeof(int), 1, file);

  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x].write(file);
  
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    int indices[3];
    indices[0] = _vertexIndices[_triangles[x].vertex(0)];
    indices[1] = _vertexIndices[_triangles[x].vertex(1)];
    indices[2] = _vertexIndices[_triangles[x].vertex(2)];

    fwrite((void*)&indices[0], sizeof(int), 3, file);
  }
}

////////////////////////////////////////////////////////////////////////////
// read in the raw data from a file stream
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::read(FILE* file)
{
  // read in dims
  int totalVertices;
  int totalTriangles;
  fread((void*)&totalVertices, sizeof(int), 1, file);
  fread((void*)&totalTriangles, sizeof(int), 1, file);

  _vertices.resize(totalVertices);
  for (int x = 0; x < totalVertices; x++)
    _vertices[x].read(file);
 
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  _triangles.resize(totalTriangles); 
  for (int x = 0; x < totalTriangles; x++)
  {
    int indices[3];
    fread((void*)&indices[0], sizeof(int), 3, file);

    _triangles[x].vertex(0) = &_vertices[indices[0]];
    _triangles[x].vertex(1) = &_vertices[indices[1]];
    _triangles[x].vertex(2) = &_vertices[indices[2]];
  }
}

////////////////////////////////////////////////////////////////////////////
// compute some coarse normals for each vertex
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNormals()
{
  _normals.clear();
  _normals.resize(_vertices.size());

  // instead of an area weighting, we'll store each normal un-normalized,
  // since the cross product magnitude will reflect the area of the triangle
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    TRIANGLE& triangle = _triangles[x];
    vector<VEC3F*>& vertices = triangle.vertices();
    VEC3F normal = cross(*vertices[1] - *vertices[0], 
                         *vertices[2] - *vertices[0]);

    int indices[] = {_vertexIndices[vertices[0]], _vertexIndices[vertices[1]], 
                     _vertexIndices[vertices[2]]};
    _normals[indices[0]] += normal;
    _normals[indices[1]] += normal;
    _normals[indices[2]] += normal;
  }
  
  // go through and normalize after everything
  for (unsigned int x = 0; x < _normals.size(); x++)
    _normals[x].normalize(); 
}

////////////////////////////////////////////////////////////////////////////
// write out an oriented point cloud to PLY
////////////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writePlyPointCloud(const string& filename)
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

  int totalVertices = _vertices.size();
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
    //float vertex[] = {_vertices[x][0], _vertices[x][1], _vertices[x][2]};
    float vertex[] = {(float)_vertices[x][0], 
                      (float)_vertices[x][1], 
                      (float)_vertices[x][2]};
    fwrite((void*)&(vertex[0]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[1]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[2]), sizeof(float), 1, file);
  }

  // close the binary file
  fclose(file);

  return true;
}

////////////////////////////////////////////////////////////////////////////
// write out an oriented point cloud to PLY
////////////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writePlyOrientedPointCloud(const string& filename)
{
  FILE* file = NULL;
  file = fopen(filename.c_str(), "w");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Couldn't open file " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  // if there are no normals, go ahead and compute them
  if (_vertices.size() != _normals.size())
    computeNormals();

  fprintf(file, "ply\n");
  //fprintf(file, "format ascii 1.0\n");
  fprintf(file, "format binary_little_endian 1.0\n");

  int totalVertices = _vertices.size();
  fprintf(file, "element vertex %i\n", totalVertices);
  fprintf(file, "property float x\n");
  fprintf(file, "property float y\n");
  fprintf(file, "property float z\n");
  fprintf(file, "property float nx\n");
  fprintf(file, "property float ny\n");
  fprintf(file, "property float nz\n");
  fprintf(file, "element face 0\n");
  fprintf(file, "property list uchar int vertex_index\n");
  fprintf(file, "end_header\n");

  /*
  for (int x = 0; x < totalVertices; x++)
    fprintf(file, "%f %f %f %f %f %f\n", 
        _vertices[x][0], _vertices[x][1], _vertices[x][2],
        _normals[x][0], _normals[x][1], _normals[x][2]);
        */

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
    //float vertex[] = {_vertices[x][0], _vertices[x][1], _vertices[x][2]};
    float vertex[] = {(float)_vertices[x][0], 
                      (float)_vertices[x][1], 
                      (float)_vertices[x][2]};
    //float normal[] = {_normals[x][0], _normals[x][1], _normals[x][2]};
    float normal[] = {(float)_normals[x][0],
                      (float)_normals[x][1], 
                      (float)_normals[x][2]};
    fwrite((void*)&(vertex[0]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[1]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[2]), sizeof(float), 1, file);
    fwrite((void*)&(normal[0]), sizeof(float), 1, file);
    fwrite((void*)&(normal[1]), sizeof(float), 1, file);
    fwrite((void*)&(normal[2]), sizeof(float), 1, file);
  }

  // close the binary file
  fclose(file);
  return true;
}

////////////////////////////////////////////////////////////////////////////
// get the signed distance field
////////////////////////////////////////////////////////////////////////////
FIELD_3D TRIANGLE_MESH::computeSignedDistanceField()
{
  const int res = 57;
  VEC3F min, max;
  boundingBox(min, max);
  VEC3F origin = (min + max) / (Real)2.0f;

  // get the dimensions
  VEC3F diff = max - min;
  Real maxLength = diff[0];
  if (diff[1] > maxLength) maxLength = diff[1];
  if (diff[2] > maxLength) maxLength = diff[2];

  maxLength *= 1.0 + 1.0 / 16.0;
  VEC3F lengths(maxLength, maxLength, maxLength);
  FIELD_3D final(res, res, res, origin, lengths);

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " UNSUPPORTED " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //initializeSignedDistanceField(final);
  //final.fastMarchingMethod();

  return final;
}

////////////////////////////////////////////////////////////////////////////
// write out geometry to Gmsh format
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeGmsh(const string& filename)
{
  cout << " Writing Gmsh file " << filename.c_str() << " ... " << flush;
  FILE* file = NULL;
  file = fopen(filename.c_str(), "w");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Couldn't open file " << filename.c_str() << "!!!" << endl;
    exit(0);
  }
  fprintf(file, "$MeshFormat\n");
  fprintf(file, "2.2 0 8\n");
  fprintf(file, "$EndMeshFormat\n");
  fprintf(file, "$Nodes\n");

  fprintf(file, "%i\n", (int)_vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    fprintf(file, "%i %f %f %f\n", x + 1, (double)_vertices[x][0], (double)_vertices[x][1], (double)_vertices[x][2]);
  fprintf(file, "$EndNodes\n");
  fprintf(file, "$Elements\n");
  
  fprintf(file, "%i\n", (int)_triangles.size());

  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    int indices[] = {_vertexIndices[_triangles[x].vertex(0)], _vertexIndices[_triangles[x].vertex(1)], _vertexIndices[_triangles[x].vertex(2)]};
    fprintf(file, "%i 2 2 1 1 %i %i %i\n", x + 1, indices[0] + 1, indices[1] + 1, indices[2] + 1);
  }
  fprintf(file, "$EndElements\n");
  fclose(file);
  cout << " done. " << endl;
}

////////////////////////////////////////////////////////////////////////////
// emulate the cell center lookup for FIELD_3D --
// assumes that _xRes, _yRes, _zRes,
//              _center, _lengths, and _dx are populated
////////////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::cellCenter(int x, int y, int z)
{
  VEC3F halfLengths = (Real)0.5 * _lengths;

  // set it to the lower corner
  VEC3F final = _center - halfLengths;

  // displace to the NNN corner
  final[0] += x * _dxs[0];
  final[1] += y * _dxs[1];
  final[2] += z * _dxs[2];

  // displace it to the cell center
  final[0] += _dxs[0] * 0.5;
  final[1] += _dxs[1] * 0.5;
  final[2] += _dxs[2] * 0.5;

  return final;
}

////////////////////////////////////////////////////////////////////////////
// compute all the slices for a low memory marching cubes
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeAllLowMemorySlices(vector<pair<int, int> >& flags)
{
  TIMER functionTimer(__FUNCTION__);

  // number of nans and infs
  int totalNans = 0;
  int totalInfs = 0;

  FIELD_2D& slab0 = _slab0;
  FIELD_2D& slab1 = _slab1;
  computeNonlinearSlice(0, slab1);
  totalNans += slab1.totalNans();
  totalInfs += slab1.totalInfs();

  // build all the vertex pairs 
  _vertexPairs.clear();
  for (int z = 0; z < _zRes - 1; z++)
  {
    // swap in the old "next" slice as the new current ont
    FIELD_2D& old = slab0;
    slab0 = slab1;
    slab1 = old;

    // compute the next needed slice on the fly
    computeNonlinearSlice(z + 1, slab1);
    totalNans += slab1.totalNans();
    totalInfs += slab1.totalInfs();

    for (int y = 0; y < _yRes - 1; y++)
      for (int x = 0; x < _xRes - 1; x++) 
      {
        int index = x + y * _xRes + z * _slabSize;

        CUBE cube;
        cube.NNN = slab0(x,y);
        cube.NNP = slab1(x,y);
        cube.NPN = slab0(x,y + 1);
        cube.NPP = slab1(x,y + 1);
        cube.PNN = slab0(x + 1,y);
        cube.PNP = slab1(x + 1,y);
        cube.PPN = slab0(x + 1,y + 1);
        cube.PPP = slab1(x + 1,y + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));

        if (flag == 0 || flag == 255) continue;

        flags.push_back(pair<int, int>(flag, index));
		
        // three vertices are added to _vertexPairs here  
        switch (flag)
#include "MARCHING_CUBES_VERTICES.include" 
      }
    if (z % (int)(_zRes / 10) == 0)
      cout << 100 * ((Real)z / _zRes) << "% " << flush;
  }
  cout << " done." << endl;
  cout << " infs: " << totalInfs << endl;
  cout << " NaNs: " << totalNans << endl;
}

////////////////////////////////////////////////////////////////////////////
// compute all the slices for a low memory marching cubes
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeAllLowMemorySlicesHuge(vector<pair<int, VEC3I> >& flags)
{
  TIMER functionTimer(__FUNCTION__);

  // number of nans and infs
  int totalNans = 0;
  int totalInfs = 0;

  FIELD_2D& slab0 = _slab0;
  FIELD_2D& slab1 = _slab1;
  computeNonlinearSlice(0, slab1);
  totalNans += slab1.totalNans();
  totalInfs += slab1.totalInfs();

  // build all the vertex pairs 
  _vertexTriplets.clear();
  for (int z = 0; z < _zRes - 1; z++)
  {
    // swap in the old "next" slice as the new current ont
    FIELD_2D& old = slab0;
    slab0 = slab1;
    slab1 = old;

    // compute the next needed slice on the fly
    computeNonlinearSlice(z + 1, slab1);
    totalNans += slab1.totalNans();
    totalInfs += slab1.totalInfs();

    for (int y = 0; y < _yRes - 1; y++)
      for (int x = 0; x < _xRes - 1; x++) 
      {
        //int index = x + y * _xRes + z * _slabSize;
        VEC3I index(x,y,z);

        CUBE cube;
        cube.NNN = slab0(x,y);
        cube.NNP = slab1(x,y);
        cube.NPN = slab0(x,y + 1);
        cube.NPP = slab1(x,y + 1);
        cube.PNN = slab0(x + 1,y);
        cube.PNP = slab1(x + 1,y);
        cube.PPN = slab0(x + 1,y + 1);
        cube.PPP = slab1(x + 1,y + 1);
		
        // construct the flag
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));

        if (flag == 0 || flag == 255) continue;

        flags.push_back(pair<int, VEC3I>(flag, index));
		
        // three vertices are added to _vertexPairs here  
        switch (flag)
//#include "MARCHING_CUBES_VERTICES.include" 
#include "MARCHING_CUBES_VERTICES.include.huge" 
      }
    if (z % (int)(_zRes / 10) == 0)
      cout << 100 * ((Real)z / _zRes) << "% " << flush;
  }
  cout << " done." << endl;
  cout << " infs: " << totalInfs << endl;
  cout << " NaNs: " << totalNans << endl;
}
