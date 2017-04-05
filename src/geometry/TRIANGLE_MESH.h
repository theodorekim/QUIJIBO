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
//////////////////////////////////////////////////////////////////////
// Triangle Mesh class
//////////////////////////////////////////////////////////////////////

#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include <map>
#include <SETTINGS.h>
#include <TRIANGLE.h>
#include <FIELD_3D.h>
#include <POLYNOMIAL_4D.h>

using namespace std;
class TRIANGLE_MESH
{
public:
  TRIANGLE_MESH();

  // OBJ mesh constructor
  TRIANGLE_MESH(const string& filename);

  // marching cubes constructor
  TRIANGLE_MESH(const FIELD_3D& field);
  
  // non-linear marching cubes constructor
  TRIANGLE_MESH(const FIELD_3D& field, const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const Real expScaling, const int maxIterations, const Real slice, const Real isosurface);

  // do a non-linear marching cubes on just two slabs at a shot
  TRIANGLE_MESH(const VEC3F& center, 
                const VEC3F& lengths, 
                const VEC3I& res, 
                const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const Real expScaling, const int maxIterations, const Real slice, const Real isosurface);
  
  // non-linear marching cubes for just the quadratic Julia set
  TRIANGLE_MESH(const FIELD_3D& field, const QUATERNION& qConst, const int maxIterations);

  // destructor
  ~TRIANGLE_MESH();

  // Normalize mesh to 1x1x1 cube centered at (0,0,0)
  void normalize();
  void normalize(Real padding);

  //get the range of the original mesh
  void getBounds(VEC3F& maxVert, VEC3F& minVert, Real& maxLength);

  const vector<TRIANGLE>& triangles() const { return _triangles; };
  const vector<VEC3F>& vertices() const { return _vertices; };
  vector<VEC3F>& vertices() { return _vertices; };
  vector<VEC3F>& texcoords() { return _texcoords; };
  int vertexIndex(int triangle, int vertex) { return _vertexIndices[_triangles[triangle].vertex(vertex)]; };
  FIELD_2D& texture() { return _texture; };

  // overloaded SURFACE functions
  virtual bool inside(const VEC3F& point);
  virtual Real distance(const VEC3F& point);
  virtual Real signedDistance(const VEC3F& point);

  virtual float distance(float* point){return 0.f;}

  // initialize a signed distance field by only setting values along
  // the surface -- replace the slow default version with
  // one that only looks at grid cells in a triangle's neighborhood
  virtual void initializeSignedDistanceField(FIELD_3D& field);

  // initialize an extension scalar to a field
  virtual void initializeExtensionScalars(const FIELD_3D& distance, FIELD_3D& scalar);

  // assign a texture coordinates with spherical projection
  void setSphericalTexcoords();
  
  // assign a texture coordinates with cylindrical projection
  void setCylindricalTexcoords(float bottomCap, float topCap);
  
  // assign a texture coordinates with swirled spherical projection,
  // just to get something more non-trivial for debugging purposes
  void setSwirledTexcoords();

  // center using vertex means
  VEC3F vertexMean() const;

  // return a bounding box
  void boundingBox(VEC3F& mins, VEC3F& maxs);

  // create a set of textures based on a solid texture
  void textureUsingSolidTexture(const FIELD_3D& solidTexture, const int textureRes = 5);
  void normalizeSolidTexture();

  // get the center coordinate of the triangle with the maximum texture value
  VEC3F maxTextureCenter();

  // perform marching cubes
  void computeMarchingCubes(const FIELD_3D& field, const bool verbose = false);

  // compute marching cubes in stages, in preparation for the non-linear solve
  void computeStagedMarchingCubes(const FIELD_3D& field, const bool verbose = false);
  
  // compute marching cubes in stages, with non-linear solve in between
  void computeNonlinearMarchingCubes(const FIELD_3D& field, const bool verbose = false);
  
  // compute marching cubes in stages, with non-linear solve in between
  void computeQuadraticMarchingCubes(const FIELD_3D& field, const bool verbose = false);

  // get the signed distance field
  FIELD_3D computeSignedDistanceField();

  ////////////////////////////////////////////////////////////////////////////
  // read/write support
  ////////////////////////////////////////////////////////////////////////////
  int sizeTriangles(){return _triangles.size();}
  int sizeVertices(){return _vertices.size();}	
  bool readOBJ(const string& filename);
  bool writeOBJ(const string& filename);
  bool writeOBJ(const string path, const int frame);
  bool writeOBJgz(const string& filename);

  // write out a point cloud to PLY
  bool writePlyPointCloud(const string& filename);
  
  // write out an oriented point cloud to PLY
  bool writePlyOrientedPointCloud(const string& filename);

  // write out LuxRender geometry
  void writeLuxRender(const char* filename, const char* path);

  // write out PBRT geometry
  void writePBRT(const char* filename);

  // write out geometry as bobj gz
  void writeBobjGz(const char* filename);

  // write out a BREP file for Sean Mauch's closest point transform (CPT) code
  void writeBRep(const char* filename);

  // write out geometry to Gmsh format
  void writeGmsh(const string& filename);

  // write the raw data to a file stream
  void write(FILE* file);
  
  // read from a raw file stream
  void read(FILE* file);

  ////////////////////////////////////////////////////////////////////////////
  // test cases
  ////////////////////////////////////////////////////////////////////////////
  void setToTet();

  void setToIcosahedron();

private:
  vector<VEC3F> _vertices;
  vector<VEC3F> _normals;
	vector<VEC3F> _texcoords;
  vector<TRIANGLE> _triangles;

  // hash table allowing lookup of vertex index from address
  map<VEC3F*, int> _vertexIndices;

  // distance field dims from marching cubes
  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;

  // hash table that allows quick lookup of if a vertex has been created
  // by marching cubes previously
  map<int, vector<int> > _vertexHash;

  ////////////////////////////////////////////////////////////////////////////
  // texturing support
  ////////////////////////////////////////////////////////////////////////////
 
  // the texture, if any
  FIELD_2D _texture;

  // GL handle for the texture
  unsigned int _glTextureHandle;

  ////////////////////////////////////////////////////////////////////////////
  // solid texturing support
  ////////////////////////////////////////////////////////////////////////////

  // per triangle texture
  vector<FIELD_2D> _solidTextures;
  
  // per triangle mean solid texture values
  vector<float> _solidTextureMeans;
  
  // GL handle for the texture
  //vector<GLuint> _glSolidTextureHandles;

  ////////////////////////////////////////////////////////////////////////////
  // Geometric tests
  ////////////////////////////////////////////////////////////////////////////

  // ray-triangle test from
  // Tomas Möller and Ben Trumbore. Fast, minimum storage ray-triangle intersection. 
  // Journal of graphics tools, 2(1):21-28, 1997
  bool intersectTriangle(const VEC3F& orig, const VEC3F& dir, int faceIndex);

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
  Real pointFaceDistanceSq(TRIANGLE& f, const VEC3F& point);

  ////////////////////////////////////////////////////////////////////////////
  // Marching Cubes support
  ////////////////////////////////////////////////////////////////////////////

  // cached marching cube distances
  Real _NNN, _NNP, _NPN, _NPP, _PNN, _PNP, _PPN, _PPP; 

  // cached marching cube deltas
  VEC3F _fieldDeltas;

  // cell center of the current grid cell being marched
  VEC3F _cellCenter;

  // what is considered "outside" by grid being marched
  Real _outside;

  // store index triplets initially for the triangles, and then set everything
  // to pointers once the vector is done being resized
  vector<int> _triangleVertices;

  // struct needed to pass data during staged marching cubes
  struct CUBE {
    union {
      struct { float NNN, NNP, NPN, NPP, PNN, PNP, PPN, PPP; };
      float data[8];
    };
  };

  // vertex pairs to compute an interpolation point between
  map<pair<int, int>, bool> _vertexPairs;
  
  // where in _vertices is the vertex that corresponds to this pair?
  map<pair<int, int>, int> _vertexPairHash;

  // the field being marching cubed
  const FIELD_3D* _toMarch;

  // polynomials for non-linear marching cubes
  POLYNOMIAL_4D _top;
  POLYNOMIAL_4D _bottom;
  Real _expScaling;
  int _maxIterations;
  Real _quaternionSlice;
  Real _isosurface;

  Real _escapeRadius;
 
  // constant for just quadratic Julia set
  QUATERNION _quadraticConst;

  // two-slab non-linear marching cubes vars
  const VEC3I _res;
  const VEC3F _lengths;
  const VEC3F _center;
  const VEC3F _dxs;
  FIELD_2D _slab0;
  FIELD_2D _slab1;

  // add a triangle to the list
  void addTriangle(int i, int j, int k, int index);
  void addStagedTriangle(int i, int j, int k, int index, const CUBE& cube, const VEC3F& center);

  // add a triangle to the list whose vertices were all precomputed
  // by computeEdgeInterpolations()
  void addVertexPairTriangle(int i, int j, int k, int index);

  // add a vertex pairs, to be interpolated later
  void addVertexPairs(int i, int j, int k, int index);
  pair<int, int> getVertexPair(int i, int index);
  
  // get the edge point
  VEC3F computeVertex(int i, int index);
  VEC3F computeStagedVertex(int i, int index, const CUBE& cube, const VEC3F& center);

  // see if the vertex has been computed before, and if not, store it
  int storeVertex(VEC3F& vertex, int index);

  // compute the linear interpolations for matching cubes
  void computeEdgeInterpolations();
  
  // compute the non-linear interpolations for matching cubes
  void computeNonlinearEdgeInterpolations();

  // compute the quadratic interpolations for matching cubes
  void computeQuadraticEdgeInterpolations();

  // compute the non-linear interpolations for matching cubes, but in a memory-stingy way
  void computeNonlinearMarchingCubesLowMemory();

  // support function for computeNonlinearMarchingCubesLowMemory(), which computes a single
  // slice of the potential function
  void computeNonlinearSlice(const int z, FIELD_2D& field);

  // compute some coarse normals for each vertex
  void computeNormals();

  // get the (x,y,z) of an index
  VEC3I getXYZ(const int index) const;
  
  // get the nonlinear function value
  Real nonlinearValue(const VEC3F& position, const bool debug = false);
  
  // get the quadratic function value
  Real quadraticValue(const VEC3F& position, const bool debug = false);

  // do a midpoint search
  VEC3F midpointSearch(const VEC3F& positiveVertex, const Real& positiveValue, const VEC3F& negativeVertex, const Real& negativeValue, const int recursion = 0);
  
  // do a midpoint search over just the quadratic
  VEC3F quadraticMidpointSearch(const VEC3F& positiveVertex, const Real& positiveValue, const VEC3F& negativeVertex, const Real& negativeValue, const int recursion = 0);

  // compute a subset of the non-linear interpolations
  //void computeInterpolationSubset(int begin, int end);

  // emulate the cell center lookup for FIELD_3D --
  // assumes that _xRes, _yRes, _zRes,
  //              _center, _lengths, and _dx are populated
  VEC3F cellCenter(int x, int y, int z);
};

#endif
