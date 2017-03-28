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

// This whole project is just a wrapper to call the blue noise
// functions in vcglib:
//
//   http://vcg.isti.cnr.it/vcglib/

#include <vcg/complex/complex.h>

// include the algorithms for updating:
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>

// include the algorithms for mesh fixing
#include <vcg/complex/algorithms/clean.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

// include the support for polygon meshes (function to convert from/to trimesh)
#include <vcg/complex/algorithms/polygon_support.h>

// include the support for polygon meshes (the component for the face )
#include <vcg/simplex/face/component_polygon.h>

// include the support for half edges
#include <vcg/complex/algorithms/update/halfedge_indexed.h>
#include<vcg/complex/algorithms/point_sampling.h>

#include "FIELD_3D.h"
#include "VECTOR3_FIELD_3D.h"
#include "QUATERNION.h"
#include "OPTIMIZE_3D.h"
#include "SIMPLE_PARSER.h"

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::VFAdj, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::VFAdj, face::FFAdj, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

//////////////////////////////////////////////////////////////////////////////////
// Dump out a PLY point cloud
//////////////////////////////////////////////////////////////////////////////////
bool writePlyPointCloud(vector<Point3f>& points, const string& filename)
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
  fprintf(file, "format binary_little_endian 1.0\n");

  int totalPoints = points.size();
  fprintf(file, "element vertex %i\n", totalPoints);
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

  Real maxMagnitude = 0;
  for (int x = 0; x < totalPoints; x++)
  {
    float vertex[] = {points[x][0], points[x][1], points[x][2]};
    fwrite((void*)&(vertex[0]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[1]), sizeof(float), 1, file);
    fwrite((void*)&(vertex[2]), sizeof(float), 1, file);

    Real magnitude = sqrt(points[x][0] * points[x][0] + points[x][1] * points[x][1] + points[x][2] * points[x][2]);

    if (magnitude > maxMagnitude)
    {
      cout << " found max: " << magnitude << endl;
      maxMagnitude = magnitude;
    }
  }

  // close the binary file
  fclose(file);

  return true;
}

template <class MeshType>
void WeightedPoissonSampling(MeshType &m, // the mesh that has to be sampled
                     std::vector<Point3f> &poissonSamples, // the vector that will contain the set of points
                     int sampleNum, // the desired number sample, if zero you must set the radius to the wanted value
                     float &radius,  // the Poisson Disk Radius (used if sampleNum==0, setted if sampleNum!=0)
                     const FIELD_3D& weights,
                     float radiusVariance=1,
                     int randSeed=123456)
{
  typedef tri::TrivialSampler<MeshType> BaseSampler;
  typedef tri::MeshSampler<MeshType> MontecarloSampler;
  tri::RequireVFAdjacency(m);

  typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam pp;
  typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam::Stat stat;
  pp.pds = &stat;
  int t0=clock();

  if(sampleNum>0) radius = tri::SurfaceSampling<MeshType,BaseSampler>::ComputePoissonDiskRadius(m,sampleNum);
  if(radius>0 && sampleNum==0) sampleNum = tri::SurfaceSampling<MeshType,BaseSampler>::ComputePoissonSampleNum(m,radius);

  pp.pds->sampleNum = sampleNum;
  poissonSamples.clear();
  MeshType MontecarloMesh;

  // First step build the sampling
  MontecarloSampler mcSampler(MontecarloMesh);
  BaseSampler pdSampler(poissonSamples);

  tri::SurfaceSampling<MeshType,MontecarloSampler>::Montecarlo(m, mcSampler, std::max(10000,sampleNum*40));
  tri::UpdateBounding<MeshType>::Box(MontecarloMesh);
  int t1=clock();
  pp.pds->montecarloTime = t1-t0;

  if(radiusVariance !=1)
  {
    pp.adaptiveRadiusFlag=true;
    pp.radiusVariance=radiusVariance;
  }
  tri::SurfaceSampling<MeshType,BaseSampler>::WeightedPoissonDiskPruning(pdSampler, MontecarloMesh, radius,pp, weights, randSeed);
  int t2=clock();
  pp.pds->totalTime = t2-t0;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) 
{
  TIMER functionTimer(__FUNCTION__);
  if (argc < 2)
  {
    cout << " USAGE: " << argv[0] << " <cfg file> " << endl;
    exit(0);
  }

  SIMPLE_PARSER parser(argv[1]);

	int loadmask;
  MyMesh pm;

  string path = parser.getString("path", "./temp");
  string objFilename = path + parser.getString("normalized obj", "dummy.obj");
  cout << " Generating blue noise on surface of " << objFilename.c_str() << endl;
	
  // load up the OBJ
	vcg::tri::io::ImporterOBJ<MyMesh>::Open(pm,objFilename.c_str(),loadmask);

  // call Poisson sampling on the numerator
  vector<Point3f> points;
  float radius = 0;
  int totalPoints = -1;
  totalPoints = parser.getInt("poisson samples", totalPoints);

  cout << " Generating " << totalPoints << " points ... " << flush;
  PoissonSampling<MyMesh>(pm,points,totalPoints,radius);
  cout << " done. " << endl;
  cout << " total points: " << points.size() << endl;

  string plyFilename = path + parser.getString("poisson output", "temp.ply");

  cout << " Writing points to " << plyFilename.c_str() << endl;
  writePlyPointCloud(points, plyFilename.c_str());

  TIMER::printTimings();
}
