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

// This whole project is just a wrapper to call Chris Batty's SDFGen:
//
//   https://github.com/christopherbatty/SDFGen

//SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
//Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
//...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
//This code is public domain. Feel free to mess with it, let me know if you like it.

#include "makelevelset3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include "FIELD_3D.h"
#include "SIMPLE_PARSER.h"

void reweightingRun(SIMPLE_PARSER& parser)
{
  cout << " DOING A REWEIGHTING RUN " << endl;

  std::string path = parser.getString("path", "./temp/");

  // get the FIELD_3D of the original OBJ
  std::string originalDistanceFilename = path + parser.getString("distance field", "temp.field3d");
  FIELD_3D originalDistance(originalDistanceFilename.c_str());
 
  // get dimensions out of the original file
  Vec3f lengths;
  Vec3f center;
  lengths[0] = originalDistance.lengths()[0];
  lengths[1] = originalDistance.lengths()[1];
  lengths[2] = originalDistance.lengths()[2];
  center[0] = originalDistance.center()[0];
  center[1] = originalDistance.center()[1];
  center[2] = originalDistance.center()[2];

  Vec3f min_box = center - lengths * 0.5;
  Vec3f max_box = center + lengths * 0.5;

  std::string filename("temp.obj");
  if(filename.size() < 5 || filename.substr(filename.size()-4) != std::string(".obj")) {
    std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
    exit(-1);
  }

  int res = -1;
  res = parser.getInt("res", res);
  float dx = -1;
 
  std::cout << "Reading data.\n";

  std::ifstream infile(filename.c_str());
  if(!infile) {
    std::cerr << "Failed to open. Terminating.\n";
    exit(-1);
  }

  int ignored_lines = 0;
  std::string line;
  std::vector<Vec3f> vertList;
  std::vector<Vec3ui> faceList;
  while(!infile.eof()) {
    std::getline(infile, line);
    if(line.substr(0,1) == std::string("v")) {
      std::stringstream data(line);
      char c;
      Vec3f point;
      data >> c >> point[0] >> point[1] >> point[2];
      vertList.push_back(point);
    }
    else if(line.substr(0,1) == std::string("f")) {
      std::stringstream data(line);
      char c;
      int v0,v1,v2;
      data >> c >> v0 >> v1 >> v2;
      faceList.push_back(Vec3ui(v0-1,v1-1,v2-1));
    }
    else {
      ++ignored_lines; 
    }
  }
  infile.close();

  // scale according to the magic number later used in main() of rootOptimizerTAO3D
  for (unsigned int x = 0; x < vertList.size(); x++)
  {
    vertList[x] *= 0.25;
    vertList[x] += Vec3f(0.5, 0.5, 0.5);
  }
  
  if(ignored_lines > 0)
    std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";

  std::cout << "Read in " << vertList.size() << " vertices and " << faceList.size() << " faces." << std::endl;

  double maxLength = lengths[0] > lengths[1] ? lengths[0] : lengths[1];
  maxLength = maxLength > lengths[2] ? maxLength : lengths[2];
  cout << " center: " << center << endl;
  cout << " lengths: " << lengths << endl;
  cout << " max length found " << maxLength << endl;
  dx = lengths[0] / res;
  cout << " dx: " << dx << endl;

  assert(dx > 0);
  Vec3ui sizes = Vec3ui((max_box - min_box)/dx);
  
  std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;
  std::cout << "Computing signed distance field.\n";
  Array3f phi_grid;
  make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);
  std::cout << "Processing complete.\n";

  // create a FIELD_3D object to save
  int xRes = phi_grid.ni;
  int yRes = phi_grid.nj;
  int zRes = phi_grid.nk;
  VEC3F myLengths(lengths[0], lengths[1], lengths[2]);
  VEC3F myCenter(center[0], center[1], center[2]);
  cout << " Writing field3d file with center: " << myCenter << " lengths: " << myLengths << endl;

  FIELD_3D field(xRes, yRes, zRes, myCenter, myLengths);
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
        field(x,y,z) = phi_grid(x,y,z);

  field.write("temp.distance.field3d");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The original version, assuming it is not a reweighting run
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) 
{
  if(argc < 2) {
    std::cout << "SDFGen - A utility for converting closed oriented triangle meshes into grid-based signed distance fields.\n";
    std::cout << "\nThe output file format is:";
    std::cout << "<ni> <nj> <nk>\n";
    std::cout << "<origin_x> <origin_y> <origin_z>\n";
    std::cout << "<dx>\n";
    std::cout << "<value_1> <value_2> <value_3> [...]\n\n";
    
    std::cout << "(ni,nj,nk) are the integer dimensions of the resulting distance field.\n";
    std::cout << "(origin_x,origin_y,origin_z) is the 3D position of the grid origin.\n";
    std::cout << "<res> is the number of grid points in a single direction.\n\n";
    std::cout << "<value_n> are the signed distance data values, in ascending order of i, then j, then k.\n";

    std::cout << "The output filename will match that of the input, with the OBJ suffix replaced with SDF.\n\n";

    std::cout << "Usage: SDFGen <config file> <is this a reweighting run? (optional)> \n\n";
    exit(-1);
  }

  SIMPLE_PARSER parser(argv[1]);

  // if we're doing a reweight run, fork now please.
  if (argc > 2)
  {
    reweightingRun(parser);
    return 0;
  }

  std::string path = parser.getString("path", "./temp/");

  std::string filename = path + parser.getString("normalized obj", "dummy.obj");
  if(filename.size() < 5 || filename.substr(filename.size()-4) != std::string(".obj")) {
    std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
    exit(-1);
  }

  int res = -1;
  res = parser.getInt("res", res);
  float dx = -1;
 
  //start with a massive inside out bound box.
  Vec3f min_box(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max()), 
    max_box(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
  
  std::cout << "Reading data.\n";

  std::ifstream infile(filename.c_str());
  if(!infile) {
    std::cerr << "Failed to open. Terminating.\n";
    exit(-1);
  }

  int ignored_lines = 0;
  std::string line;
  std::vector<Vec3f> vertList;
  std::vector<Vec3ui> faceList;
  while(!infile.eof()) {
    std::getline(infile, line);
    if(line.substr(0,1) == std::string("v")) {
      std::stringstream data(line);
      char c;
      Vec3f point;
      data >> c >> point[0] >> point[1] >> point[2];
      vertList.push_back(point);
      update_minmax(point, min_box, max_box);
    }
    else if(line.substr(0,1) == std::string("f")) {
      std::stringstream data(line);
      char c;
      int v0,v1,v2;
      data >> c >> v0 >> v1 >> v2;
      faceList.push_back(Vec3ui(v0-1,v1-1,v2-1));
    }
    else {
      ++ignored_lines; 
    }
  }
  infile.close();
  
  if(ignored_lines > 0)
    std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";

  std::cout << "Read in " << vertList.size() << " vertices and " << faceList.size() << " faces." << std::endl;

  // find the center
  Vec3f center = (min_box + max_box) * 0.5;
  Vec3f lengths = max_box - min_box;
  double maxLength = lengths[0] > lengths[1] ? lengths[0] : lengths[1];
  maxLength = maxLength > lengths[2] ? maxLength : lengths[2];
  cout << " center: " << center << endl;
  cout << " lengths: " << lengths << endl;
  cout << " max length found " << maxLength << endl;
  lengths = Vec3f(maxLength, maxLength, maxLength);

  min_box = center - lengths * 0.5;
  max_box = center + lengths * 0.5;

  dx = maxLength / res;
  cout << " dx: " << dx << endl;

  //Add padding around the box.
  int padding = 1;
  Vec3f unit(1,1,1);
  min_box -= padding*dx*unit;
  max_box += padding*dx*unit;
  cout << " min: " << min_box[0] << " " << min_box[1] << " " << min_box[2] << endl;
  cout << " max: " << max_box[0] << " " << max_box[1] << " " << max_box[2] << endl;

  // now re-do dx to include the padding as well
  lengths = max_box - min_box;
  dx = lengths[0] / res;
  cout << " lengths: " << lengths[0] << " " << lengths[1] << " " << lengths[2] << endl;
  cout << " dx: " << dx << endl;

  assert(dx > 0);

  Vec3ui sizes = Vec3ui((max_box - min_box)/dx);
  
  std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

  std::cout << "Computing signed distance field.\n";
  Array3f phi_grid;
  make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);

  std::cout << "Processing complete.\n";

  // create a FIELD_3D object to save
  int xRes = phi_grid.ni;
  int yRes = phi_grid.nj;
  int zRes = phi_grid.nk;
  VEC3F myLengths(lengths[0], lengths[1], lengths[2]);
  VEC3F myCenter(center[0], center[1], center[2]);
  cout << " Writing field3d file with center: " << myCenter << " lengths: " << myLengths << endl;

  FIELD_3D field(xRes, yRes, zRes, myCenter, myLengths);
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
        field(x,y,z) = phi_grid(x,y,z);

  string outname = path + parser.getString("distance field", "dummy.field3d");
  field.write(outname.c_str());

  return 0;
}
