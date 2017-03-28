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
#include <cmath>

#include "VEC3.h"
#include "TRIANGLE_MESH.h"
#include "SIMPLE_PARSER.h"

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc != 2)
  {
    cout << endl;
    cout << " USAGE: " << argv[0] << " <cfg file> " << endl;
    cout << endl;
    return 0;
  }

  SIMPLE_PARSER parser(argv[1]);

  string path = parser.getString("path", "./temp/");

  string inputFilename = path + parser.getString("obj filename", "dummy.obj");
  string outputFilename = path + parser.getString("msh filename", "dummy.msh");
  string normalizedOutputFilename = path + parser.getString("normalized obj", "dummy.normalized.obj");

  Real padding = 1.0 / 16.0;
  padding = parser.getFloat("padding", padding);
  cout << " Using padding " << padding << endl;

  TRIANGLE_MESH inputMesh(inputFilename);
  inputMesh.normalize(padding);

  inputMesh.writeOBJ(normalizedOutputFilename);
  inputMesh.writeGmsh(outputFilename);
  return 1;
}
