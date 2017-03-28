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

// This whole project is just a wrapper to call BEM++
//   http://www.bempp.org/

#include <cmath>
#include "VEC3.h"
#include "TRIANGLE_MESH.h"
#include "SIMPLE_PARSER.h"
#include "TIMER.h"
#include <cstdio>

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  TIMER functionTimer(__FUNCTION__);
  if (argc != 2)
  {
    cout << endl;
    cout << " USAGE: " << argv[0] << " <cfg file> " << endl;
    cout << endl;
    return 0;
  }

  SIMPLE_PARSER parser(argv[1]);

  string path = parser.getString("path", "./temp/");

  char pwd[1024];
  getcwd(pwd, 1024);
  printf("pwd -> %s\n", pwd);
  string fullPath = path;

  // strip off the leading "."
  fullPath.erase(0,1);
  fullPath = string(pwd) + fullPath;
  cout << " Full path: " << fullPath.c_str() << endl;

  string inputFilename = fullPath + parser.getString("msh filename", "dummy.msh");
  string distanceFilename = fullPath + parser.getString("distance field", "distance.field3D");
  string outputFilename = fullPath + parser.getString("potential field", "dummy.field3d");

  // TODO: If you're trying the advanced version, you will need to modify these
  // to point to your own build of tutorial_dirichlet in BEM++
  string cdCall("cd /Users/tedkim/bempp/build/bempp/examples");
  string exportCall("export DYLD_LIBRARY_PATH=/Users/tedkim/bempp/lib");

  // call the BEM solver
  char buffer[1024];
  sprintf(buffer,"./tutorial_dirichlet %s %s %s", inputFilename.c_str(), distanceFilename.c_str(), outputFilename.c_str());
  string bemCall(buffer);

  cout << " BEM call: " << bemCall.c_str() << endl;
  string concatCall = cdCall + string(";") + exportCall + string(";") + bemCall;
  cout << " Concatenated call: " << concatCall << endl;

  TIMER bemTimer("BEM Call");
  system(concatCall.c_str());
  bemTimer.stop();

  TIMER::printTimings();

  return 1;
}
