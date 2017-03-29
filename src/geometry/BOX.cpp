/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
// BOX.cpp: implementation of the BOX class.
//
//////////////////////////////////////////////////////////////////////

#include "BOX.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BOX::BOX(float xPlus, float xMinus, 
         float yPlus, float yMinus, 
         float zPlus, float zMinus) :
  _xPlus(xPlus), _xMinus(xMinus), _yPlus(yPlus),
  _yMinus(yMinus), _zPlus(zPlus), _zMinus(zMinus)
{
}


BOX::BOX(const VEC3F& center, const VEC3F& sizes)
{
  _xPlus  = center[0] + sizes[0] * 0.5;
  _xMinus = center[0] - sizes[0] * 0.5;
  _yPlus  = center[1] + sizes[1] * 0.5;
  _yMinus = center[1] - sizes[1] * 0.5;
  _zPlus  = center[2] + sizes[2] * 0.5;
  _zMinus = center[2] - sizes[2] * 0.5;
}

BOX::~BOX()
{
}

static void drawNormal(VEC3F& n)
{
    glNormal3f(n[0], n[1], n[2]);
}

static void drawVertex(VEC3F& v)
{
    glVertex3f(v[0], v[1], v[2]);
}

void BOX::draw()
{
  VEC3F v000(_xMinus, _yMinus, _zMinus); 
  VEC3F v100(_xPlus, _yMinus, _zMinus); 
  VEC3F v010(_xMinus, _yPlus, _zMinus); 
  VEC3F v110(_xPlus, _yPlus, _zMinus); 
  VEC3F v001(_xMinus, _yMinus, _zPlus); 
  VEC3F v101(_xPlus, _yMinus, _zPlus); 
  VEC3F v011(_xMinus, _yPlus, _zPlus); 
  VEC3F v111(_xPlus, _yPlus, _zPlus); 

  glBegin(GL_QUADS);
    // x plus
    VEC3F normal = cross(v000 - v100, v000 - v110);
    normal.normalize();
    normal *= -1.0;
    drawNormal(normal);
    drawVertex(v010);
    drawVertex(v110);
    drawVertex(v100);
    drawVertex(v000);

    // x minus
    normal = cross(v001 - v101, v001 - v111);
    normal.normalize();
    drawNormal(normal);
    drawVertex(v001);
    drawVertex(v101);
    drawVertex(v111);
    drawVertex(v011);

    // y minus
    normal = cross(v000 - v100, v000 - v101);
    normal.normalize();
    drawNormal(normal);
    drawVertex(v000);
    drawVertex(v100);
    drawVertex(v101);
    drawVertex(v001);

    // y plus
    normal = cross(v010 - v110, v010 - v111);
    normal.normalize();
    normal *= -1.0;
    drawNormal(normal);
    drawVertex(v011);
    drawVertex(v111);
    drawVertex(v110);
    drawVertex(v010);

    // z plus
    normal = cross(v000 - v010, v000 - v011);
    normal.normalize();
    normal *= -1.0;
    drawNormal(normal);
    drawVertex(v001);
    drawVertex(v011);
    drawVertex(v010);
    drawVertex(v000);

    // z minus
    normal = cross(v100 - v110, v100 - v111);
    normal.normalize();
    drawNormal(normal);
    drawVertex(v100);
    drawVertex(v110);
    drawVertex(v111);
    drawVertex(v101);
  glEnd();
}

void BOX::translateX(float delta)
{
  _xPlus += delta;
  _xMinus += delta;
}

void BOX::translateY(float delta)
{
  _yPlus += delta;
  _yMinus += delta;
}

void BOX::translateZ(float delta)
{
  _zPlus += delta;
  _zMinus += delta;
}

void BOX::scaleX(float delta)
{
  _xPlus += delta;
  _xMinus -= delta;
}

void BOX::scaleY(float delta)
{
  _yPlus += delta;
  _yMinus -= delta;
}

void BOX::scaleZ(float delta)
{
  _zPlus += delta;
  _zMinus -= delta;
}
