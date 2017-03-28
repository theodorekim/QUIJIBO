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
#include "SPHERE.h"
#include "FIELD_3D.h"
#include "TIMER.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SPHERE::SPHERE(const VEC3F& center, const float radius) :
  _radius(radius)
{
  _translation = center;
  _rotation = MATRIX3::I();
}

SPHERE::SPHERE(float x, float y, float z, float radius) :
  _radius(radius)
{
  _translation[0] = x;
  _translation[1] = y;
  _translation[2] = z;
  _rotation = MATRIX3::I();
}

SPHERE::SPHERE(const VEC3F& translation, const MATRIX3& rotation, float radius) :
  _radius(radius)
{
  _translation = translation;
  _rotation = rotation;
}

SPHERE::~SPHERE()
{

}

bool SPHERE::inside(float x, float y, float z)
{
  float translate[] = {x - _translation[0], y - _translation[1], z - _translation[2]};
  float magnitude = translate[0] * translate[0] + 
                    translate[1] * translate[1] + 
                    translate[2] * translate[2];

  return (magnitude < _radius * _radius);
}

void SPHERE::scale(const Real alpha)
{
  _translation *= alpha;
  _radius *= alpha;
}

void SPHERE::boundingBox(VEC3F& mins, VEC3F& maxs) const
{
  mins = _translation;
  maxs= _translation;

  // add just a little bit of padding
  mins -= _radius;
  maxs += _radius;
}

Real SPHERE::distance(const VEC3F& point) const
{
  VEC3F diff = point - _translation;

  Real centerDistance = norm(diff);

  return centerDistance - _radius;
}

///////////////////////////////////////////////////////////////////////
// get the cells that are inside this sphere
///////////////////////////////////////////////////////////////////////
void SPHERE::cellsInside(const FIELD_3D& field, vector<int>& inside) const
{
  TIMER functionTimer(__FUNCTION__);
  const Real dx = field.dx();
  VEC3I centerIndices;
  field.cellIndex(_translation, centerIndices);
  int radiusCells = _radius / dx + 1;

  VEC3I startIndices = centerIndices;
  VEC3I endIndices = centerIndices;
  startIndices -= radiusCells;
  endIndices += radiusCells;

  startIndices[0] = startIndices[0] < 0 ? 0 : startIndices[0];
  startIndices[1] = startIndices[1] < 0 ? 0 : startIndices[1];
  startIndices[2] = startIndices[2] < 0 ? 0 : startIndices[2];

  int xRes = field.xRes();
  int yRes = field.yRes();
  int zRes = field.zRes();
  endIndices[0] = endIndices[0] > xRes ? xRes : endIndices[0];
  endIndices[1] = endIndices[1] > yRes ? yRes : endIndices[1];
  endIndices[2] = endIndices[2] > zRes ? zRes : endIndices[2];

  inside.clear();
  //Real radiusSq = _radius * _radius;

  // fatten the radius to check against
  Real radiusSq = (_radius + dx) * (_radius + dx);

  TIMER functionTimer2("cellsInside grid iteration");
  for (int z = startIndices[2]; z < endIndices[2]; z++)
    for (int y = startIndices[1]; y < endIndices[1]; y++)
      for (int x = startIndices[0]; x < endIndices[0]; x++)
      {
        VEC3F diff = field.cellCenter(x,y,z) - _translation;
        if (diff * diff <= radiusSq)
        {
          int index = x + y * xRes + z * xRes * yRes;
          inside.push_back(index);
        }
      }
}

///////////////////////////////////////////////////////////////////////
// get the cells that are inside this sphere
///////////////////////////////////////////////////////////////////////
void SPHERE::smartCellsInside(const FIELD_3D& field, vector<int>& inside) const
{
  TIMER functionTimer(__FUNCTION__);
  const Real dx = field.dx();
  VEC3I centerIndices;
  field.cellIndex(_translation, centerIndices);
  int radiusCells = _radius / dx + 1;

  VEC3I startIndices = centerIndices;
  VEC3I endIndices = centerIndices;
  startIndices -= radiusCells;
  endIndices += radiusCells;

  startIndices[0] = startIndices[0] < 0 ? 0 : startIndices[0];
  startIndices[1] = startIndices[1] < 0 ? 0 : startIndices[1];
  startIndices[2] = startIndices[2] < 0 ? 0 : startIndices[2];

  int xRes = field.xRes();
  int yRes = field.yRes();
  int zRes = field.zRes();
  endIndices[0] = endIndices[0] > xRes ? xRes : endIndices[0];
  endIndices[1] = endIndices[1] > yRes ? yRes : endIndices[1];
  endIndices[2] = endIndices[2] > zRes ? zRes : endIndices[2];

  map<int, bool> seenAlready;
  for (unsigned int x = 0; x < inside.size(); x++)
    seenAlready[inside[x]] = true;

  vector<int> toCheck;
  for (int z = startIndices[2]; z < endIndices[2]; z++)
    for (int y = startIndices[1]; y < endIndices[1]; y++)
      for (int x = startIndices[0]; x < endIndices[0]; x++)
      {
        int index = x + y * xRes + z * xRes * yRes;
        if (seenAlready.find(index) == seenAlready.end())
          toCheck.push_back(index);
      }

  inside.clear();
  const Real radiusSq = _radius * _radius;
  for (unsigned int x = 0; x < toCheck.size(); x++)
  {
    const int& index = toCheck[x];
    VEC3F diff = field.cellCenter(index) - _translation;
    if (diff * diff < radiusSq)
      inside.push_back(index);
  }
}

///////////////////////////////////////////////////////////////////////
// does this sphere overlap another one?
///////////////////////////////////////////////////////////////////////
bool SPHERE::overlap(const SPHERE& compare) const
{
  VEC3F diff = compare._translation - _translation;

  float distance = diff.magnitude();

  return (distance < _radius + compare._radius);
}

///////////////////////////////////////////////////////////////////////
// does this sphere overlap another one?
///////////////////////////////////////////////////////////////////////
float SPHERE::distance(const SPHERE& compare) const
{
  VEC3F diff = compare._translation - _translation;
  return diff.magnitude();
}
