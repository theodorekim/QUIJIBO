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
#ifndef SPHERE_H
#define SPHERE_H

#include "VEC3.h"
#include "MATRIX3.h"
#include "FIELD_3D.h"

class SPHERE
{
public:
	SPHERE(const VEC3F& center, const float radius);
	SPHERE(float x, float y, float z, float radius);
	SPHERE(const VEC3F& translation, const MATRIX3& rotation, float radius);
	virtual ~SPHERE();

  bool inside(float x, float y, float z);
  void draw();
  void scale(const Real alpha);
  float& radius() { return _radius; };
  const float radius() const { return _radius; };
  const VEC3F& center() const { return _translation; };
  const VEC3F& translation() const { return _translation; };
  void boundingBox(VEC3F& mins, VEC3F& maxs) const;

  Real distance(const VEC3F& point) const;

  // get the cells that are inside this sphere
  void cellsInside(const FIELD_3D& field, vector<int>& inside) const;
 
  // get the cells that are inside this sphere, but re-use the ones already in "inside"
  void smartCellsInside(const FIELD_3D& field, vector<int>& inside) const;

  // does this sphere overlap another one?
  bool overlap(const SPHERE& compare) const;

  // what's the distance between two spheres?
  float distance(const SPHERE& compare) const;

private:
  float _radius;
  VEC3F _translation;
  MATRIX3 _rotation;
};

#endif
