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
// BOX.h: interface for the BOX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BOX_H
#define BOX_H

#include "TRIANGLE.h"

class BOX
{
public:
	BOX(float xPlus = 0.55f, float xMinus = 0.45f, 
      float yPlus = 0.55f, float yMinus = 0.45f, 
      float zPlus = 0.55f, float zMinus = 0.45f);
  BOX(const VEC3F& center, const VEC3F& sizes);
	virtual ~BOX();

  virtual void draw();

  void translateX(float fraction);
  void translateY(float fraction);
  void translateZ(float fraction);

  void scaleX(float scale);
  void scaleY(float scale);
  void scaleZ(float scale);
  void scale(float scale)
  {
    scaleX(scale);
    scaleY(scale);
    scaleZ(scale);
  }

  VEC3F boxMins() { return VEC3F(_xMinus, _yMinus, _zMinus); };
  VEC3F boxMaxs() { return VEC3F(_xPlus, _yPlus, _zPlus); };

  VEC3F center() { return VEC3F((_xPlus + _xMinus) * 0.5,
                                (_yPlus + _yMinus) * 0.5,
                                (_zPlus + _zMinus) * 0.5); };
  VEC3F lengths() { return VEC3F((_xPlus - _xMinus),
                                 (_yPlus - _yMinus),
                                 (_zPlus - _zMinus)); };

  void setCenter(const VEC3F& newCenter);

private:
  float _xPlus, _xMinus;
  float _yPlus, _yMinus;
  float _zPlus, _zMinus;

  inline float mymin(float i, float j) { return (i < j) ? i : j; };
};

#endif
