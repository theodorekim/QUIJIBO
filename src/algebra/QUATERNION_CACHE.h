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
#ifndef QUATERNION_CACHE_H
#define QUATERNION_CACHE_H

#include "QUATERNION.h"

// A cache class for arrays of quaternion values
class QUATERNION_CACHE {
public:
  QUATERNION_CACHE() : 
    _xRes(-1), _yRes(-1), _zRes(-1), _slabSize(-1), _totalCells(-1)
  {
  }

  QUATERNION_CACHE(const int& xRes, const int& yRes, const int& zRes) : 
    _xRes(xRes), _yRes(yRes), _zRes(zRes), _slabSize(xRes * yRes),
    _totalCells(xRes * yRes * zRes)
  {
    _data.resize(_totalCells);
  }

  const vector<QUATERNION>& operator()(int x, int y, int z) const {
    return _data[x + y * _xRes + z * _slabSize];
  };
  vector<QUATERNION>& operator()(int x, int y, int z) {
    return _data[x + y * _xRes + z * _slabSize];
  };

  const int xRes() { return _xRes; };
  const int yRes() { return _yRes; };
  const int zRes() { return _zRes; };
  const int totalCells() { return _totalCells; };

private:
  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;
  int _totalCells;

  vector<vector<QUATERNION> > _data;
};

#endif
