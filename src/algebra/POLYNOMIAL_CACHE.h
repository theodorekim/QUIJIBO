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
#ifndef POLYNOMIAL_CACHE_H
#define POLYNOMIAL_CACHE_H

#include "QUATERNION.h"

// A cache class for partially evaluated quaternionic polynomials.
// Avoids redoing all the power operations that were performed when
// computing the initial fractal
class POLYNOMIAL_CACHE {
public:
  POLYNOMIAL_CACHE() : 
    _xRes(-1), _yRes(-1), _zRes(-1), _slabSize(-1), _totalCells(-1)
  {
  }

  POLYNOMIAL_CACHE(const int& xRes, const int& yRes, const int& zRes) : 
    _xRes(xRes), _yRes(yRes), _zRes(zRes), _slabSize(xRes * yRes),
    _totalCells(xRes * yRes * zRes)
  {
    _forwards.resize(_totalCells);
    _backwards.resize(_totalCells);
  }

  const vector<QUATERNION>& forward(int x, int y, int z) const {
    return _forwards[x + y * _xRes + z * _slabSize];
  };
  vector<QUATERNION>& forward(int x, int y, int z) {
    return _forwards[x + y * _xRes + z * _slabSize];
  };
  const vector<QUATERNION>& backward(int x, int y, int z) const {
    return _backwards[x + y * _xRes + z * _slabSize];
  };
  vector<QUATERNION>& backward(int x, int y, int z) {
    return _backwards[x + y * _xRes + z * _slabSize];
  };

  const vector<vector<QUATERNION> >& forwards() const { return _forwards; };
  const vector<vector<QUATERNION> >& backwards() const { return _backwards; };

  VEC3F dims() const { return VEC3F(_xRes, _yRes, _zRes); };

private:
  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;
  int _totalCells;

  vector<vector<QUATERNION> > _forwards;
  vector<vector<QUATERNION> > _backwards;
};

#endif
