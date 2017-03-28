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

// This is a port of the Toolkit for Advanced Optimiation (TAO):
//
//   http://www.mcs.anl.gov/research/projects/tao/
//
// but with all the PETSc dependnecies removed.
// TAO Without PETSc -> T W/O Pe -> TOUPEE

#ifndef _MORE_THUENTE_H_
#define _MORE_THUENTE_H_

#include "LINE_SEARCH.h"

class MORE_THUENTE : public LINE_SEARCH
{
public:
  MORE_THUENTE(void (*functionGradient)(const VECTOR&, Real&, VECTOR&));
  virtual void LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s);
  void mcstep(Real& stx, Real& fx, Real& dx,
              Real& sty, Real& fy, Real& dy,
              Real& stp, Real& fp, Real& dp);

private:
  // these appear to be More-Thuente-specific
  int _bracket;
  int _infoc;
  VECTOR _x;
  VECTOR _work;
};

#endif
