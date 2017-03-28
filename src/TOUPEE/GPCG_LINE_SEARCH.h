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

#ifndef _GPCG_LINE_SEARCH_H_
#define _GPCG_LINE_SEARCH_H_

#include "LINE_SEARCH.h"

// TAO sure doesn't provide a lot of information on where this one came from.
class GPCG_LINE_SEARCH : public LINE_SEARCH
{
public:
  GPCG_LINE_SEARCH(void (*functionGradient)(const VECTOR&, Real&, VECTOR&));

  ////////////////////////////////////////////////////////////////////////////////////////
  // Originally:
  //    static PetscErrorCode TaoLineSearchApply_GPCG(TaoLineSearch ls,Vec x,PetscReal *f,
  //                                                  Vec g,Vec step_direction)
  // In file: src/linesearch/impls/gpcglinesearch/gpcglinesearch.c
  ////////////////////////////////////////////////////////////////////////////////////////
  virtual void LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s);

private:
  Real _maxstep;	     // maximum step size
  int _bracket;
  int _infoc;

  VECTOR _x;
  VECTOR _W1;
  VECTOR _W2;
  VECTOR _Gold;
};

#endif
