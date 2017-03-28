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

#ifndef _UNIT_LINE_SEARCH_H_
#define _UNIT_LINE_SEARCH_H_

#include "LINE_SEARCH.h"

class UNIT_LINE_SEARCH : public LINE_SEARCH
{
public:
  UNIT_LINE_SEARCH(void (*functionGradient)(const VECTOR&, Real&, VECTOR&)) : 
    LINE_SEARCH(functionGradient)
  {
  };

  ////////////////////////////////////////////////////////////////////////////////////////
  // Originally:
  //    static PetscErrorCode TaoLineSearchApply_Unit(TaoLineSearch ls,Vec x,PetscReal *f,
  //                                                  Vec g,Vec step_direction)
  // In file: src/linesearch/impls/unit/unit.c
  ////////////////////////////////////////////////////////////////////////////////////////
  virtual void LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s)
  {
    Real ftry;
    //PetscReal startf = *f;
    Real startf = f;

    //PetscFunctionBegin;
    
    // Take unit step (newx = startx + 1.0*step_direction)
    //VecAXPY(x,1.0,step_direction);
    x += s;

    //TaoLineSearchComputeObjectiveAndGradient(ls,x,&ftry,g);
    ComputeObjectiveAndGradient(x,ftry,g);

    if (startf < ftry)
    {
      //printf("Tao Apply Unit Step, FINCREASE: F old:= %12.10e, F new: %12.10e\n",startf,ftry);
    }
    //*f = ftry;
    f = ftry;
    //ls->step = 1.0;
    _step = 1.0;
    //ls->reason=TAOLINESEARCH_SUCCESS;
    _reason = LINESEARCH_SUCCESS;
  };
};

#endif
