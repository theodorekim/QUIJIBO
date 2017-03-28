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

#ifndef _LINE_SEARCH_H_
#define _LINE_SEARCH_H_

#include "SETTINGS.h"

class LINE_SEARCH
{
public:
  LINE_SEARCH(void (*functionGradient)(const VECTOR&, Real&, VECTOR&));
  void VecStepBoundInfo(VECTOR& x, VECTOR& xl, VECTOR& xu, VECTOR& dx, Real& boundmin, Real& wolfemin, Real& boundmax);
  void ComputeObjectiveAndGradient(VECTOR& x, Real& f, VECTOR& g);
  void VecBoundGradientProjection(const VECTOR& G, const VECTOR& X, const VECTOR& XL, const VECTOR& XU, VECTOR& GP);
  VECTOR median(const VECTOR& v0, const VECTOR& v1, const VECTOR& v2);
  void LineSearchApply(VECTOR& x, Real& f, VECTOR& g, VECTOR& s);
  virtual void LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s) = 0;

  enum LINE_SEARCH_TERMINATION_REASON
  {
    LINESEARCH_FAILED_INFORNAN = -1,
    LINESEARCH_FAILED_BADPARAMETER = -2,
    LINESEARCH_FAILED_ASCENT = -3,
    LINESEARCH_CONTINUE_ITERATING = 0,
    LINESEARCH_SUCCESS = 1,
    LINESEARCH_SUCCESS_USER = 2,
    LINESEARCH_HALTED_OTHER = 3,
    LINESEARCH_HALTED_MAXFCN = 4,
    LINESEARCH_HALTED_UPPERBOUND = 5,
    LINESEARCH_HALTED_LOWERBOUND = 6,
    LINESEARCH_HALTED_RTOL = 7,
    LINESEARCH_HALTED_USER = 8
  };

public:
  //PETSCHEADER(struct _TaoLineSearchOps);
  //void *userctx_func;
  //void *userctx_grad;
  //void *userctx_funcgrad;
  //void *userctx_funcgts;
  
  bool _setupcalled;
  bool _usegts;
  bool _usetaoroutines;
  bool _hasobjective;
  bool _hasgradient;
  bool _hasobjectiveandgradient;
  void *_data;

  // bounds used for some line searches
  VECTOR _lower;
  VECTOR _upper;
  int _bounded;

  VECTOR _start_x;
  VECTOR _stepdirection;
  Real _f_fullstep;
  Real _new_f;
  VECTOR _new_x;
  VECTOR _new_g;

  Real _step;
  Real _initstep;

  int _max_funcs;
  int _nfeval;  // function evaluations
  int _ngeval;  // gradient evaluations
  int _nfgeval; // function-grad evaluations
  LINE_SEARCH_TERMINATION_REASON _reason;

  Real _rtol;	 /* relative tol for acceptable step (rtol>0) */
  Real _ftol;	 /* tol for sufficient decr. condition (ftol>0) */
  Real _gtol;	 /* tol for curvature condition (gtol>0)*/
  Real _stepmin;	 /* lower bound for step */
  Real _stepmax;	 /* upper bound for step */
  bool _viewls;    /* print out information if true */

  // the objective function to minimize
  void (*_functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient);
};

#endif
