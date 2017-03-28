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

#ifndef _TAO_SOLVER_H
#define _TAO_SOLVER_H

#include "LINE_SEARCH.h"

class TAO_SOLVER 
{
public:
  TAO_SOLVER();

  // solution - contains the initial guess, and will contain the final solution
  // upperBound - the max values "solution" should be allowed to take
  // lowerBound - the min values "solution" should be allowed to take
  // functionGradient - pointer to the function being optimized
  void optimizeBLMVM(VECTOR& solution, VECTOR& upperBound, VECTOR& lowerBound, 
                     void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient),
                     void (*monitor)() = NULL);
  
  // solution - contains the initial guess, and will contain the final solution
  // functionGradient - pointer to the function being optimized
  void optimizeLMVM(VECTOR& solution, void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient));

  // solution - contains the initial guess, and will contain the final solution
  // functionGradient - pointer to the function being optimized
  // buildHessian - pointer to the Hessian of the function being optimized
  void optimizeNLS(VECTOR& solution, 
                   void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient),
                   void (*buildHessian)(const VECTOR& state, MATRIX& hessian));

  // solution - contains the initial guess, and will contain the final solution
  // functionGradient - pointer to the function being optimized
  // buildHessian - pointer to the Hessian of the function being optimized
  //
  // This is the "verbose" version with all the quasi-trust-region cruft from TAO
  // still inlined for debugging purposes. The easier to read version is "optimizeNLS"
  void verboseOptimizeNLS(VECTOR& solution, 
                          void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient),
                          void (*buildHessian)(const VECTOR& state, MATRIX& hessian));

  enum SOLVER_STATE 
  {
    CONVERGED_FATOL          =  1, // f(X)-f(X*) <= fatol
    CONVERGED_FRTOL          =  2, // |F(X) - f(X*)|/|f(X)| < frtol
    CONVERGED_GATOL          =  3, // ||g(X)|| < gatol
    CONVERGED_GRTOL          =  4, // ||g(X)|| / f(X)  < grtol
    CONVERGED_GTTOL          =  5, // ||g(X)|| / ||g(X0)|| < gttol
    CONVERGED_STEPTOL        =  6, // step size small
    CONVERGED_MINF          =  7, // F < F_min
    CONVERGED_USER          =  8, // User defined
    // diverged
    DIVERGED_MAXITS         = -2,
    DIVERGED_NAN            = -4,
    DIVERGED_MAXFCN         = -5,
    DIVERGED_LS_FAILURE     = -6,
    DIVERGED_TR_REDUCTION   = -7,
    DIVERGED_USER           = -8, // User defined
    // keep going
    CONTINUE_ITERATING      =  0
  };
  enum LINE_SEARCH_TYPE
  {
    USING_GPCG,
    USING_MORE_THUENTE,
    USING_UNIT,
    USING_ARMIJO
  };
  
  enum KSPConvergedReason {
    // converged
    KSP_CONVERGED_RTOL_NORMAL        =  1,
    KSP_CONVERGED_ATOL_NORMAL        =  9,
    KSP_CONVERGED_RTOL               =  2,
    KSP_CONVERGED_ATOL               =  3,
    KSP_CONVERGED_ITS                =  4,
    KSP_CONVERGED_CG_NEG_CURVE       =  5,
    KSP_CONVERGED_CG_CONSTRAINED     =  6,
    KSP_CONVERGED_STEP_LENGTH        =  7,
    KSP_CONVERGED_HAPPY_BREAKDOWN    =  8,
    // diverged
    KSP_DIVERGED_NULL                = -2,
    KSP_DIVERGED_ITS                 = -3,
    KSP_DIVERGED_DTOL                = -4,
    KSP_DIVERGED_BREAKDOWN           = -5,
    KSP_DIVERGED_BREAKDOWN_BICG      = -6,
    KSP_DIVERGED_NONSYMMETRIC        = -7,
    KSP_DIVERGED_INDEFINITE_PC       = -8,
    KSP_DIVERGED_NAN                 = -9,
    KSP_DIVERGED_INDEFINITE_MAT      = -10,
   
    KSP_CONVERGED_ITERATING          =  0
  };

  void DefaultConvergenceTest(const int functionEvals, const int functionGradEvals, SOLVER_STATE& reason);
  VECTOR median(const VECTOR& v0, const VECTOR& v1, const VECTOR& v2);
  void SetTolerances(Real fatol, Real frtol, Real gatol, Real grtol, Real gttol);
  Real& maxScore() { return _maxScore; };

  // set which line search you want
  void useMoreThuente()    { _whichLineSearch = USING_MORE_THUENTE; };
  void useUnitLineSearch() { _whichLineSearch = USING_UNIT; };
  void useArmijo()         { _whichLineSearch = USING_ARMIJO; };
  void useGPCG()           { _whichLineSearch = USING_GPCG; };

private:
  int _max_it;
  int _max_funcs;
  Real _fatol;
  Real _frtol;
  Real _gatol;
  Real _grtol;
  Real _gttol;
  Real _catol;
  Real _crtol;
  Real _xtol;
  Real _steptol;
  Real _fmin;
  int _hist_max;
  int _hist_len;
  Real _residual;
  int _niter;
  Real _gnorm0;
  Real _fc;
  Real _step;
  Real _cnorm;
  Real _cnorm0;
  Real _maxScore;

  LINE_SEARCH_TYPE _whichLineSearch;
};

#endif
