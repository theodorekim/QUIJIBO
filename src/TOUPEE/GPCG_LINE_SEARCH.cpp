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

#include "GPCG_LINE_SEARCH.h"

#include <iostream>
#include <cstdio>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
GPCG_LINE_SEARCH::GPCG_LINE_SEARCH(void (*functionGradient)(const VECTOR&, Real&, VECTOR&)) :
  LINE_SEARCH(functionGradient)
{
  _ftol		  = 0.05;
  _rtol		  = 0.0;
  _gtol		  = 0.0;
  _stepmin		= 1.0e-20;
  _stepmax		= 1.0e+20;
  _nfeval		  = 0; 
  _max_funcs	= 30;
  _step       = 1.0;

  _bracket		= 0; 
  _infoc      = 1;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void GPCG_LINE_SEARCH::LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s)
{
  //TAOLINESEARCH_GPCG_CTX *neP = (TAOLINESEARCH_GPCG_CTX *)ls->data;
  //PetscErrorCode  ierr;
  //PetscInt i;
  int i;
  //PetscBool g_computed=PETSC_FALSE; /* to prevent extra gradient computation */
  bool g_computed = false; // to prevent extra gradient computation
  //PetscReal d1,finit,actred,prered,rho, gdx;
  Real d1,finit,actred,prered,rho, gdx;

  //PetscFunctionBegin;
  // ls->stepmin - lower bound for step
  // ls->stepmax - upper bound for step
  // ls->rtol 	  - relative tolerance for an acceptable step
  // ls->ftol 	  - tolerance for sufficient decrease condition
  // ls->gtol 	  - tolerance for curvature condition
  // ls->nfeval	  - number of function evaluations
  // ls->nfeval	  - number of function/gradient evaluations
  // ls->max_funcs  - maximum number of function evaluations

  //ls->reason = TAOLINESEARCH_CONTINUE_ITERATING;
  _reason = LINESEARCH_CONTINUE_ITERATING;
  //ls->step = ls->initstep;
  _step = _initstep;
  //if (!neP->W2) 
  if (_W2.size() == 0) 
  {
    //ierr = VecDuplicate(x,&neP->W2); CHKERRQ(ierr);
    _W2 = x;
    //ierr = VecDuplicate(x,&neP->W1); CHKERRQ(ierr);
    _W1 = x;
    //ierr = VecDuplicate(x,&neP->Gold); CHKERRQ(ierr);
    _Gold = x;
    //neP->x = x;
    _x = x;
    //ierr = PetscObjectReference((PetscObject)neP->x); CHKERRQ(ierr);
  }
  // If X has changed, remake work vectors
  // TK: Just do the copy every time. If this becomes the bottleneck,
  // re-evaluate.
  //else if (x != neP->x) 
  else 
  {
    //ierr = VecDestroy(&neP->x); CHKERRQ(ierr);
    //ierr = VecDestroy(&neP->W1); CHKERRQ(ierr);
    //ierr = VecDestroy(&neP->W2); CHKERRQ(ierr);
    //ierr = VecDestroy(&neP->Gold); CHKERRQ(ierr);
    //ierr = VecDuplicate(x,&neP->W1); CHKERRQ(ierr);
    _W1 = x;
    //ierr = VecDuplicate(x,&neP->W2); CHKERRQ(ierr);
    _W2 = x;
    //ierr = VecDuplicate(x,&neP->Gold); CHKERRQ(ierr);
    _Gold = x;
    //ierr = PetscObjectDereference((PetscObject)neP->x); CHKERRQ(ierr);
    //neP->x = x;
    _x = x;
    //ierr = PetscObjectReference((PetscObject)neP->x); CHKERRQ(ierr);
  }

  //ierr = VecDot(g,s,&gdx); CHKERRQ(ierr);
  //gdx = g * s;
  gdx = g.dot(s);
  //ierr = VecCopy(x,neP->W2); CHKERRQ(ierr);
  _W2 = x;
  //ierr = VecCopy(g,neP->Gold); CHKERRQ(ierr);
  _Gold = g;
  //if (ls->bounded) 
  if (_bounded) 
  {
	  // Compute the smallest steplength that will make one nonbinding variable
	  // equal the bound
    //ierr = VecStepBoundInfo(x,ls->lower,ls->upper,s,&rho,&actred,&d1); CHKERRQ(ierr);
    VecStepBoundInfo(x,_lower,_upper,s,rho,actred,d1);
    //ls->step = PetscMin(ls->step,d1);
    _step = std::min(_step,d1);
  }
  rho=0; actred=0;

  //if (ls->step < 0) 
  if (_step < 0) 
  {
    //ierr = PetscInfo1(ls,"Line search error: initial step parameter %G < 0\n",ls->step); CHKERRQ(ierr);
    printf("Line search error: initial step parameter %G < 0\n", (double)_step);
    //ls->reason = TAOLINESEARCH_HALTED_OTHER;
    _reason = LINESEARCH_HALTED_OTHER;
    //PetscFunctionReturn(0);
    return;
  }

  // Initialization
  finit = f;
  //for (i=0; i< ls->max_funcs; i++) 
  for (i=0; i< _max_funcs; i++) 
  {
    // Force the step to be within the bounds
    //ls->step = PetscMax(ls->step,ls->stepmin);
    _step = std::max(_step,_stepmin);
    //ls->step = PetscMin(ls->step,ls->stepmax);
    _step = std::min(_step,_stepmax);

    //ierr = VecCopy(x,neP->W2); CHKERRQ(ierr);
    _W2 = x;
    //ierr = VecAXPY(neP->W2,ls->step,s); CHKERRQ(ierr);
    _W2 += _step * s;
    //if (ls->bounded) 
    if (_bounded) 
    	//ierr = VecMedian(ls->lower,neP->W2,ls->upper,neP->W2); CHKERRQ(ierr);
    	_W2 = median(_lower,_W2,_upper);

    // Gradient is not needed here.  Unless there is a separate
    //   gradient routine, compute it here anyway to prevent recomputing at
    //   the end of the line search
    //if (ls->hasobjective) 
    //{
    //  ierr = TaoLineSearchComputeObjective(ls,neP->W2,f); CHKERRQ(ierr);
    //  g_computed=PETSC_FALSE;
    //} 
    //else if (ls->usegts)
    //{
    //  ierr = TaoLineSearchComputeObjectiveAndGTS(ls,neP->W2,f,&gdx); CHKERRQ(ierr);
    //  g_computed=PETSC_FALSE;
    //} 
    //else 
    {
      // TK: Only supporting objective and gradient at the moment
      //ierr = TaoLineSearchComputeObjectiveAndGradient(ls,neP->W2,f,g); CHKERRQ(ierr);
      ComputeObjectiveAndGradient(_W2,f,g);
      //g_computed=PETSC_TRUE;
      g_computed = true;
    }

    if (0 == i) 
    	//ls->f_fullstep = *f;
    	_f_fullstep = f;

    //actred = *f - finit;
    actred = f - finit;
    //ierr = VecCopy(neP->W2,neP->W1); CHKERRQ(ierr);
    _W1 = _W2;
    //ierr = VecAXPY(neP->W1,-1.0,x); CHKERRQ(ierr);    /* W1 = W2 - X */
    _W1 += -1.0 * x;
    //ierr = VecDot(neP->W1,neP->Gold,&prered); CHKERRQ(ierr);
    //prered = _W1 * _Gold;
    prered = _W1.dot(_Gold);
    
    if (fabs(prered)<1.0e-100) prered=1.0e-12;
    rho = actred/prered;
    // If sufficient progress has been obtained, accept the
    // point.  Otherwise, backtrack. 

    //if (rho > ls->ftol){
    if (rho > _ftol)
    {
      break;
    }
    else
      //ls->step = (ls->step)/2;
      _step = (_step)/2;

    // Convergence testing
    //if (ls->step <= ls->stepmin || ls->step >= ls->stepmax)
    if (_step <= _stepmin || _step >= _stepmax)
    {
      //ls->reason = TAOLINESEARCH_HALTED_OTHER;
      _reason = LINESEARCH_HALTED_OTHER;
      //ierr = PetscInfo(ls,"Rounding errors may prevent further progress.  May not be a step satisfying\n"); CHKERRQ(ierr);
      printf("Rounding errors may prevent further progress.  May not be a step satisfying\n");
      //ierr = PetscInfo(ls,"sufficient decrease and curvature conditions. Tolerances may be too small.\n"); CHKERRQ(ierr);
      printf("sufficient decrease and curvature conditions. Tolerances may be too small.\n");

      break;
    }
    //if (ls->step == ls->stepmax) 
    if (_step == _stepmax) 
    {
      //ierr = PetscInfo1(ls,"Step is at the upper bound, stepmax (%G)\n",ls->stepmax); CHKERRQ(ierr);
      printf("Step is at the upper bound, stepmax (%G)\n",(double)_stepmax);
      //ls->reason = TAOLINESEARCH_HALTED_UPPERBOUND;
      _reason = LINESEARCH_HALTED_UPPERBOUND;
      break;
    }
    //if (ls->step == ls->stepmin) 
    if (_step == _stepmin) 
    {
      //ierr = PetscInfo1(ls,"Step is at the lower bound, stepmin (%G)\n",ls->stepmin); CHKERRQ(ierr);
      printf("Step is at the lower bound, stepmin (%G)\n",(double)_stepmin);
      //ls->reason = TAOLINESEARCH_HALTED_LOWERBOUND;
      _reason = LINESEARCH_HALTED_LOWERBOUND;
      break;
    }
    //if ((ls->nfeval+ls->nfgeval) >= ls->max_funcs) 
    if ((_nfeval+_nfgeval) >= _max_funcs) 
    {
      //ierr = PetscInfo2(ls,"Number of line search function evals (%D) > maximum (%D)\n",ls->nfeval+ls->nfgeval,ls->max_funcs); CHKERRQ(ierr);
      printf("Number of line search function evals (%i) > maximum (%i)\n",_nfeval+_nfgeval,_max_funcs);
      //ls->reason = TAOLINESEARCH_HALTED_MAXFCN;
      _reason = LINESEARCH_HALTED_MAXFCN;
      break;
    }
    //if ((neP->bracket) && (ls->stepmax - ls->stepmin <= ls->rtol*ls->stepmax))
    if ((_bracket) && (_stepmax - _stepmin <= _rtol*_stepmax))
    {
      //ierr = PetscInfo1(ls,"Relative width of interval of uncertainty is at most rtol (%G)\n",ls->rtol); CHKERRQ(ierr);
      printf("Relative width of interval of uncertainty is at most rtol (%G)\n",(double)_rtol);
      //ls->reason = TAOLINESEARCH_HALTED_RTOL;
      _reason = LINESEARCH_HALTED_RTOL;
    	break;
    }
  }
  //ierr = PetscInfo2(ls,"%D function evals in line search, step = %G\n",ls->nfeval+ls->nfgeval,ls->step); CHKERRQ(ierr);
  printf("%i function evals in line search, step = %G\n",_nfeval+_nfgeval,(double)_step);

  // TK: This will never happen currently, as only objective + gradient is supported.
  // set new solution vector and compute gradient if necessary
  //ierr = VecCopy(neP->W2, x); CHKERRQ(ierr);
  //if (!g_computed) {
  //  ierr = TaoLineSearchComputeGradient(ls,x,g); CHKERRQ(ierr);
  //}
  //PetscFunctionReturn(0);
}
