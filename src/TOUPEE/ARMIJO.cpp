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

#include "ARMIJO.h"
#include <cmath>
#include <cstdio>

using namespace std;

#define REPLACE_FIFO 1
#define REPLACE_MRU  2

#define REFERENCE_MAX  1
#define REFERENCE_AVE  2
#define REFERENCE_MEAN 3

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
ARMIJO::ARMIJO(void (*functionGradient)(const VECTOR&, Real&, VECTOR&)) : 
  LINE_SEARCH(functionGradient)
{
  _memory = VECTOR(1);
  _alpha = 1.0;
  _beta = 0.5;
  _beta_inf = 0.5;
  _sigma = 1e-4;
  _memorySize = 1;
  _referencePolicy = REFERENCE_MAX;
  _replacementPolicy = REPLACE_MRU;
  _nondescending = false;

  _memorySetup = false;

  _initstep = 1.0;
}

////////////////////////////////////////////////////////////////////////////////////////
// Originally:
//    static PetscErrorCode TaoLineSearchApply_Armijo(TaoLineSearch ls,Vec x,PetscReal *f,
//                                                  Vec g,Vec step_direction)
// In file: src/linesearch/impls/armijo/armijo.c
////////////////////////////////////////////////////////////////////////////////////////
/* @ TaoApply_Armijo - This routine performs a linesearch. It
     backtracks until the (nonmonotone) Armijo conditions are satisfied.

     Input Parameters:
  +  tao - TaoSolver context
  .  X - current iterate (on output X contains new iterate, X + step*S)
  .  S - search direction
  .  f - merit function evaluated at X
  .  G - gradient of merit function evaluated at X
  .  W - work vector
  -  step - initial estimate of step length

     Output parameters:
  +  f - merit function evaluated at new iterate, X + step*S
  .  G - gradient of merit function evaluated at new iterate, X + step*S
  .  X - new iterate
  -  step - final step length

     Info is set to one of:
  .   0 - the line search succeeds; the sufficient decrease
     condition and the directional derivative condition hold

     negative number if an input parameter is invalid
  -   -1 -  step < 0 

     positive number > 1 if the line search otherwise terminates
  +    1 -  Step is at the lower bound, stepmin.
@ */
void ARMIJO::LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s)
{
  //TAOLINESEARCH_ARMIJO_CTX *armP = (TAOLINESEARCH_ARMIJO_CTX *)ls->data;
  //PetscErrorCode ierr;
  //PetscInt i;
  int i;

  //PetscReal fact, ref, gdx;
  Real fact, ref, gdx;
  
  //PetscInt idx;
  int idx;
  
  //PetscBool g_computed=PETSC_FALSE;
  bool g_computed = false;

  //PetscFunctionBegin;

  //ls->reason = TAOLINESEARCH_CONTINUE_ITERATING;
  _reason = LINESEARCH_CONTINUE_ITERATING;

  //if (!armP->work) 
  if (_work.size() != x.size())
  {
    //ierr = VecDuplicate(x,&armP->work); CHKERRQ(ierr);
    _work = x;
    //armP->x = x;
    _x = x;
    //ierr = PetscObjectReference((PetscObject)armP->x); CHKERRQ(ierr);
  }
  // If x has changed, then recreate work
  //else if (x != armP->x) 
  else if ((_x - x).norm() < 1e-8) 
  {
    //ierr = VecDestroy(&armP->work); CHKERRQ(ierr);
    //ierr = VecDuplicate(x,&armP->work); CHKERRQ(ierr);
    _work = x;
    //ierr = PetscObjectDereference((PetscObject)armP->x); CHKERRQ(ierr);
    //armP->x = x;
    _x = x;
    //ierr = PetscObjectReference((PetscObject)armP->x); CHKERRQ(ierr);
  }

  /* Check linesearch parameters */
  //if (armP->alpha < 1) 
  if (_alpha < 1) 
  {
    //ierr = PetscInfo1(ls,"Armijo line search error: alpha (%G) < 1\n", armP->alpha); CHKERRQ(ierr);
    printf("Armijo line search error: alpha (%G) < 1\n", (double)_alpha);
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  } 
  //else if ((armP->beta <= 0) || (armP->beta >= 1)) 
  else if ((_beta <= 0) || (_beta >= 1)) 
  {
    //ierr = PetscInfo1(ls,"Armijo line search error: beta (%G) invalid\n", armP->beta); CHKERRQ(ierr);
    printf("Armijo line search error: beta (%G) invalid\n", (double)_beta);
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  } 
  //else if ((armP->beta_inf <= 0) || (armP->beta_inf >= 1)) 
  else if ((_beta_inf <= 0) || (_beta_inf >= 1)) 
  {
    //ierr = PetscInfo1(ls,"Armijo line search error: beta_inf (%G) invalid\n", armP->beta_inf); CHKERRQ(ierr);
    printf("Armijo line search error: beta_inf (%G) invalid\n", (double)_beta_inf);
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  } 
  //else if ((armP->sigma <= 0) || (armP->sigma >= 0.5)) 
  else if ((_sigma <= 0) || (_sigma >= 0.5)) 
  {
    //ierr = PetscInfo1(ls,"Armijo line search error: sigma (%G) invalid\n", armP->sigma); CHKERRQ(ierr);
    printf("Armijo line search error: sigma (%G) invalid\n", (double)_sigma);
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  } 
  //else if (armP->memorySize < 1) 
  else if (_memorySize < 1) 
  {
    //ierr = PetscInfo1(ls,"Armijo line search error: memory_size (%D) < 1\n", armP->memorySize); CHKERRQ(ierr);
    printf("Armijo line search error: memory_size (%i) < 1\n", _memorySize);
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  } 
  //else if ((armP->referencePolicy != REFERENCE_MAX) &&
  //         (armP->referencePolicy != REFERENCE_AVE) &&
  //         (armP->referencePolicy != REFERENCE_MEAN)) 
  else if ((_referencePolicy != REFERENCE_MAX) &&
           (_referencePolicy != REFERENCE_AVE) &&
           (_referencePolicy != REFERENCE_MEAN)) 
  {
    //ierr = PetscInfo(ls,"Armijo line search error: reference_policy invalid\n"); CHKERRQ(ierr);
    printf("Armijo line search error: reference_policy invalid\n");
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  } 
  //else if ((armP->replacementPolicy != REPLACE_FIFO) && 
  //    (armP->replacementPolicy != REPLACE_MRU)) 
  else if ((_replacementPolicy != REPLACE_FIFO) && 
           (_replacementPolicy != REPLACE_MRU)) 
  {
    //ierr = PetscInfo(ls,"Armijo line search error: replacement_policy invalid\n"); CHKERRQ(ierr);
    printf("Armijo line search error: replacement_policy invalid\n");
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }
  //else if (PetscIsInfOrNanReal(*f)) 
  else if (isinf(f) || isnan(f)) 
  {
    //ierr = PetscInfo(ls,"Armijo line search error: initial function inf or nan\n"); CHKERRQ(ierr);
    printf("Armijo line search error: initial function inf or nan\n");
    //ls->reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }
  //if (ls->reason != TAOLINESEARCH_CONTINUE_ITERATING) 
  if (_reason != LINESEARCH_CONTINUE_ITERATING) 
  {
    //PetscFunctionReturn(0);
    return;
  }

  // Check to see of the memory has been allocated.  If not, allocate
  //   the historical array and populate it with the initial function
  //   values.
  //if (armP->memory == PETSC_NULL) 
  if (_memory.size() != _memorySize)
  {
    //ierr = PetscMalloc(sizeof(PetscReal)*armP->memorySize, &armP->memory ); CHKERRQ(ierr);
    //_memory.resizeAndWipe(_memorySize);
    _memory.resize(_memorySize);
    _memory.setZero();
  }

  //if (!armP->memorySetup) 
  if (!_memorySetup) 
  {
    //for (i = 0; i < armP->memorySize; i++)
    for (i = 0; i < _memorySize; i++)
      //armP->memory[i] = armP->alpha*(*f);
      _memory[i] = _alpha * f;

    //armP->current = 0;
    _current = 0;
    //armP->lastReference = armP->memory[0];
    _lastReference = _memory[0];
    //armP->memorySetup=PETSC_TRUE;
    _memorySetup = true;
  }
  
  // Calculate reference value (MAX)
  //ref = armP->memory[0];
  ref = _memory[0];
  idx = 0;

  //for (i = 1; i < armP->memorySize; i++) 
  for (i = 1; i < _memorySize; i++) 
  {
    //if (armP->memory[i] > ref)
    if (_memory[i] > ref) 
    {
      //ref = armP->memory[i];
      ref = _memory[i];
      idx = i;
    }
  }

  //if (armP->referencePolicy == REFERENCE_AVE) 
  if (_referencePolicy == REFERENCE_AVE) 
  {
    ref = 0;
    //for (i = 0; i < armP->memorySize; i++)
    for (i = 0; i < _memorySize; i++)
      //ref += armP->memory[i];
      ref += _memory[i];
    //ref = ref / armP->memorySize;
    ref = ref / _memorySize;
    //ref = PetscMax(ref, armP->memory[armP->current]);
    ref = std::max((Real)ref, (Real)_memory[_current]);
  } 
  //else if (armP->referencePolicy == REFERENCE_MEAN) 
  else if (_referencePolicy == REFERENCE_MEAN) 
    //ref = PetscMin(ref, 0.5*(armP->lastReference + armP->memory[armP->current]));
    ref = std::min(ref, (Real)(0.5*(_lastReference + _memory[_current])));
  //ierr = VecDot(g,s,&gdx); CHKERRQ(ierr);
  gdx = g.dot(s);

  //if (PetscIsInfOrNanReal(gdx)) 
  if (isinf(gdx) || isnan(gdx)) 
  {
    //ierr = PetscInfo1(ls,"Initial Line Search step * g is Inf or Nan (%G)\n",gdx); CHKERRQ(ierr);
    printf("Initial Line Search step * g is Inf or Nan (%G)\n",(double)gdx);
    //ls->reason=TAOLINESEARCH_FAILED_INFORNAN;
    _reason = LINESEARCH_FAILED_INFORNAN;
    //PetscFunctionReturn(0);
    return;
  }
  if (gdx >= 0.0) 
  {
    //ierr = PetscInfo1(ls,"Initial Line Search step is not descent direction (g's=%G)\n",gdx); CHKERRQ(ierr);
    printf("Initial Line Search step is not descent direction (g's=%G)\n",(double)gdx);
    //ls->reason = TAOLINESEARCH_FAILED_ASCENT;
    _reason = LINESEARCH_FAILED_ASCENT;
    //PetscFunctionReturn(0);
    return;
  }
  
  //if (armP->nondescending) 
  if (_nondescending) 
    //fact = armP->sigma; 
    fact = _sigma; 
  else 
    //fact = armP->sigma * gdx;
    fact = _sigma * gdx;
  //ls->step = ls->initstep;
  _step = _initstep;
  //while (ls->step >= ls->stepmin && (ls->nfeval+ls->nfgeval) < ls->max_funcs) 
  while (_step >= _stepmin && (_nfeval + _nfgeval) < _max_funcs) 
  {
    // Calculate iterate
    //ierr = VecCopy(x,armP->work); CHKERRQ(ierr);
    _work = x;
    //ierr = VecAXPY(armP->work,ls->step,s); CHKERRQ(ierr);
    _work += _step * s;
    //if (ls->bounded) 
    if (_bounded) 
    {
      //ierr = VecMedian(ls->lower,armP->work,ls->upper,armP->work); CHKERRQ(ierr);
      _work = median(_lower, _work, _upper);
    }

    // TK: only supporting objective + gradient at the moment
    /*
    // Calculate function at new iterate
    if (ls->hasobjective) 
    {
      ierr = TaoLineSearchComputeObjective(ls,armP->work,f); CHKERRQ(ierr);
      //_nfeval++;
      g_computed=PETSC_FALSE;
    } 
    else if (ls->usegts) 
    {
      ierr = TaoLineSearchComputeObjectiveAndGTS(ls,armP->work,f,&gdx); CHKERRQ(ierr);
      //_nfeval++;
      g_computed=PETSC_FALSE;
    } 
    else 
    */
    {
      //ierr = TaoLineSearchComputeObjectiveAndGradient(ls,armP->work,f,g); CHKERRQ(ierr);
      ComputeObjectiveAndGradient(_work,f,g);
      _nfgeval++;

      //g_computed=PETSC_TRUE;
      g_computed = true;
    }

    //if (ls->step == ls->initstep) 
    if (_step == _initstep) 
      //ls->f_fullstep = *f;
      _f_fullstep = f;

    //if (PetscIsInfOrNanReal(*f))
    if (isinf(f) || isnan(f))
      //ls->step *= armP->beta_inf;
      _step *= _beta_inf;
    else 
    {
      // Check descent condition
      //if (armP->nondescending && *f <= ref - ls->step*fact*ref)
      if (_nondescending && f <= ref - _step*fact*ref)
      {
        break;
      }
      //if (!armP->nondescending && *f <= ref + ls->step*fact) {
      if (!_nondescending && f <= ref + _step*fact) 
      {
        break;
      }

      //ls->step *= armP->beta;
      _step *= _beta;
    }
  }

  // Check termination
  //if (PetscIsInfOrNanReal(*f)) 
  if (isinf(f) || isnan(f)) 
  {
    //ierr = PetscInfo(ls, "Function is inf or nan.\n"); CHKERRQ(ierr);
    printf( "Function is inf or nan.\n");
    //ls->reason = TAOLINESEARCH_FAILED_INFORNAN;
    _reason = LINESEARCH_FAILED_INFORNAN;
  } 
  //else if (ls->step < ls->stepmin) 
  else if (_step < _stepmin) 
  {
    //ierr = PetscInfo(ls, "Step length is below tolerance.\n"); CHKERRQ(ierr);
    printf("Step length is below tolerance.\n");
    //ls->reason = TAOLINESEARCH_HALTED_RTOL;
    _reason = LINESEARCH_HALTED_RTOL;
  } 
  //else if ((ls->nfeval+ls->nfgeval) >= ls->max_funcs) 
  else if ((_nfeval+_nfgeval) >= _max_funcs) 
  {
    //ierr = PetscInfo2(ls, "Number of line search function evals (%D) > maximum allowed (%D)\n",ls->nfeval+ls->nfgeval, ls->max_funcs); CHKERRQ(ierr);
    printf("Number of line search function evals (%i) > maximum allowed (%i)\n",_nfeval+_nfgeval, _max_funcs);
    //ls->reason = TAOLINESEARCH_HALTED_MAXFCN;
    _reason = LINESEARCH_HALTED_MAXFCN;
  } 
  //if (ls->reason) 
  if (_reason) 
    //PetscFunctionReturn(0);
    return;

  // Successful termination, update memory 
  //armP->lastReference = ref;
  _lastReference = ref;
  //if (armP->replacementPolicy == REPLACE_FIFO) 
  if (_replacementPolicy == REPLACE_FIFO) 
  {
    //armP->memory[armP->current++] = *f;
    _memory[_current++] = f;
    //if (armP->current >= armP->memorySize) 
    if (_current >= _memorySize) 
    {
      //armP->current = 0;
      _current = 0;
    }
  } 
  else 
  {
    //armP->current = idx;
    _current = idx;
    //armP->memory[idx] = *f;
    _memory[idx] = f;
  }

  // Update iterate and compute gradient
  //ierr = VecCopy(armP->work,x); CHKERRQ(ierr);
  x = _work;

  // TK: Not supporting isolated gradient computations right now
  //if (!g_computed) 
  //  ierr = TaoLineSearchComputeGradient(ls, x, g); CHKERRQ(ierr);

  // Finish computations
  //ierr = PetscInfo2(ls, "%D function evals in line search, step = %G\n",ls->nfeval, ls->step); CHKERRQ(ierr);
  printf("%i function evals in line search, step = %G\n",_nfeval, (double)_step);

  //PetscFunctionReturn(0);
}
