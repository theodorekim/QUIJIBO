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

#include "TAO_SOLVER.h"
#include "BLMVM.h"
#include "MORE_THUENTE.h"
#include "UNIT_LINE_SEARCH.h"
#include "ARMIJO.h"
#include "GPCG_LINE_SEARCH.h"
#include <iostream>
#include <cstdio>

using namespace std;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TAO_SOLVER::TAO_SOLVER()
{
  //_max_it     = 10000;
  //_max_funcs   = 10000;
  
  // this is the default for BLMVM
  _max_it      = 2000;
  _max_funcs   = 4000;
  _fatol       = 1e-8;
  _frtol       = 1e-8;
  _gatol       = 1e-8;
  _grtol       = 1e-8;
  _gttol       = 0.0;
  _catol       = 0.0;
  _crtol       = 0.0;
  _xtol        = 0.0;
  _steptol       = 0.0;
  _fmin        = -1e100;
  _hist_max = 0;
  _hist_len = 0;
  _cnorm = 0.0;
  _cnorm0 = 0.0;

  _maxScore = 1;
  _whichLineSearch = USING_MORE_THUENTE;
}

///////////////////////////////////////////////////////////////////////
// solution - contains the initial guess, and will contain the 
//            final solution
// functionGradient - pointer to the function being optimized
///////////////////////////////////////////////////////////////////////
void TAO_SOLVER::optimizeBLMVM(VECTOR& solution, VECTOR& upperBound, VECTOR& lowerBound, 
                               void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient),
                               void (*monitor)())
{
  cout << "======================================================================== " << endl;
  cout << " BEGIN TAO BLMVM" << endl;
  cout << "======================================================================== " << endl;

  // create the BLMVM object
  BLMVM lmvm;

  //TaoSolverTerminationReason reason = TAO_CONTINUE_ITERATING;
  //TaoLineSearchTerminationReason ls_status = TAOLINESEARCH_CONTINUE_ITERATING;

  //PetscReal f, fold, gdx, gnorm;
  //PetscReal stepsize = 1.0,delta;
  //PetscInt iter = 0;
  Real f; 
  Real fold; 
  Real gdx;
  Real gnorm;
  //Real delta;
  Real stepsize = 1.0;
  int iter = 0;
  
  //PetscFunctionBegin;
  
  //  Project initial point onto bounds

  // This is irrelevant, because no bounds function was ever set.
  //TaoComputeVariableBounds(tao);
  VECTOR& XL = lowerBound;
  VECTOR& XU = upperBound;
  assert(XL.size() == solution.size());
  assert(XU.size() == solution.size());
  solution = median(XL, XU, solution);

  //TaoLineSearchSetVariableBounds(tao->linesearch,tao->XL,tao->XU);
  LINE_SEARCH* ls = NULL;
  switch (_whichLineSearch)
  {
    case USING_UNIT:
      ls = new UNIT_LINE_SEARCH(functionGradient);
      break;
    case USING_ARMIJO:
      ls = new ARMIJO(functionGradient);
      break;
    case USING_GPCG:
      ls = new GPCG_LINE_SEARCH(functionGradient);
      break;
    default:
      ls = new MORE_THUENTE(functionGradient);
      break;
  }
  
  ls->_lower = XL;
  ls->_upper = XU;
  ls->_bounded = 1;

  // Check convergence criteria
  //TaoComputeObjectiveAndGradient(tao, tao->solution,&f,blmP->unprojected_gradient);
  VECTOR unprojected_gradient(solution.size());
  VECTOR gradient(solution.size());

  // the generic function-gradient call
  functionGradient(solution, f, unprojected_gradient);
  gradient = unprojected_gradient;  

  // at least for the non-negativity stuff I'm interested in, this doesn't do anything.
  lmvm.VecBoundGradientProjection(unprojected_gradient,solution, XL,XU,gradient);

  //VecNorm(tao->gradient,NORM_2,&gnorm);
  //if (PetscIsInfOrNanReal(f) || PetscIsInfOrNanReal(gnorm)) {
  //  SETERRQ(PETSC_COMM_SELF,1, "User provided compute function generated Inf pr NaN");
  //}
  //gnorm = gradient.norm2();
  gnorm = gradient.norm();
  //cout << " g norm: " << gnorm << endl;

  // not sure if this does anything I need to worry about
  //TaoMonitor(tao, iter, f, gnorm, 0.0, stepsize, &reason);
  //if (reason != TAO_CONTINUE_ITERATING) {
  //  PetscFunctionReturn(0);
  //}
  // it does assign this, which I'll just go ahead and do directly
  _gnorm0 = gnorm;

  // Set initial scaling for the function
  Real& delta = lmvm._delta;
  if (f != 0.0) {
    delta = 2.0*fabs(f) / (gnorm*gnorm);
  }
  else {
    delta = 2.0 / (gnorm*gnorm);
  }
  //delta = (delta > lmvm._delta_min) ? delta : lmvm._delta_min;
  //delta = (delta < lmvm._delta_max) ? delta : lmvm._delta_max;
  lmvm.SetDelta(delta);

  // Set counter for gradient/reset steps
  int blmP_grad = 0;
  int blmP_reset = 0;

  MATRIX blmP_M;
  VECTOR blmP_Xold;
  VECTOR blmP_Gold;
  VECTOR stepdirection;
  SOLVER_STATE reason = CONTINUE_ITERATING;

  // print out preliminary stats
  printf("\niter = %i\t", iter); 
  printf("objective = %g (%g%%)\t", (double)f, (double)(100 * fabs(f / _maxScore)));
  printf("residual: %g\n",(double)gnorm);

  // Have not converged; continue with Newton method
  while (reason == CONTINUE_ITERATING) 
  {
    // Compute direction
    //MatLMVMUpdate(blmP->M, tao->solution, tao->gradient);
    lmvm.Update(blmP_M, solution, gradient);

    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //MatLMVMSolve(blmP->M, blmP->unprojected_gradient, tao->stepdirection);
    lmvm.MatSolve(blmP_M, unprojected_gradient, stepdirection);

    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //VecBoundGradientProjection(tao->stepdirection,tao->solution,tao->XL,tao->XU,tao->gradient);
    lmvm.VecBoundGradientProjection(stepdirection,solution,XL,XU,gradient);

    // Check for success (descent direction)
    //VecDot(blmP->unprojected_gradient, tao->gradient, &gdx);
    //gdx = unprojected_gradient * gradient;
    gdx = unprojected_gradient.dot(gradient);
    if (gdx <= 0) {
      // Step is not descent or solve was not successful Use steepest descent direction (scaled)
      //++blmP->grad;
      blmP_grad++;

      if (f != 0.0) {
        //delta = 2.0*PetscAbsScalar(f) / (gnorm*gnorm);
        delta = 2.0*fabs(f) / (gnorm*gnorm);
      }
      else {
        delta = 2.0 / (gnorm*gnorm);
      }
      //MatLMVMSetDelta(blmP->M,delta);
      lmvm.SetDelta(delta);
      //MatLMVMReset(blmP->M);
      lmvm.MatReset(blmP_M);
      //MatLMVMUpdate(blmP->M, tao->solution, blmP->unprojected_gradient);
      lmvm.Update(blmP_M, solution, unprojected_gradient);
      //MatLMVMSolve(blmP->M,blmP->unprojected_gradient, tao->stepdirection);
      lmvm.MatSolve(blmP_M, unprojected_gradient, stepdirection);
    } 
    //VecScale(tao->stepdirection,-1.0);
    stepdirection *= -1.0;

    // Perform the linesearch
    fold = f;
    //VecCopy(tao->solution, blmP->Xold);
    blmP_Xold = solution;
    //VecCopy(blmP->unprojected_gradient, blmP->Gold);
    blmP_Gold = unprojected_gradient;

    //TaoLineSearchSetInitialStepLength(tao->linesearch,1.0);
    ls->_initstep = 1.0;
    //TaoLineSearchApply(tao->linesearch, tao->solution, &f, blmP->unprojected_gradient, tao->stepdirection, &stepsize, &ls_status);
    ls->LineSearchApply(solution, f, unprojected_gradient, stepdirection);//, &stepsize, &ls_status);

    // don't need to do this -- it just copies ls->_nfeval etc. to the TAO object,
    // which we're not imitating here anyway.
    //TaoAddLineSearchCounts(tao);

    //if (ls_status != TAOLINESEARCH_SUCCESS && ls_status != TAOLINESEARCH_SUCCESS_USER) 
    //if (ls->_reason != LINESEARCH_SUCCESS && ls._reason != LINESEARCH_SUCCESS_USER && ls._reason != LINESEARCH_CONTINUE_ITERATING) 
    if (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS && ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER)
    {
      cout << " TAO: BLMVM, line search FAILED" << endl;
      cout << " reason: " << ls->_reason << endl;
      // Linesearch failed 
      // Reset factors and use scaled (projected) gradient step
      //++blmP->reset;
      blmP_reset++;

      f = fold;
      //VecCopy(blmP->Xold, tao->solution);
      solution = blmP_Xold  ;
      //VecCopy(blmP->Gold, blmP->unprojected_gradient);
      unprojected_gradient = blmP_Gold;

      if (f != 0.0) 
      {
        //delta = 2.0* PetscAbsScalar(f) / (gnorm*gnorm);
        delta = 2.0 * fabs(f) / (gnorm*gnorm);
      }
      else 
      {
        delta = 2.0 / (gnorm*gnorm);
      }
      //MatLMVMSetDelta(blmP->M,delta);
      lmvm.SetDelta(delta);
      //MatLMVMReset(blmP->M);
      lmvm.MatReset(blmP_M);
      //MatLMVMUpdate(blmP->M, tao->solution, blmP->unprojected_gradient);
      lmvm.Update(blmP_M, solution, unprojected_gradient);
      //MatLMVMSolve(blmP->M, blmP->unprojected_gradient, tao->stepdirection);
      lmvm.MatSolve(blmP_M, unprojected_gradient, stepdirection);
      //VecScale(tao->stepdirection, -1.0);
      stepdirection *= -1.0;

      // This may be incorrect; linesearch has values fo stepmax and stepmin
      //  that should be reset.
      //TaoLineSearchSetInitialStepLength(tao->linesearch,1.0);
      ls->_initstep = 1.0;
      //TaoLineSearchApply(tao->linesearch,tao->solution,&f, blmP->unprojected_gradient, tao->stepdirection,  &stepsize, &ls_status);
      ls->LineSearchApply(solution, f, unprojected_gradient, stepdirection);//, &stepsize, &ls_status);
      //TaoAddLineSearchCounts(tao);

      if (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS && ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER) 
      {
        reason = DIVERGED_LS_FAILURE;
        break;
      }
    }
    else
    {
      //cout << " TAO: BLMVM, line search SUCCESS" << endl;
    }

    // Check for termination
    //VecBoundGradientProjection(blmP->unprojected_gradient, tao->solution, tao->XL, tao->XU, tao->gradient);
    lmvm.VecBoundGradientProjection(unprojected_gradient,solution,XL,XU,gradient);
    //VecNorm(tao->gradient, NORM_2, &gnorm);
    //gnorm = gradient.norm2();
    gnorm = gradient.norm();

    if (isinf(f) || isnan(f) || isinf(gnorm) || isnan(gnorm)) 
    {
      //SETERRQ(PETSC_COMM_SELF,1, "User provided compute function generated Not-a-Number");
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << "User provided compute function generated Not-a-Number" << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

    }
    iter++;

    // Monitor call goes here
    _residual = gnorm;
    _niter = iter;
    _fc = f;
    _step = stepsize;
    DefaultConvergenceTest(ls->_nfeval, ls->_nfgeval, reason);

    // print out some progress
    printf("\niter = %i\t",iter); 
    printf("objective = %g (%g%%)\t", (double)f, (double)(100 * fabs(f / _maxScore)));
    printf("residual: %g\n",(double)gnorm);

    if (monitor != NULL)
      monitor();
  }
  delete ls;

  cout << "======================================================================== " << endl;
  cout << " END TAO BLMVM" << endl;
  cout << "======================================================================== " << endl;
  //exit(0);
}

#define LMVM_BFGS                0
#define LMVM_SCALED_GRADIENT     1
#define LMVM_GRADIENT            2
///////////////////////////////////////////////////////////////////////
// solution - contains the initial guess, and will contain the 
//            final solution
// functionGradient - pointer to the function being optimized
///////////////////////////////////////////////////////////////////////
void TAO_SOLVER::optimizeLMVM(VECTOR& solution, void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient))
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << " BEGIN LMVM TAO CLASS TEST" << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  // create the BLMVM object - this has a superset of the functions we need,
  // so there is no separate LMVM class
  BLMVM lmvm;

  //TaoSolverTerminationReason reason = TAO_CONTINUE_ITERATING;
  //TaoLineSearchTerminationReason ls_status = TAOLINESEARCH_CONTINUE_ITERATING;

  //PetscReal f, fold, gdx, gnorm;
  //PetscReal stepsize = 1.0,delta;
  //PetscInt iter = 0;
  Real f; 
  Real fold; 
  Real gdx;
  Real gnorm;
  //Real delta;
  Real stepsize = 1.0;
  int iter = 0;
  
  //PetscFunctionBegin;
  
  LINE_SEARCH* ls = NULL;
  switch (_whichLineSearch)
  {
    case USING_UNIT:
      ls = new UNIT_LINE_SEARCH(functionGradient);
      break;
    case USING_ARMIJO:
      ls = new ARMIJO(functionGradient);
      break;
    case USING_GPCG:
      ls = new GPCG_LINE_SEARCH(functionGradient);
      break;
    default:
      ls = new MORE_THUENTE(functionGradient);
      break;
  }
  ls->_bounded = 0;

  // Check convergence criteria
  //TaoComputeObjectiveAndGradient(tao, tao->solution, &f, tao->gradient);
  VECTOR gradient(solution.size());
  functionGradient(solution, f, gradient);

  //VecNorm(tao->gradient,NORM_2,&gnorm);
  //if (PetscIsInfOrNanReal(f) || PetscIsInfOrNanReal(gnorm)) {
  //  SETERRQ(PETSC_COMM_SELF,1, "User provided compute function generated Inf pr NaN");
  //}
  //gnorm = gradient.norm2();
  gnorm = gradient.norm();
  //cout << " g norm: " << gnorm << endl;

  // not sure if this does anything I need to worry about
  //TaoMonitor(tao, iter, f, gnorm, 0.0, stepsize, &reason);
  //if (reason != TAO_CONTINUE_ITERATING) {
  //  PetscFunctionReturn(0);
  //}
  // it does assign this, which I'll just go ahead and do directly
  _gnorm0 = gnorm;

  // Set initial scaling for the function
  Real& delta = lmvm._delta;
  if (f != 0.0) {
    delta = 2.0*fabs(f) / (gnorm*gnorm);
  }
  else {
    delta = 2.0 / (gnorm*gnorm);
  }
  //delta = (delta > lmvm._delta_min) ? delta : lmvm._delta_min;
  //delta = (delta < lmvm._delta_max) ? delta : lmvm._delta_max;
  lmvm.SetDelta(delta);

  // Set counter for gradient/reset steps
  int lmP_bfgs  = 0;
  int lmP_sgrad = 0;
  int lmP_grad  = 0;

  MATRIX lmP_M;
  VECTOR lmP_Xold;
  VECTOR lmP_Gold;
  VECTOR lmP_D;
  SOLVER_STATE reason = CONTINUE_ITERATING;
  int stepType;
  Real step = 1.0;

  // print out preliminary stats
  printf("\niter = %i\t", iter); 
  printf("objective = %g\t", (double)f);
  printf("residual: %g\n",(double)gnorm);

  // Have not converged; continue with Newton method
  while (reason == CONTINUE_ITERATING) 
  {
    // Compute direction
    //MatLMVMUpdate(blmP->M, tao->solution, tao->gradient);
    lmvm.Update(lmP_M, solution, gradient);

    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //MatLMVMSolve(lmP->M, tao->gradient, lmP->D);
    lmvm.MatSolve(lmP_M, gradient, lmP_D);

    //++lmP->bfgs;
    lmP_bfgs++;

    // Check for success (descent direction)
    // VecDot(lmP->D, tao->gradient, &gdx);
    gdx = lmP_D.dot(gradient);
    if (gdx <= 0 || isinf(gdx) || isnan(gdx)) {
      // Step is not descent or direction produced not a number
      //   We can assert bfgsUpdates > 1 in this case because
      //   the first solve produces the scaled gradient direction,
      //   which is guaranteed to be descent
      //  
      //   Use steepest descent direction (scaled)
      //++lmP->grad;
      lmP_grad++;

      if (f != 0.0) {
        //delta = 2.0*PetscAbsScalar(f) / (gnorm*gnorm);
        delta = 2.0*fabs(f) / (gnorm*gnorm);
      }
      else {
        delta = 2.0 / (gnorm*gnorm);
      }
      //MatLMVMSetDelta(blmP->M,delta);
      lmvm.SetDelta(delta);
      //MatLMVMReset(blmP->M);
      lmvm.MatReset(lmP_M);
      //MatLMVMUpdate(lmP->M, tao->solution, tao->gradient);
      lmvm.Update(lmP_M, solution, gradient);
      //MatLMVMSolve(lmP->M,tao->gradient, lmP->D);
      lmvm.MatSolve(lmP_M, gradient, lmP_D);
      
      //lmP->bfgs = 1;
      lmP_bfgs = 1;
      //++lmP->sgrad;
      lmP_sgrad++;
      //stepType = LMVM_SCALED_GRADIENT;
      stepType = LMVM_SCALED_GRADIENT;

      // TODO: here
    }
    else
    {
      //if (1 == lmP->bfgs) 
      if (1 == lmP_bfgs) 
      {
        //  The first BFGS direction is always the scaled gradient
        //++lmP->sgrad;
        lmP_sgrad++;
        stepType = LMVM_SCALED_GRADIENT;
      }
      else {
        //++lmP->bfgs;
        lmP_bfgs++;
        stepType = LMVM_BFGS;
      }
    }
    //VecScale(lmP->D, -1.0);
    lmP_D *= -1.0;

    // Perform the linesearch
    fold = f;
    //VecCopy(tao->solution, lmP->Xold);
    lmP_Xold = solution;
    //VecCopy(tao->gradient, lmP->Gold);
    lmP_Gold = gradient;

    //TaoLineSearchSetInitialStepLength(tao->linesearch,1.0);
    ls->_initstep = 1.0;
    //TaoLineSearchApply(tao->linesearch, tao->solution, &f, tao->gradient, lmP->D, &step,&ls_status);
    ls->LineSearchApply(solution, f, gradient, lmP_D);//, &stepsize, &ls_status);

    // don't need to do this -- it just copies ls->_nfeval etc. to the TAO object,
    // which we're not imitating here anyway.
    //TaoAddLineSearchCounts(tao);

    //while (ls_status != TAOLINESEARCH_SUCCESS && 
    // ls_status != TAOLINESEARCH_SUCCESS_USER
    // && (stepType != LMVM_GRADIENT)) {
    while (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS && 
           ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER &&
           stepType != LMVM_GRADIENT)
    {
      cout << " TAO: LMVM, line search FAILED" << endl;
      // Linesearch failed 
      //  Reset factors and use scaled gradient step 
      f = fold;
      //VecCopy(lmP->Xold, tao->solution);
      solution = lmP_Xold  ;
      //VecCopy(lmP->Gold, tao->gradient);
      gradient = lmP_Gold;

      switch(stepType)
      {
        case LMVM_BFGS:
          //  Failed to obtain acceptable iterate with BFGS step
          //  Attempt to use the scaled gradient direction

          if (f != 0.0) {
            //delta = 2.0 * PetscAbsScalar(f) / (gnorm*gnorm);
            delta = 2.0 * fabs(f) / (gnorm*gnorm);
          }
          else {
            delta = 2.0 / (gnorm*gnorm);
          }
          //MatLMVMSetDelta(lmP->M, delta);
          lmvm.SetDelta(delta);
          //MatLMVMReset(lmP->M);
          lmvm.MatReset(lmP_M);
          //MatLMVMUpdate(lmP->M, tao->solution, tao->gradient);
          lmvm.Update(lmP_M, solution, gradient);
          //MatLMVMSolve(lmP->M, tao->gradient, lmP->D);
          lmvm.MatSolve(lmP_M, gradient, lmP_D);

          // On a reset, the direction cannot be not a number; it is a 
          // scaled gradient step.  No need to check for this condition.
          //lmP->bfgs = 1;
          lmP_bfgs = 1;
          //++lmP->sgrad;
          lmP_sgrad++;
          stepType = LMVM_SCALED_GRADIENT;
          break;

        case LMVM_SCALED_GRADIENT:
          // The scaled gradient step did not produce a new iterate;
          // attempt to use the gradient direction.
          // Need to make sure we are not using a different diagonal scaling
          //MatLMVMSetDelta(lmP->M, 1.0);
          Real one = 1.0;
          lmvm.SetDelta(one);
          //MatLMVMReset(lmP->M);
          lmvm.MatReset(lmP_M);
          //MatLMVMUpdate(lmP->M, tao->solution, tao->gradient);
          lmvm.Update(lmP_M, solution, gradient);
          //MatLMVMSolve(lmP->M, tao->gradient, lmP->D);
          lmvm.MatSolve(lmP_M, gradient, lmP_D);

          //lmP->bfgs = 1;
          lmP_bfgs = 1;
          //++lmP->grad;
          lmP_grad++;
          stepType = LMVM_GRADIENT;
          break;
        }
      //VecScale(lmP->D, -1.0);
      lmP_D *= (Real)-1.0;
      //TaoLineSearchApply(tao->linesearch, tao->solution, &f, tao->gradient, lmP->D, &step, &ls_status);
      ls->LineSearchApply(solution, f, gradient, lmP_D);//, &stepsize, &ls_status);
    }
    if (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS &&
        ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER)
    {
      //  Failed to find an improving point
      f = fold;
      //VecCopy(lmP->Xold, tao->solution);
      solution = lmP_Xold;
      //VecCopy(lmP->Gold, tao->gradient);
      gradient = lmP_Gold;
      step = 0.0;
      reason = DIVERGED_LS_FAILURE;
    }

    // Check for termination
    gnorm = gradient.norm();

    if (isinf(f) || isnan(f) || isinf(gnorm) || isnan(gnorm)) 
    {
      //SETERRQ(PETSC_COMM_SELF,1, "User provided compute function generated Not-a-Number");
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << "User provided compute function generated Not-a-Number" << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

    }
    iter++;

    // Monitor call goes here
    _residual = gnorm;
    _niter = iter;
    _fc = f;
    _step = stepsize;
    DefaultConvergenceTest(ls->_nfeval, ls->_nfgeval, reason);

    // print out some progress
    printf("\niter = %i\t",iter); 
    //printf("objective = %G \t", (double)f);
    printf("objective = %g\t", (double)f);
    printf("residual: %g\n",(double)gnorm);
  }
  delete ls;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << " END LMVM TAO CLASS TEST" << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //exit(0);
}

////////////////////////////////////////////////////////////////////////////////////////
// return a vector with the median values from all three vectors
////////////////////////////////////////////////////////////////////////////////////////
VECTOR TAO_SOLVER::median(const VECTOR& v0, const VECTOR& v1, const VECTOR& v2)
{
  assert(v0.size() == v1.size());
  assert(v0.size() == v2.size());

  VECTOR final(v2.size());
  vector<Real> sortable(3);
  for (int x = 0; x < v2.size(); x++)
  {
    sortable[0] = v0[x];
    sortable[1] = v1[x];
    sortable[2] = v2[x];
    sort(sortable.begin(), sortable.end());

    final[x] = sortable[1];
  }

  return final;
}

////////////////////////////////////////////////////////////////////////////////////////
// Originally:
//
// PetscErrorCode TaoDefaultConvergenceTest(TaoSolver tao,void *dummy)
////////////////////////////////////////////////////////////////////////////////////////
void TAO_SOLVER::DefaultConvergenceTest(const int functionEvals, const int functionGradEvals, SOLVER_STATE& reason)
{
  int niter = _niter;
  int nfuncs = std::max(functionEvals, functionGradEvals);

  int max_funcs = _max_funcs;
  Real gnorm    = _residual;
  Real gnorm0   = _gnorm0;
  Real f        = _fc;
  Real steptol  = _steptol;
  Real trradius = _step;
  Real gatol    = _gatol;
  Real grtol    = _grtol;
  Real gttol    = _gttol;
  Real fatol    = _fatol;
  Real frtol    = _frtol;
  Real catol    = _catol;
  Real crtol    = _crtol;
  Real fmin     = _fmin;
  Real cnorm    = _cnorm;
  Real cnorm0   = _cnorm0;
  Real gnorm2;
  
  //TaoSolverTerminationReason reason=tao->reason;
  //PetscErrorCode ierr;

  //PetscFunctionBegin;
  //PetscValidHeaderSpecific(tao, TAOSOLVER_CLASSID,1);
  
  if (reason != CONTINUE_ITERATING) 
  {
    //PetscFunctionReturn(0);
    return;
  }
  gnorm2 = gnorm * gnorm;

  if (isinf(f) || isnan(f)) 
  {
    printf("Failed to converged, function value is Inf or NaN\n");
    reason = DIVERGED_NAN;
  } 
  else if (f <= fmin && cnorm <=catol) 
  {
    printf("Converged due to function value %G < minimum function value %G\n", (double)f,(double)fmin);
    reason = CONVERGED_MINF;
  } 
  else if (gnorm2 <= fatol && cnorm <=catol) 
  {
    printf("Converged due to estimated f(X) - f(X*) = %G < %G\n",(double)gnorm2,(double)fatol);
    reason = CONVERGED_FATOL;
  } 
  else if (f != 0 && gnorm2 / fabs(f)<= frtol && cnorm/std::max(cnorm0,(Real)1.0) <= crtol) 
  {
    printf("Converged due to estimated |f(X)-f(X*)|/f(X) = %G < %G\n",(double)(gnorm2/fabs(f)),(double)frtol); 
    reason = CONVERGED_FRTOL;
  } 
  else if (gnorm<= gatol && cnorm <=catol) 
  {
    printf("Converged due to residual norm ||g(X)||=%G < %G\n",(double)gnorm,(double)gatol); 
    reason = CONVERGED_GATOL;
  } 
  else if ( f!=0 && fabs(gnorm/f) <= grtol && cnorm <= crtol) 
  {
    printf("Converged due to residual ||g(X)||/|f(X)| =%G < %G\n",(double)(gnorm/f),(double)grtol); 
    reason = CONVERGED_GRTOL;
  } 
  else if (gnorm0 != 0 && gnorm/gnorm0 <= gttol && cnorm <= crtol) 
  {
    printf("Converged due to relative residual norm ||g(X)||/||g(X0)|| = %G < %G\n",(double)(gnorm/gnorm0),(double)gttol); 
    reason = CONVERGED_GTTOL;
  } 
  else if (nfuncs > max_funcs)
  {
    printf("Exceeded maximum number of function evaluations: %i > %i\n", nfuncs,max_funcs); 
    reason = DIVERGED_MAXFCN;
  } 
  //else if ( tao->lsflag != 0 )
  //else if (ls->_reason != 0 )
  else if (reason == DIVERGED_LS_FAILURE)
  {
    printf("Tao Line Search failure.\n"); 
    reason = DIVERGED_LS_FAILURE;
  } 
  else if (trradius < steptol && niter > 0)
  {
    printf("Trust region/step size too small: %G < %G\n", (double)trradius,(double)steptol); 
    reason = CONVERGED_STEPTOL;
  } 
  //else if (niter > tao->max_it) 
  else if (niter > _max_it) 
  {
    printf("Exceeded maximum number of iterations: %i > %i\n",niter,_max_it); 
    reason = DIVERGED_MAXITS;
  } 
  else 
  {
    reason = CONTINUE_ITERATING;
  }
  //tao->reason = reason;

  //PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////////
// Originally:
/*@
  TaoSetTolerances - Sets parameters used in TAO convergence tests

  Logically collective on TaoSolver

  Input Parameters:
+ tao - the TaoSolver context
. fatol - absolute convergence tolerance
. frtol - relative convergence tolerance
. gatol - stop if norm of gradient is less than this
. grtol - stop if relative norm of gradient is less than this
- gttol - stop if norm of gradient is reduced by this factor

  Options Database Keys:
+ -tao_fatol <fatol> - Sets fatol
. -tao_frtol <frtol> - Sets frtol
. -tao_gatol <gatol> - Sets gatol
. -tao_grtol <grtol> - Sets grtol
- -tao_gttol <gttol> - Sets gttol

  Stopping Criteria:
$ f(X) - f(X*) (estimated)            <= fatol 
$ |f(X) - f(X*)| (estimated) / |f(X)| <= frtol
$ ||g(X)||                            <= gatol
$ ||g(X)|| / |f(X)|                   <= grtol
$ ||g(X)|| / ||g(X0)||                <= gttol

  Notes: 
  Use PETSC_DEFAULT to leave one or more tolerances unchanged.

  Level: beginner

.seealso: TaoGetTolerances()

@*/
////////////////////////////////////////////////////////////////////////////////////////
void TAO_SOLVER::SetTolerances(Real fatol, Real frtol, Real gatol, Real grtol, Real gttol)
{
  if (fatol<0) {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    printf("Tried to set negative fatol -- ignored.\n");
  } else {
    _fatol = std::max((Real)0,fatol);
  }

  if (frtol<0) {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    printf("Tried to set negative frtol -- ignored.\n");
  } else {
    _frtol = std::max((Real)0,frtol);
  }

  if (gatol<0) {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    printf("Tried to set negative gatol -- ignored.\n");
  } else {
    _gatol = std::max((Real)0,gatol);
  }

  if (grtol<0) {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    printf("Tried to set negative grtol -- ignored.\n");
  } else {
    _grtol = std::max((Real)0,grtol);
  }

  if (gttol<0) {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    printf("Tried to set negative gttol -- ignored.\n");
  } else {
    _gttol = std::max((Real)0,gttol);
  }
}

#define NLS_NEWTON    0
#define NLS_BFGS    1
#define NLS_SCALED_GRADIENT   2
#define NLS_GRADIENT    3

#define NLS_KSP_CG	0
#define NLS_KSP_NASH	1
#define NLS_KSP_STCG	2
#define NLS_KSP_GLTR	3
#define NLS_KSP_PETSC	4
#define NLS_KSP_TYPES	5

#define NLS_PC_NONE	0
#define NLS_PC_AHESS	1
#define NLS_PC_BFGS	2
#define NLS_PC_PETSC	3
#define NLS_PC_TYPES	4

#define BFGS_SCALE_AHESS	0
#define BFGS_SCALE_PHESS	1
#define BFGS_SCALE_BFGS		2
#define BFGS_SCALE_TYPES	3

#define NLS_INIT_CONSTANT         0
#define NLS_INIT_DIRECTION        1
#define NLS_INIT_INTERPOLATION    2
#define NLS_INIT_TYPES            3

#define NLS_UPDATE_STEP           0
#define NLS_UPDATE_REDUCTION      1
#define NLS_UPDATE_INTERPOLATION  2
#define NLS_UPDATE_TYPES          3

////////////////////////////////////////////////////////////////////////////////////////
// solution - contains the initial guess, and will contain the final solution
// functionGradient - pointer to the function being optimized
// buildHessian - pointer to the Hessian of the function being optimized
////////////////////////////////////////////////////////////////////////////////////////
void TAO_SOLVER::optimizeNLS(VECTOR& solution, 
                  void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient),
                  void (*buildHessian)(const VECTOR& state, MATRIX& hessian))
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << " BEGIN NLS TAO CLASS TEST" << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  // initializations from TaoCreate_NLS(TaoSolver tao)
  Real nlsP_sval   = 0.0;
  Real nlsP_imin   = 1.0e-4;
  Real nlsP_imax   = 1.0e+2;
  Real nlsP_imfac  = 1.0e-1;

  Real nlsP_pmin   = 1.0e-12;
  Real nlsP_pmax   = 1.0e+2;
  Real nlsP_pgfac  = 1.0e+1;
  Real nlsP_psfac  = 4.0e-1;
  Real nlsP_pmgfac = 1.0e-1;
  Real nlsP_pmsfac = 1.0e-1;

  // build the line search object
  LINE_SEARCH* ls = NULL;
  switch (_whichLineSearch)
  {
    case USING_UNIT:
      ls = new UNIT_LINE_SEARCH(functionGradient);
      break;
    case USING_ARMIJO:
      ls = new ARMIJO(functionGradient);
      break;
    case USING_GPCG:
      ls = new GPCG_LINE_SEARCH(functionGradient);
      break;
    default:
      ls = new MORE_THUENTE(functionGradient);
      break;
  }
  ls->_bounded = 0;

  // TK: Actual NLS code from TAO starts here

  //PetscErrorCode ierr;
  //TAO_NLS *nlsP = (TAO_NLS *)tao->data;
  //PC pc;

  //KSPConvergedReason ksp_reason;
  SOLVER_STATE reason = CONTINUE_ITERATING;
  
  //PetscReal fmin, ftrial, f_full, prered, actred, kappa, sigma;
  Real f_full;
  Real f, fold, gdx, gnorm, pert;
  //PetscReal step = 1.0;
  Real step = 1.0;

  //MatStructure matflag;

  //PetscInt stepType;
  int stepType;
  //PetscInt iter = 0;
  int iter = 0;
  int N;
  //PetscInt needH;
  int needH;
  
  //PetscFunctionBegin;

  // TK: Not possible - can't even pass these in to the function.
  //if (tao->XL || tao->XU || tao->ops->computebounds) {
  //  PetscPrintf(((PetscObject)tao)->comm,"WARNING: Variable bounds have been set but will be ignored by nls algorithm\n");
  //}

  // Initialized variables
  //pert = nlsP->sval;
  pert = nlsP_sval;

  // Check convergence criteria
  //TaoComputeObjectiveAndGradient(tao, tao->solution, &f, tao->gradient);
  VECTOR gradient(solution.size());
  functionGradient(solution, f, gradient);
  //VecNorm(tao->gradient,NORM_2,&gnorm);
  gnorm = gradient.norm();
  //if (PetscIsInfOrNanReal(f) || PetscIsInfOrNanReal(gnorm)) 
  if (isinf(f) || isnan(f) || isinf(gnorm) || isnan(gnorm)) 
  {
    //SETERRQ(PETSC_COMM_SELF,1, "User provided compute function generated Inf or NaN");
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "User provided compute function generated Inf or NaN" << endl;
  }
  needH = 1;

  // Monitor call goes here
  //TaoMonitor(tao, iter, f, gnorm, 0.0, 1.0, &reason);
  _residual = gnorm;
  _niter = iter;
  _fc = f;
  _step = 1.0;
  DefaultConvergenceTest(0,0, reason);
  
  // print out preliminary stats
  printf("\niter = %i\t", iter); 
  printf("objective = %g\t", (double)f);
  printf("residual: %g\n",(double)gnorm);

  if (reason != CONTINUE_ITERATING) 
  {
    //PetscFunctionReturn(0);
    return;
  }

  // Set counter for gradient/reset steps
  //nlsP->newt = 0;
  int nlsP_newt = 0;
  //nlsP->grad = 0;
  int nlsP_grad = 0;

  N = solution.size();
  MATRIX hessian(N,N);
  MATRIX hessian_pre(N,N);

  // search direction found from the matrix solve
  VECTOR nlsP_D;

  // line search vectors
  VECTOR nlsP_Xold;
  VECTOR nlsP_Gold;

  // Have not converged; continue with Newton method
  //while (reason == TAO_CONTINUE_ITERATING) 
  while (reason == CONTINUE_ITERATING) 
  {
    ++iter;

    // Compute the Hessian
    if (needH) 
    {
      //TaoComputeHessian(tao, tao->solution, &tao->hessian, &tao->hessian_pre, &matflag);
      buildHessian(solution, hessian);

      // not doing anything fancy to create a preconditioned version of hessian here
      hessian_pre = hessian;
      needH = 0;
    }

    // Shift the Hessian matrix
    if (pert > 0) 
    {
      //MatShift(tao->hessian, pert);
      hessian = hessian + MATRIX::Identity(N, N) * pert;
      //if (tao->hessian != tao->hessian_pre) 
      //{
      //  MatShift(tao->hessian_pre, pert);
      //}
      hessian_pre = hessian_pre + MATRIX::Identity(N, N) * pert;
    }

    // Solve the Newton system of equations
    ///////////////////////////////////////////////////////////////////////////
    // SOLVER GOES HERE
    // 
    // TK: For now, just do a direct solve. Don't need anything fancy here
    // until the system gets bigger
    ///////////////////////////////////////////////////////////////////////////
    //KSPSetOperators(tao->ksp,tao->hessian,tao->hessian_pre,matflag);
    //KSPSolve(tao->ksp, tao->gradient, nlsP->D);
    //nlsP_D = hessian.colPivHouseholderQr().solve(gradient);
    nlsP_D = hessian * gradient;
    //KSPGetIterationNumber(tao->ksp, &kspits);
    //tao->ksp_its += kspits;

    //VecScale(nlsP->D, -1.0);
    nlsP_D *= (Real)-1.0;

    // Check for success (descent direction)
    //VecDot(nlsP->D, tao->gradient, &gdx);
    gdx = nlsP_D.dot(gradient);

    //if ((gdx >= 0.0) || PetscIsInfOrNanReal(gdx)) 
    if ((gdx >= 0.0) || isinf(gdx) || isnan(gdx)) 
    {
      // Newton step is not descent or direction produced Inf or NaN
      // Update the perturbation for next time
      if (pert <= 0.0) 
      {
        // Initialize the perturbation
        //pert = PetscMin(nlsP->imax, PetscMax(nlsP->imin, nlsP->imfac * gnorm));
        pert = std::min(nlsP_imax, std::max(nlsP_imin, nlsP_imfac * gnorm));
      }
      else {
        // Increase the perturbation
        //pert = PetscMin(nlsP->pmax, PetscMax(nlsP->pgfac * pert, nlsP->pmgfac * gnorm));
        pert = std::min(nlsP_pmax, std::max(nlsP_pgfac * pert, nlsP_pmgfac * gnorm));
      }

      // We don't have the bfgs matrix around and updated
      //   Must use gradient direction in this case
      //VecCopy(tao->gradient, nlsP->D);
      nlsP_D = gradient;
      //VecScale(nlsP->D, -1.0);
      nlsP_D *= -1.0;
      //++nlsP->grad;
      nlsP_grad++;
      stepType = NLS_GRADIENT;
    }
    else {
      // TK: Not using a KSP method here, so the other options are meaningless

      // Newton step computation is good; decrease perturbation
      //pert = PetscMin(nlsP->psfac * pert, nlsP->pmsfac * gnorm);
      pert = std::min(nlsP_psfac * pert, nlsP_pmsfac * gnorm);
      //if (pert < nlsP->pmin) 
      if (pert < nlsP_pmin) 
      {
        pert = 0.0;
      }

      //++nlsP->newt;
      nlsP_newt++;
      stepType = NLS_NEWTON;
    }

    // Perform the linesearch
    fold = f;
    //VecCopy(tao->solution, nlsP->Xold);
    nlsP_Xold = solution;
    //VecCopy(tao->gradient, nlsP->Gold);
    nlsP_Gold = gradient;

    //TaoLineSearchApply(tao->linesearch, tao->solution, &f, tao->gradient, nlsP->D, &step, &ls_reason);
    ls->LineSearchApply(solution, f, gradient, nlsP_D);//, &stepsize, &ls_status);
    //TaoAddLineSearchCounts(tao);

    //while (ls_reason != TAOLINESEARCH_SUCCESS &&
    //       ls_reason != TAOLINESEARCH_SUCCESS_USER &&
    //       stepType != NLS_GRADIENT) 
    while (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS &&
           ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER &&
           stepType != NLS_GRADIENT) 
    {      
      // Linesearch failed
      f = fold;
      //VecCopy(nlsP->Xold, tao->solution);
      solution = nlsP_Xold;
      //VecCopy(nlsP->Gold, tao->gradient);
      gradient = nlsP_Gold;

      // Failed to obtain acceptable iterate with Newton 1step
      // Update the perturbation for next time
      if (pert <= 0.0) {
        // Initialize the perturbation
        //pert = PetscMin(nlsP->imax, PetscMax(nlsP->imin, nlsP->imfac * gnorm));
        pert = std::min(nlsP_imax, std::max(nlsP_imin, nlsP_imfac * gnorm));
      }
      else {
        // Increase the perturbation
        //pert = PetscMin(nlsP->pmax, PetscMax(nlsP->pgfac * pert, nlsP->pmgfac * gnorm));
        pert = std::min(nlsP_pmax, std::max(nlsP_pgfac * pert, nlsP_pmgfac * gnorm));
      }

      // We don't have the bfgs matrix around and being updated
      // Must use gradient direction in this case
      //VecCopy(tao->gradient, nlsP->D);
      nlsP_D = gradient;
      //++nlsP->grad;
      nlsP_grad++;
      stepType = NLS_GRADIENT;

      //VecScale(nlsP->D, -1.0);
      nlsP_D *= (Real)-1.0;

      //TaoLineSearchApply(tao->linesearch, tao->solution, &f, tao->gradient, nlsP->D, &step, &ls_reason);
      ls->LineSearchApply(solution, f, gradient, nlsP_D);
      //TaoLineSearchGetFullStepObjective(tao->linesearch, &f_full);
      f_full = ls->_f_fullstep;
      //TaoAddLineSearchCounts(tao);
    }

    //if (ls_reason != TAOLINESEARCH_SUCCESS &&
    //    ls_reason != TAOLINESEARCH_SUCCESS_USER) 
    if (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS &&
        ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER) 
    {
      // Failed to find an improving point
      f = fold;
      //VecCopy(nlsP->Xold, tao->solution);
      solution = nlsP_Xold;
      //VecCopy(nlsP->Gold, tao->gradient);
      gradient = nlsP_Gold;
      step = 0.0;
      //reason = TAO_DIVERGED_LS_FAILURE;
      reason = DIVERGED_LS_FAILURE;
      //tao->reason = TAO_DIVERGED_LS_FAILURE;
      break;
    }

    //  Check for termination
    //VecNorm(tao->gradient, NORM_2, &gnorm);
    gnorm = gradient.norm();
    //if (PetscIsInfOrNanReal(f) || PetscIsInfOrNanReal(gnorm)) 
    if (isinf(f) || isnan(f) || isinf(gnorm) || isnan(gnorm)) 
    {
      //SETERRQ(PETSC_COMM_SELF,1,"User provided compute function generated Not-a-Number");
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << "User provided compute function generated Not-a-Number" << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    }
    needH = 1;

    //TaoMonitor(tao, iter, f, gnorm, 0.0, step, &reason);
    _residual = gnorm;
    _niter = iter;
    _fc = f;
    _step = step;
    DefaultConvergenceTest(ls->_nfeval, ls->_nfgeval, reason);
    printf("\niter = %i\t",iter); 
    printf("objective = %g\t", (double)f);
    printf("residual: %g\n",(double)gnorm);
  }
  delete ls;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << " END NLS TAO CLASS TEST" << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  //PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////////
// solution - contains the initial guess, and will contain the final solution
// functionGradient - pointer to the function being optimized
// buildHessian - pointer to the Hessian of the function being optimized
////////////////////////////////////////////////////////////////////////////////////////
void TAO_SOLVER::verboseOptimizeNLS(VECTOR& solution, 
                                    void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient),
                                    void (*buildHessian)(const VECTOR& state, MATRIX& hessian))
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << " BEGIN NLS TAO CLASS TEST" << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  // initializations from TaoCreate_NLS(TaoSolver tao)
  Real nlsP_sval   = 0.0;
  Real nlsP_imin   = 1.0e-4;
  Real nlsP_imax   = 1.0e+2;
  Real nlsP_imfac  = 1.0e-1;

  Real nlsP_pmin   = 1.0e-12;
  Real nlsP_pmax   = 1.0e+2;
  Real nlsP_pgfac  = 1.0e+1;
  Real nlsP_psfac  = 4.0e-1;
  Real nlsP_pmgfac = 1.0e-1;
  Real nlsP_pmsfac = 1.0e-1;

  /*
  //  Default values for trust-region radius update based on steplength 
  Real nlsP_nu1 = 0.25;
  Real nlsP_nu2 = 0.50;
  Real nlsP_nu3 = 1.00;
  Real nlsP_nu4 = 1.25;

  Real nlsP_omega1 = 0.25;
  Real nlsP_omega2 = 0.50;
  Real nlsP_omega3 = 1.00;
  Real nlsP_omega4 = 2.00;
  Real nlsP_omega5 = 4.00;

  //  Default values for trust-region radius update based on reduction
  Real nlsP_eta1 = 1.0e-4;
  Real nlsP_eta2 = 0.25;
  Real nlsP_eta3 = 0.50;
  Real nlsP_eta4 = 0.90;

  Real nlsP_alpha1 = 0.25;
  Real nlsP_alpha2 = 0.50;
  Real nlsP_alpha3 = 1.00;
  Real nlsP_alpha4 = 2.00;
  Real nlsP_alpha5 = 4.00;

  //  Default values for trust-region radius update based on interpolation
  Real nlsP_mu1 = 0.10;
  Real nlsP_mu2 = 0.50;

  Real nlsP_gamma1 = 0.25;
  Real nlsP_gamma2 = 0.50;
  Real nlsP_gamma3 = 2.00;
  Real nlsP_gamma4 = 4.00;

  Real nlsP_theta = 0.05;

  //  Default values for trust region initialization based on interpolation 
  Real nlsP_mu1_i = 0.35;
  Real nlsP_mu2_i = 0.50;

  Real nlsP_gamma1_i = 0.0625;
  Real nlsP_gamma2_i = 0.5;
  Real nlsP_gamma3_i = 2.0;
  Real nlsP_gamma4_i = 5.0;
  
  Real nlsP_theta_i = 0.25;

  //  Remaining parameters
  Real nlsP_min_radius = 1.0e-10;
  Real nlsP_max_radius = 1.0e10;
  Real nlsP_epsilon = 1.0e-6;
  */

  //int nlsP_ksp_type        = NLS_KSP_STCG;
  //int nlsP_ksp_type        = NLS_KSP_CG;
  int nlsP_pc_type         = NLS_PC_BFGS;
  //int nlsP_bfgs_scale_type = BFGS_SCALE_PHESS;
  //int nlsP_init_type       = NLS_INIT_INTERPOLATION;
  //int nlsP_update_type     = NLS_UPDATE_STEP;

  // build the line search object
  LINE_SEARCH* ls = NULL;
  switch (_whichLineSearch)
  {
    case USING_UNIT:
      ls = new UNIT_LINE_SEARCH(functionGradient);
      break;
    case USING_ARMIJO:
      ls = new ARMIJO(functionGradient);
      break;
    case USING_GPCG:
      ls = new GPCG_LINE_SEARCH(functionGradient);
      break;
    default:
      ls = new MORE_THUENTE(functionGradient);
      break;
  }
  ls->_bounded = 0;

  // TK: Actual NLS code from TAO starts here

  //PetscErrorCode ierr;
  //TAO_NLS *nlsP = (TAO_NLS *)tao->data;
  //PC pc;

  //KSPConvergedReason ksp_reason;
  SOLVER_STATE reason = CONTINUE_ITERATING;
  
  //PetscReal fmin, ftrial, f_full, prered, actred, kappa, sigma;
  //Real fmin, ftrial, f_full, prered, actred, kappa, sigma;
  Real f_full;
  //PetscReal tau, tau_1, tau_2, tau_max, tau_min, max_radius;
  //Real tau, tau_1, tau_2, tau_max, tau_min, max_radius;
  //PetscReal f, fold, gdx, gnorm, pert;
  Real f, fold, gdx, gnorm, pert;
  //PetscReal step = 1.0;
  Real step = 1.0;

  //PetscReal delta;
  //Real delta;
  //PetscReal norm_d = 0.0, e_min;
  //Real norm_d = 0.0, e_min;

  //MatStructure matflag;

  //PetscInt stepType;
  int stepType = NLS_GRADIENT;
  //PetscInt iter = 0;
  int iter = 0;
  //PetscInt bfgsUpdates = 0;
  //int bfgsUpdates = 0;
  //PetscInt n,N,kspits;
  //int n,N,kspits;
  int N;
  //PetscInt needH;
  int needH;
  
  //PetscInt i_max = 5;
  //int i_max = 5;
  //PetscInt j_max = 1;
  //int j_max = 1;
  //PetscInt i, j;
  //int i, j;

  //PetscFunctionBegin;

  // TK: Not possible - can't even pass these in to the function.
  //if (tao->XL || tao->XU || tao->ops->computebounds) {
  //  PetscPrintf(((PetscObject)tao)->comm,"WARNING: Variable bounds have been set but will be ignored by nls algorithm\n");
  //}

  // Initialized variables
  //pert = nlsP->sval;
  pert = nlsP_sval;

  /*
  int nlsP_ksp_atol = 0;
  int nlsP_ksp_rtol = 0;
  int nlsP_ksp_dtol = 0;
  int nlsP_ksp_ctol = 0;
  int nlsP_ksp_negc = 0;
  int nlsP_ksp_iter = 0;
  int nlsP_ksp_othr = 0;
  */

  // TK: For now, let's just assume we're using CG under the hood,
  // and not trying to do a trust region method.
  // (would use IpOpt for that anyway)
  /*
  // Modify the linear solver to a trust region method if desired
  switch(nlsP->ksp_type) {
    case NLS_KSP_CG:
      KSPSetType(tao->ksp, KSPCG);
      if (tao->ksp->ops->setfromoptions) {
        (*tao->ksp->ops->setfromoptions)(tao->ksp);
      }
      break;

    case NLS_KSP_NASH:
      KSPSetType(tao->ksp, KSPNASH);
      if (tao->ksp->ops->setfromoptions) {
        (*tao->ksp->ops->setfromoptions)(tao->ksp);
      }
      break;

    case NLS_KSP_STCG:
      KSPSetType(tao->ksp, KSPSTCG);
      if (tao->ksp->ops->setfromoptions) {
        (*tao->ksp->ops->setfromoptions)(tao->ksp);
      }
      break;

    case NLS_KSP_GLTR:
      KSPSetType(tao->ksp, KSPGLTR);
      if (tao->ksp->ops->setfromoptions) {
        (*tao->ksp->ops->setfromoptions)(tao->ksp);
      }
      break;

    default:
      // Use the method set by the ksp_type
      break;
  }
  */

  /*
  // Initialize trust-region radius when using nash, stcg, or gltr
  //   Will be reset during the first iteration
  if (NLS_KSP_NASH == nlsP->ksp_type) {
      KSPNASHSetRadius(tao->ksp,nlsP->max_radius);
  } else if (NLS_KSP_STCG == nlsP->ksp_type) {
      KSPSTCGSetRadius(tao->ksp,nlsP->max_radius);
  } else if (NLS_KSP_GLTR == nlsP->ksp_type) {
      KSPGLTRSetRadius(tao->ksp,nlsP->max_radius);
  }
  
  if (NLS_KSP_NASH == nlsP->ksp_type ||
      NLS_KSP_STCG == nlsP->ksp_type || 
      NLS_KSP_GLTR == nlsP->ksp_type) {
    tao->trust = tao->trust0;

    if (tao->trust < 0.0) {
      SETERRQ(PETSC_COMM_SELF,1, "Initial radius negative");
    }

    // Modify the radius if it is too large or small
    tao->trust = PetscMax(tao->trust, nlsP->min_radius);
    tao->trust = PetscMin(tao->trust, nlsP->max_radius);
  }
  */

  /*
  // Get vectors we will need
  if (NLS_PC_BFGS == nlsP->pc_type && !nlsP->M) 
  {
    VecGetLocalSize(tao->solution,&n);
    VecGetSize(tao->solution,&N);
    MatCreateLMVM(((PetscObject)tao)->comm,n,N,&nlsP->M);
    MatLMVMAllocateVectors(nlsP->M,tao->solution);
  }
  */

  // Check convergence criteria
  //TaoComputeObjectiveAndGradient(tao, tao->solution, &f, tao->gradient);
  VECTOR gradient(solution.size());
  functionGradient(solution, f, gradient);
  //VecNorm(tao->gradient,NORM_2,&gnorm);
  gnorm = gradient.norm();
  //if (PetscIsInfOrNanReal(f) || PetscIsInfOrNanReal(gnorm)) 
  if (isinf(f) || isnan(f) || isinf(gnorm) || isnan(gnorm)) 
  {
    //SETERRQ(PETSC_COMM_SELF,1, "User provided compute function generated Inf or NaN");
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "User provided compute function generated Inf or NaN" << endl;
  }
  needH = 1;

  // Monitor call goes here
  //TaoMonitor(tao, iter, f, gnorm, 0.0, 1.0, &reason);
  _residual = gnorm;
  _niter = iter;
  _fc = f;
  _step = 1.0;
  DefaultConvergenceTest(0,0, reason);
  
  // print out preliminary stats
  printf("\niter = %i\t", iter); 
  printf("objective = %g\t", (double)f);
  printf("residual: %g\n",(double)gnorm);

  if (reason != CONTINUE_ITERATING) 
  {
    //PetscFunctionReturn(0);
    return;
  }

  /*
  // create vectors for the limited memory preconditioner
  if ((NLS_PC_BFGS == nlsP->pc_type) && 
      (BFGS_SCALE_BFGS != nlsP->bfgs_scale_type)) {
    if (!nlsP->Diag) {
      VecDuplicate(tao->solution,&nlsP->Diag);
    }
  }

  // Modify the preconditioner to use the bfgs approximation
  KSPGetPC(tao->ksp, &pc);
  switch(nlsP->pc_type) {
    case NLS_PC_NONE:
      PCSetType(pc, PCNONE);
      if (pc->ops->setfromoptions) {
        (*pc->ops->setfromoptions)(pc);
      }
      break;

    case NLS_PC_AHESS:
      PCSetType(pc, PCJACOBI);
      if (pc->ops->setfromoptions) {
        (*pc->ops->setfromoptions)(pc);
      }
      PCJacobiSetUseAbs(pc);
      break;

    case NLS_PC_BFGS:
      PCSetType(pc, PCSHELL);
      if (pc->ops->setfromoptions) {
        (*pc->ops->setfromoptions)(pc);
      }
      PCShellSetName(pc, "bfgs");
      PCShellSetContext(pc, nlsP->M);
      PCShellSetApply(pc, MatLMVMSolveShell);
      break;

    default:
      // Use the pc method set by pc_type
      break;
  }

  // Initialize trust-region radius.  The initialization is only performed 
  //   when we are using Nash, Steihaug-Toint or the Generalized Lanczos method.
  if (NLS_KSP_NASH == nlsP->ksp_type ||
      NLS_KSP_STCG == nlsP->ksp_type || 
      NLS_KSP_GLTR == nlsP->ksp_type) {
    switch(nlsP->init_type) {
      case NLS_INIT_CONSTANT:
        // Use the initial radius specified
        break;

      case NLS_INIT_INTERPOLATION:
        // Use the initial radius specified
        max_radius = 0.0;
    
        for (j = 0; j < j_max; ++j) {
          fmin = f;
          sigma = 0.0;
    
          if (needH) {
            TaoComputeHessian(tao, tao->solution, &tao->hessian, &tao->hessian_pre, &matflag);
            needH = 0;
          }
    
          for (i = 0; i < i_max; ++i) {
            VecCopy(tao->solution,nlsP->W);
            VecAXPY(nlsP->W,-tao->trust/gnorm,tao->gradient);
            TaoComputeObjective(tao, nlsP->W, &ftrial);
            if (PetscIsInfOrNanReal(ftrial)) {
              tau = nlsP->gamma1_i;
            }
            else {
              if (ftrial < fmin) {
                fmin = ftrial;
                sigma = -tao->trust / gnorm;
              }
        
              MatMult(tao->hessian, tao->gradient, nlsP->D);
              VecDot(tao->gradient, nlsP->D, &prered);
    
              prered = tao->trust * (gnorm - 0.5 * tao->trust * prered / (gnorm * gnorm));
              actred = f - ftrial;
              if ((PetscAbsScalar(actred) <= nlsP->epsilon) && 
                  (PetscAbsScalar(prered) <= nlsP->epsilon)) {
                kappa = 1.0;
              }
              else {
                kappa = actred / prered;
              }
    
              tau_1 = nlsP->theta_i * gnorm * tao->trust / (nlsP->theta_i * gnorm * tao->trust + (1.0 - nlsP->theta_i) * prered - actred);
              tau_2 = nlsP->theta_i * gnorm * tao->trust / (nlsP->theta_i * gnorm * tao->trust - (1.0 + nlsP->theta_i) * prered + actred);
              tau_min = PetscMin(tau_1, tau_2);
              tau_max = PetscMax(tau_1, tau_2);
    
              if (PetscAbsScalar(kappa - 1.0) <= nlsP->mu1_i) {
                // Great agreement
                max_radius = PetscMax(max_radius, tao->trust);
    
                if (tau_max < 1.0) {
                  tau = nlsP->gamma3_i;
                }
                else if (tau_max > nlsP->gamma4_i) {
                  tau = nlsP->gamma4_i;
                }
                else if (tau_1 >= 1.0 && tau_1 <= nlsP->gamma4_i && tau_2 < 1.0) {
                  tau = tau_1;
                }
                else if (tau_2 >= 1.0 && tau_2 <= nlsP->gamma4_i && tau_1 < 1.0) {
                  tau = tau_2;
                }
                else {
                  tau = tau_max;
                }
              }
              else if (PetscAbsScalar(kappa - 1.0) <= nlsP->mu2_i) {
                // Good agreement
                max_radius = PetscMax(max_radius, tao->trust);
    
                if (tau_max < nlsP->gamma2_i) {
                  tau = nlsP->gamma2_i;
                }
                else if (tau_max > nlsP->gamma3_i) {
                  tau = nlsP->gamma3_i;
                }
                else {
                  tau = tau_max;
                }
              }
              else {
                // Not good agreement
                if (tau_min > 1.0) {
                  tau = nlsP->gamma2_i;
                }
                else if (tau_max < nlsP->gamma1_i) {
                  tau = nlsP->gamma1_i;
                }
                else if ((tau_min < nlsP->gamma1_i) && (tau_max >= 1.0)) {
                  tau = nlsP->gamma1_i;
                }
                else if ((tau_1 >= nlsP->gamma1_i) && (tau_1 < 1.0) &&
                         ((tau_2 < nlsP->gamma1_i) || (tau_2 >= 1.0))) {
                  tau = tau_1;
                }
                else if ((tau_2 >= nlsP->gamma1_i) && (tau_2 < 1.0) &&
                         ((tau_1 < nlsP->gamma1_i) || (tau_2 >= 1.0))) {
                  tau = tau_2;
                }
                else {
                  tau = tau_max;
                }
              }
            }
            tao->trust = tau * tao->trust;
          }
    
          if (fmin < f) {
            f = fmin;
            VecAXPY(tao->solution,sigma,tao->gradient);
            TaoComputeGradient(tao,tao->solution,tao->gradient);

            VecNorm(tao->gradient,NORM_2,&gnorm);
            if (PetscIsInfOrNanReal(gnorm)) {
              SETERRQ(PETSC_COMM_SELF,1, "User provided compute gradient generated Inf or NaN");
            }
            needH = 1;
    
            TaoMonitor(tao, iter, f, gnorm, 0.0, 1.0, &reason);
            if (reason != TAO_CONTINUE_ITERATING) {
              PetscFunctionReturn(0);
            }
          }
        }
        tao->trust = PetscMax(tao->trust, max_radius);

        // Modify the radius if it is too large or small
        tao->trust = PetscMax(tao->trust, nlsP->min_radius);
        tao->trust = PetscMin(tao->trust, nlsP->max_radius);
        break;

      default:
        // Norm of the first direction will initialize radius
        tao->trust = 0.0;
        break;
    }
  } 

  // Set initial scaling for the BFGS preconditioner
  //   This step is done after computing the initial trust-region radius
  //   since the function value may have decreased
  if (NLS_PC_BFGS == nlsP->pc_type) {
    if (f != 0.0) {
      delta = 2.0 * PetscAbsScalar(f) / (gnorm*gnorm);
    }
    else {
      delta = 2.0 / (gnorm*gnorm);
    }
    MatLMVMSetDelta(nlsP->M,delta);
  }
  */

  // Set counter for gradient/reset steps
  //nlsP->newt = 0;
  int nlsP_newt = 0;
  //nlsP->bfgs = 0;
  //int nlsP_bfgs = 0;
  //nlsP->sgrad = 0;
  //int nlsP_sgrad = 0;
  //nlsP->grad = 0;
  int nlsP_grad = 0;

  N = solution.size();
  MATRIX hessian(N,N);
  MATRIX hessian_pre(N,N);

  // search direction found from the matrix solve
  VECTOR nlsP_D;

  // line search vectors
  VECTOR nlsP_Xold;
  VECTOR nlsP_Gold;

  // Have not converged; continue with Newton method
  //while (reason == TAO_CONTINUE_ITERATING) 
  while (reason == CONTINUE_ITERATING) 
  {
    ++iter;

    // Compute the Hessian
    if (needH) 
    {
      //TaoComputeHessian(tao, tao->solution, &tao->hessian, &tao->hessian_pre, &matflag);
      buildHessian(solution, hessian);
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " Hessian: " << hessian << endl;

      // not doing anything fancy to create a preconditioned version of hessian here
      hessian_pre = hessian;
      needH = 0;
    }

    /*
    if ((NLS_PC_BFGS == nlsP->pc_type) && 
        (BFGS_SCALE_AHESS == nlsP->bfgs_scale_type)) {
      // Obtain diagonal for the bfgs preconditioner 
      MatGetDiagonal(tao->hessian, nlsP->Diag);
      VecAbs(nlsP->Diag);
      VecReciprocal(nlsP->Diag);
      MatLMVMSetScale(nlsP->M,nlsP->Diag);
    }
    */
 
    // Shift the Hessian matrix
    if (pert > 0) 
    {
      //MatShift(tao->hessian, pert);
      hessian = hessian + MATRIX::Identity(N, N) * pert;
      //if (tao->hessian != tao->hessian_pre) 
      //{
      //  MatShift(tao->hessian_pre, pert);
      //}
      hessian_pre = hessian_pre + MATRIX::Identity(N, N) * pert;
    }

    /*
    if (NLS_PC_BFGS == nlsP->pc_type) {
      if (BFGS_SCALE_PHESS == nlsP->bfgs_scale_type) {
        // Obtain diagonal for the bfgs preconditioner
        MatGetDiagonal(tao->hessian, nlsP->Diag);
        VecAbs(nlsP->Diag);
        VecReciprocal(nlsP->Diag);
         MatLMVMSetScale(nlsP->M,nlsP->Diag);
      }
      // Update the limited memory preconditioner
      MatLMVMUpdate(nlsP->M, tao->solution, tao->gradient);
      ++bfgsUpdates;
    }
    */

    // Solve the Newton system of equations
    //KSPSetOperators(tao->ksp,tao->hessian,tao->hessian_pre,matflag);
    /*
    if (NLS_KSP_NASH == nlsP->ksp_type ||
        NLS_KSP_STCG == nlsP->ksp_type || 
        NLS_KSP_GLTR == nlsP->ksp_type) {

      if (NLS_KSP_NASH == nlsP->ksp_type) {
        KSPNASHSetRadius(tao->ksp,nlsP->max_radius);
      } else if (NLS_KSP_STCG == nlsP->ksp_type) {
        KSPSTCGSetRadius(tao->ksp,nlsP->max_radius);
      } else if (NLS_KSP_GLTR == nlsP->ksp_type) {
        KSPGLTRSetRadius(tao->ksp,nlsP->max_radius);
      }
  
      KSPSolve(tao->ksp, tao->gradient, nlsP->D);
      KSPGetIterationNumber(tao->ksp,&kspits);
      tao->ksp_its+=kspits;

      if (NLS_KSP_NASH == nlsP->ksp_type) {
        KSPNASHGetNormD(tao->ksp,&norm_d);
      } else if (NLS_KSP_STCG == nlsP->ksp_type) {
        KSPSTCGGetNormD(tao->ksp,&norm_d);
      } else if (NLS_KSP_GLTR == nlsP->ksp_type) {
        KSPGLTRGetNormD(tao->ksp,&norm_d);
      }

      if (0.0 == tao->trust) {
        // Radius was uninitialized; use the norm of the direction
        if (norm_d > 0.0) {
          tao->trust = norm_d;

          // Modify the radius if it is too large or small
          tao->trust = PetscMax(tao->trust, nlsP->min_radius);
          tao->trust = PetscMin(tao->trust, nlsP->max_radius);
        }
        else {
          // The direction was bad; set radius to default value and re-solve 
          // the trust-region subproblem to get a direction 
          tao->trust = tao->trust0;

          // Modify the radius if it is too large or small
          tao->trust = PetscMax(tao->trust, nlsP->min_radius);
          tao->trust = PetscMin(tao->trust, nlsP->max_radius);

          if (NLS_KSP_NASH == nlsP->ksp_type) {
            KSPNASHSetRadius(tao->ksp,nlsP->max_radius);
          } else if (NLS_KSP_STCG == nlsP->ksp_type) {
            KSPSTCGSetRadius(tao->ksp,nlsP->max_radius);
          } else if (NLS_KSP_GLTR == nlsP->ksp_type) {
            KSPGLTRSetRadius(tao->ksp,nlsP->max_radius);
          }
        
          KSPSolve(tao->ksp, tao->gradient, nlsP->D);
          KSPGetIterationNumber(tao->ksp,&kspits);
          tao->ksp_its+=kspits;
          if (NLS_KSP_NASH == nlsP->ksp_type) {
            KSPNASHGetNormD(tao->ksp,&norm_d);
          } else if (NLS_KSP_STCG == nlsP->ksp_type) {
            KSPSTCGGetNormD(tao->ksp,&norm_d);
          } else if (NLS_KSP_GLTR == nlsP->ksp_type) {
            KSPGLTRGetNormD(tao->ksp,&norm_d);
          }

          if (norm_d == 0.0) {
            SETERRQ(PETSC_COMM_SELF,1, "Initial direction zero");
          }
        }
      }
    }
    else 
    */
    ///////////////////////////////////////////////////////////////////////////
    // SOLVER GOES HERE
    ///////////////////////////////////////////////////////////////////////////
    // TK: For now, just do a direct solve. Don't need anything fancy here
    // until the system gets bigger
    {
      //KSPSolve(tao->ksp, tao->gradient, nlsP->D);
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " D: " << nlsP_D << endl;
      //cout << " gradient: " << gradient << endl;
      //cout << " Hessian before solve: " << hessian << endl;
      nlsP_D = hessian.colPivHouseholderQr().solve(gradient);
      //nlsP_D = hessian * gradient;
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << " D after solve: " << nlsP_D << endl;
      //KSPGetIterationNumber(tao->ksp, &kspits);
      //tao->ksp_its += kspits;
    }
    //VecScale(nlsP->D, -1.0);
    nlsP_D *= (Real)-1.0;

    /*
    KSPGetConvergedReason(tao->ksp, &ksp_reason);
    if ((KSP_DIVERGED_INDEFINITE_PC == ksp_reason) &&
        (NLS_PC_BFGS == nlsP->pc_type) && (bfgsUpdates > 1)) {
      // Preconditioner is numerically indefinite; reset the 
      // approximate if using BFGS preconditioning.

      if (f != 0.0) {
        delta = 2.0 * PetscAbsScalar(f) / (gnorm*gnorm);
      }
      else {
        delta = 2.0 / (gnorm*gnorm);
      }
      MatLMVMSetDelta(nlsP->M,delta);
      MatLMVMReset(nlsP->M);
      MatLMVMUpdate(nlsP->M, tao->solution, tao->gradient);
      bfgsUpdates = 1;
    }

    if (KSP_CONVERGED_ATOL == ksp_reason) {
      ++nlsP->ksp_atol;
    }
    else if (KSP_CONVERGED_RTOL == ksp_reason) {
      ++nlsP->ksp_rtol;
    }
    else if (KSP_CONVERGED_CG_CONSTRAINED == ksp_reason) {
      ++nlsP->ksp_ctol;
    }
    else if (KSP_CONVERGED_CG_NEG_CURVE == ksp_reason) {
      ++nlsP->ksp_negc;
    }
    else if (KSP_DIVERGED_DTOL == ksp_reason) {
      ++nlsP->ksp_dtol;
    }
    else if (KSP_DIVERGED_ITS == ksp_reason) {
      ++nlsP->ksp_iter;
    }
    else {
      ++nlsP->ksp_othr;
    }
    */ 

    // Check for success (descent direction)
    //VecDot(nlsP->D, tao->gradient, &gdx);
    gdx = nlsP_D.dot(gradient);

    //if ((gdx >= 0.0) || PetscIsInfOrNanReal(gdx)) 
    if ((gdx >= 0.0) || isinf(gdx) || isnan(gdx)) 
    {
      // Newton step is not descent or direction produced Inf or NaN
      // Update the perturbation for next time
      if (pert <= 0.0) 
      {
        // Initialize the perturbation
        //pert = PetscMin(nlsP->imax, PetscMax(nlsP->imin, nlsP->imfac * gnorm));
        pert = std::min(nlsP_imax, std::max(nlsP_imin, nlsP_imfac * gnorm));
        /*
        if (NLS_KSP_GLTR == nlsP->ksp_type) {
          KSPGLTRGetMinEig(tao->ksp,&e_min);
          pert = PetscMax(pert, -e_min);
        }
        */
      }
      else {
        // Increase the perturbation
        //pert = PetscMin(nlsP->pmax, PetscMax(nlsP->pgfac * pert, nlsP->pmgfac * gnorm));
        pert = std::min(nlsP_pmax, std::max(nlsP_pgfac * pert, nlsP_pmgfac * gnorm));
      }

      //if (NLS_PC_BFGS != nlsP->pc_type) 
      if (NLS_PC_BFGS != nlsP_pc_type) 
      {
        // We don't have the bfgs matrix around and updated
        //   Must use gradient direction in this case
        //VecCopy(tao->gradient, nlsP->D);
        nlsP_D = gradient;
        //VecScale(nlsP->D, -1.0);
        nlsP_D *= -1.0;
        //++nlsP->grad;
        nlsP_grad++;
        stepType = NLS_GRADIENT;
      }
      /*
      else {
        // Attempt to use the BFGS direction
        MatLMVMSolve(nlsP->M, tao->gradient, nlsP->D);
        VecScale(nlsP->D, -1.0);

        // Check for success (descent direction)
        VecDot(tao->gradient, nlsP->D, &gdx);
        if ((gdx >= 0) || PetscIsInfOrNanReal(gdx)) {
          // BFGS direction is not descent or direction produced not a number
          // We can assert bfgsUpdates > 1 in this case because
          // the first solve produces the scaled gradient direction,
          // which is guaranteed to be descent

          // Use steepest descent direction (scaled)
          if (f != 0.0) {
            delta = 2.0 * PetscAbsScalar(f) / (gnorm*gnorm);
          }
          else {
            delta = 2.0 / (gnorm*gnorm);
          }
          MatLMVMSetDelta(nlsP->M, delta);
          MatLMVMReset(nlsP->M);
          MatLMVMUpdate(nlsP->M, tao->solution, tao->gradient);
          MatLMVMSolve(nlsP->M, tao->gradient, nlsP->D);
          VecScale(nlsP->D, -1.0);
  
          bfgsUpdates = 1;
          ++nlsP->sgrad;
          stepType = NLS_SCALED_GRADIENT;
        }
        else {
          if (1 == bfgsUpdates) {
            // The first BFGS direction is always the scaled gradient
            ++nlsP->sgrad;
            stepType = NLS_SCALED_GRADIENT;
          }
          else {
            ++nlsP->bfgs;
            stepType = NLS_BFGS;
          }
        }
      }
      */
    }
    else {
      /*
      // Computed Newton step is descent
      switch (ksp_reason) {
      case KSP_DIVERGED_NAN:
      case KSP_DIVERGED_BREAKDOWN:
      case KSP_DIVERGED_INDEFINITE_MAT:
      case KSP_DIVERGED_INDEFINITE_PC:
      case KSP_CONVERGED_CG_NEG_CURVE:
        // Matrix or preconditioner is indefinite; increase perturbation
        if (pert <= 0.0) {
          // Initialize the perturbation
          pert = PetscMin(nlsP->imax, PetscMax(nlsP->imin, nlsP->imfac * gnorm));
          if (NLS_KSP_GLTR == nlsP->ksp_type) {
            KSPGLTRGetMinEig(tao->ksp, &e_min);
            pert = PetscMax(pert, -e_min);
          }
        }
        else {
          // Increase the perturbation
          pert = PetscMin(nlsP->pmax, PetscMax(nlsP->pgfac * pert, nlsP->pmgfac * gnorm));
        }
        break;

      default:
      */
        // TK: Not using a KSP method here, so the other options are meaningless

        // Newton step computation is good; decrease perturbation
        //pert = PetscMin(nlsP->psfac * pert, nlsP->pmsfac * gnorm);
        pert = std::min(nlsP_psfac * pert, nlsP_pmsfac * gnorm);
        //if (pert < nlsP->pmin) 
        if (pert < nlsP_pmin) 
        {
          pert = 0.0;
        }
      //  break; 
      //}

      //++nlsP->newt;
      nlsP_newt++;
      stepType = NLS_NEWTON;
    }

    // Perform the linesearch
    fold = f;
    //VecCopy(tao->solution, nlsP->Xold);
    nlsP_Xold = solution;
    //VecCopy(tao->gradient, nlsP->Gold);
    nlsP_Gold = gradient;

    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout << " solution: " << solution << endl;
    //cout << " gradient: " << gradient << endl;
    //cout << " D: " << nlsP_D << endl;

    //TaoLineSearchApply(tao->linesearch, tao->solution, &f, tao->gradient, nlsP->D, &step, &ls_reason);
    ls->LineSearchApply(solution, f, gradient, nlsP_D);//, &stepsize, &ls_status);
    //TaoAddLineSearchCounts(tao);

    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout << " Step type: " << stepType << endl;
    //while (ls_reason != TAOLINESEARCH_SUCCESS &&
    //       ls_reason != TAOLINESEARCH_SUCCESS_USER &&
    //       stepType != NLS_GRADIENT) 
    while (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS &&
           ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER &&
           stepType != NLS_GRADIENT) 
    {      
      // Linesearch failed
      f = fold;
      //VecCopy(nlsP->Xold, tao->solution);
      solution = nlsP_Xold;
      //VecCopy(nlsP->Gold, tao->gradient);
      gradient = nlsP_Gold;

      switch(stepType) 
      {
        case NLS_NEWTON:
          // Failed to obtain acceptable iterate with Newton 1step
          // Update the perturbation for next time
          if (pert <= 0.0) {
            // Initialize the perturbation
            //pert = PetscMin(nlsP->imax, PetscMax(nlsP->imin, nlsP->imfac * gnorm));
            pert = std::min(nlsP_imax, std::max(nlsP_imin, nlsP_imfac * gnorm));
            /*
            if (NLS_KSP_GLTR == nlsP->ksp_type) 
            {
              KSPGLTRGetMinEig(tao->ksp,&e_min);
              pert = PetscMax(pert, -e_min);
            }
            */
          }
          else {
            // Increase the perturbation
            //pert = PetscMin(nlsP->pmax, PetscMax(nlsP->pgfac * pert, nlsP->pmgfac * gnorm));
            pert = std::min(nlsP_pmax, std::max(nlsP_pgfac * pert, nlsP_pmgfac * gnorm));
          }

          //if (NLS_PC_BFGS != nlsP->pc_type) 
          {
            // We don't have the bfgs matrix around and being updated
            // Must use gradient direction in this case
            //VecCopy(tao->gradient, nlsP->D);
            nlsP_D = gradient;
            //++nlsP->grad;
            nlsP_grad++;
            stepType = NLS_GRADIENT;
          }
          /*
          else {
            // Attempt to use the BFGS direction
            MatLMVMSolve(nlsP->M, tao->gradient, nlsP->D);
            // Check for success (descent direction)
            VecDot(tao->solution, nlsP->D, &gdx);
            if ((gdx <= 0) || PetscIsInfOrNanReal(gdx)) 
            {
              // BFGS direction is not descent or direction produced not a number
              // We can assert bfgsUpdates > 1 in this case
              // Use steepest descent direction (scaled)
              if (f != 0.0) {
                delta = 2.0 * PetscAbsScalar(f) / (gnorm*gnorm);
              }
              else {
                delta = 2.0 / (gnorm*gnorm);
              }
              MatLMVMSetDelta(nlsP->M, delta);
              MatLMVMReset(nlsP->M);
              MatLMVMUpdate(nlsP->M, tao->solution, tao->gradient);
              MatLMVMSolve(nlsP->M, tao->gradient, nlsP->D);
    
              bfgsUpdates = 1;
              ++nlsP->sgrad;
              stepType = NLS_SCALED_GRADIENT;
            }
            else {
              if (1 == bfgsUpdates) {
                // The first BFGS direction is always the scaled gradient
                ++nlsP->sgrad;
                stepType = NLS_SCALED_GRADIENT;
              }
              else {
                ++nlsP->bfgs;
                stepType = NLS_BFGS;
              }
            }
          }
          */
          break;

      /*
        case NLS_BFGS:
          // Can only enter if pc_type == NLS_PC_BFGS
          // Failed to obtain acceptable iterate with BFGS step
          // Attempt to use the scaled gradient direction
          if (f != 0.0) {
            delta = 2.0 * PetscAbsScalar(f) / (gnorm*gnorm);
          }
          else {
            delta = 2.0 / (gnorm*gnorm);
          }
          MatLMVMSetDelta(nlsP->M, delta);
          MatLMVMReset(nlsP->M);
          MatLMVMUpdate(nlsP->M, tao->solution, tao->gradient);
          MatLMVMSolve(nlsP->M, tao->gradient, nlsP->D);

          bfgsUpdates = 1;
          ++nlsP->sgrad;
          stepType = NLS_SCALED_GRADIENT;
          break;

        case NLS_SCALED_GRADIENT:
          // Can only enter if pc_type == NLS_PC_BFGS
          // The scaled gradient step did not produce a new iterate;
          // attemp to use the gradient direction.
          // Need to make sure we are not using a different diagonal scaling
    
          MatLMVMSetScale(nlsP->M,0);
          MatLMVMSetDelta(nlsP->M,1.0);
          MatLMVMReset(nlsP->M);
          MatLMVMUpdate(nlsP->M, tao->solution, tao->gradient);
          MatLMVMSolve(nlsP->M, tao->gradient, nlsP->D);
    
          bfgsUpdates = 1;
          ++nlsP->grad;
          stepType = NLS_GRADIENT;
          break;
          */
        default:
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " Switch fell through here. That's not good." << endl;
          exit(0);
          break;
      }
      //VecScale(nlsP->D, -1.0);
      nlsP_D *= (Real)-1.0;

      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " solution: " << solution << endl;
      cout << " gradient: " << gradient << endl;
      cout << " D: " << nlsP_D << endl;
      //TaoLineSearchApply(tao->linesearch, tao->solution, &f, tao->gradient, nlsP->D, &step, &ls_reason);
      ls->LineSearchApply(solution, f, gradient, nlsP_D);
      //TaoLineSearchGetFullStepObjective(tao->linesearch, &f_full);
      f_full = ls->_f_fullstep;
      //TaoAddLineSearchCounts(tao);
    }

    //if (ls_reason != TAOLINESEARCH_SUCCESS &&
    //    ls_reason != TAOLINESEARCH_SUCCESS_USER) 
    if (ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS &&
        ls->_reason != LINE_SEARCH::LINESEARCH_SUCCESS_USER) 
    {
      // Failed to find an improving point
      f = fold;
      //VecCopy(nlsP->Xold, tao->solution);
      solution = nlsP_Xold;
      //VecCopy(nlsP->Gold, tao->gradient);
      gradient = nlsP_Gold;
      step = 0.0;
      //reason = TAO_DIVERGED_LS_FAILURE;
      reason = DIVERGED_LS_FAILURE;
      //tao->reason = TAO_DIVERGED_LS_FAILURE;
      break;
    }

    /*
    // Update trust region radius
    if (NLS_KSP_NASH == nlsP->ksp_type ||
        NLS_KSP_STCG == nlsP->ksp_type || 
        NLS_KSP_GLTR == nlsP->ksp_type) 
    {
      switch(nlsP->update_type) {
        case NLS_UPDATE_STEP:
          if (stepType == NLS_NEWTON) {
            if (step < nlsP->nu1) {
              // Very bad step taken; reduce radius
              tao->trust = nlsP->omega1 * PetscMin(norm_d, tao->trust);
            }
            else if (step < nlsP->nu2) {
              // Reasonably bad step taken; reduce radius
              tao->trust = nlsP->omega2 * PetscMin(norm_d, tao->trust);
            }
            else if (step < nlsP->nu3) {
              //  Reasonable step was taken; leave radius alone
              if (nlsP->omega3 < 1.0) {
                tao->trust = nlsP->omega3 * PetscMin(norm_d, tao->trust);
              }
              else if (nlsP->omega3 > 1.0) {
                tao->trust = PetscMax(nlsP->omega3 * norm_d, tao->trust);  
              }
            }
            else if (step < nlsP->nu4) {
              //  Full step taken; increase the radius
              tao->trust = PetscMax(nlsP->omega4 * norm_d, tao->trust);  
            }
            else {
              //  More than full step taken; increase the radius
              tao->trust = PetscMax(nlsP->omega5 * norm_d, tao->trust);  
            }
          }
          else {
            //  Newton step was not good; reduce the radius
            tao->trust = nlsP->omega1 * PetscMin(norm_d, tao->trust);
          }
          break;

        case NLS_UPDATE_REDUCTION:
          if (stepType == NLS_NEWTON) {
            //  Get predicted reduction 

            if (NLS_KSP_STCG == nlsP->ksp_type) {
              KSPSTCGGetObjFcn(tao->ksp,&prered);
            } else if (NLS_KSP_NASH == nlsP->ksp_type)  {
              KSPNASHGetObjFcn(tao->ksp,&prered);
            } else {
              KSPGLTRGetObjFcn(tao->ksp,&prered);
            }

            if (prered >= 0.0) {
              //  The predicted reduction has the wrong sign.  This cannot
              //  happen in infinite precision arithmetic.  Step should
              //  be rejected!
              tao->trust = nlsP->alpha1 * PetscMin(tao->trust, norm_d);
            }
            else {
              if (PetscIsInfOrNanReal(f_full)) {
                tao->trust = nlsP->alpha1 * PetscMin(tao->trust, norm_d);
              }
              else {
                //  Compute and actual reduction
                actred = fold - f_full;
                prered = -prered;
                if ((PetscAbsScalar(actred) <= nlsP->epsilon) && 
                    (PetscAbsScalar(prered) <= nlsP->epsilon)) {
                  kappa = 1.0;
                }
                else {
                  kappa = actred / prered;
                }
    
                //  Accept of reject the step and update radius
                if (kappa < nlsP->eta1) {
                  //  Very bad step 
                  tao->trust = nlsP->alpha1 * PetscMin(tao->trust, norm_d);
                }
                else if (kappa < nlsP->eta2) {
                  //  Marginal bad step
                  tao->trust = nlsP->alpha2 * PetscMin(tao->trust, norm_d);
                }
                else if (kappa < nlsP->eta3) {
                  //  Reasonable step
                  if (nlsP->alpha3 < 1.0) {
                    tao->trust = nlsP->alpha3 * PetscMin(norm_d, tao->trust);
                  }
                  else if (nlsP->alpha3 > 1.0) {
                    tao->trust = PetscMax(nlsP->alpha3 * norm_d, tao->trust);  
                  }
                }
                else if (kappa < nlsP->eta4) {
                  //  Good step
                  tao->trust = PetscMax(nlsP->alpha4 * norm_d, tao->trust);
                }
                else {
                  //  Very good step
                  tao->trust = PetscMax(nlsP->alpha5 * norm_d, tao->trust);
                }
              }
            }
          }
          else {
            //  Newton step was not good; reduce the radius
            tao->trust = nlsP->alpha1 * PetscMin(norm_d, tao->trust);
          }
          break;

        default:
          if (stepType == NLS_NEWTON) {

            if (NLS_KSP_STCG == nlsP->ksp_type) {
                KSPSTCGGetObjFcn(tao->ksp,&prered);
            } else if (NLS_KSP_NASH == nlsP->ksp_type)  {
                KSPNASHGetObjFcn(tao->ksp,&prered);
            } else {
                KSPGLTRGetObjFcn(tao->ksp,&prered);
            }
            if (prered >= 0.0) {
              //  The predicted reduction has the wrong sign.  This cannot
              //  happen in infinite precision arithmetic.  Step should
              //  be rejected!
              tao->trust = nlsP->gamma1 * PetscMin(tao->trust, norm_d);
            }
            else {
              if (PetscIsInfOrNanReal(f_full)) {
                tao->trust = nlsP->gamma1 * PetscMin(tao->trust, norm_d);
              }
              else {
                actred = fold - f_full;
                prered = -prered;
                if ((PetscAbsScalar(actred) <= nlsP->epsilon) && 
                    (PetscAbsScalar(prered) <= nlsP->epsilon)) {
                  kappa = 1.0;
                }
                else {
                  kappa = actred / prered;
                }

                tau_1 = nlsP->theta * gdx / (nlsP->theta * gdx - (1.0 - nlsP->theta) * prered + actred);
                tau_2 = nlsP->theta * gdx / (nlsP->theta * gdx + (1.0 + nlsP->theta) * prered - actred);
                tau_min = PetscMin(tau_1, tau_2);
                tau_max = PetscMax(tau_1, tau_2);

                if (kappa >= 1.0 - nlsP->mu1) {
                  //  Great agreement 
                  if (tau_max < 1.0) {
                    tao->trust = PetscMax(tao->trust, nlsP->gamma3 * norm_d);
                  }
                  else if (tau_max > nlsP->gamma4) {
                    tao->trust = PetscMax(tao->trust, nlsP->gamma4 * norm_d);
                  }
                  else {
                    tao->trust = PetscMax(tao->trust, tau_max * norm_d);
                  }
                }
                else if (kappa >= 1.0 - nlsP->mu2) {
                  //  Good agreement
                  if (tau_max < nlsP->gamma2) {
                    tao->trust = nlsP->gamma2 * PetscMin(tao->trust, norm_d);
                  }
                  else if (tau_max > nlsP->gamma3) {
                    tao->trust = PetscMax(tao->trust, nlsP->gamma3 * norm_d);
                  }
                  else if (tau_max < 1.0) {
                    tao->trust = tau_max * PetscMin(tao->trust, norm_d);
                  }
                  else {
                    tao->trust = PetscMax(tao->trust, tau_max * norm_d);
                  }
                }
                else {
                  //  Not good agreement
                  if (tau_min > 1.0) {
                    tao->trust = nlsP->gamma2 * PetscMin(tao->trust, norm_d);
                  }
                  else if (tau_max < nlsP->gamma1) {
                    tao->trust = nlsP->gamma1 * PetscMin(tao->trust, norm_d);
                  }
                  else if ((tau_min < nlsP->gamma1) && (tau_max >= 1.0)) {
                    tao->trust = nlsP->gamma1 * PetscMin(tao->trust, norm_d);
                  }
                  else if ((tau_1 >= nlsP->gamma1) && (tau_1 < 1.0) &&
                           ((tau_2 < nlsP->gamma1) || (tau_2 >= 1.0))) {
                    tao->trust = tau_1 * PetscMin(tao->trust, norm_d);
                  }
                  else if ((tau_2 >= nlsP->gamma1) && (tau_2 < 1.0) &&
                           ((tau_1 < nlsP->gamma1) || (tau_2 >= 1.0))) {
                    tao->trust = tau_2 * PetscMin(tao->trust, norm_d);
                  }
                  else {
                    tao->trust = tau_max * PetscMin(tao->trust, norm_d);
                  }
                }
              } 
            }
          }
          else {
            //  Newton step was not good; reduce the radius
            tao->trust = nlsP->gamma1 * PetscMin(norm_d, tao->trust);
          }
          break;
      }

      //  The radius may have been increased; modify if it is too large
      tao->trust = PetscMin(tao->trust, nlsP->max_radius);
    }
    */

    //  Check for termination
    //VecNorm(tao->gradient, NORM_2, &gnorm);
    gnorm = gradient.norm();
    //if (PetscIsInfOrNanReal(f) || PetscIsInfOrNanReal(gnorm)) 
    if (isinf(f) || isnan(f) || isinf(gnorm) || isnan(gnorm)) 
    {
      //SETERRQ(PETSC_COMM_SELF,1,"User provided compute function generated Not-a-Number");
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << "User provided compute function generated Not-a-Number" << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    }
    needH = 1;

    //TaoMonitor(tao, iter, f, gnorm, 0.0, step, &reason);
    _residual = gnorm;
    _niter = iter;
    _fc = f;
    _step = step;
    DefaultConvergenceTest(ls->_nfeval, ls->_nfgeval, reason);
    printf("\niter = %i\t",iter); 
    printf("objective = %g\t", (double)f);
    printf("residual: %g\n",(double)gnorm);
  }
  delete ls;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << " END NLS TAO CLASS TEST" << endl;
  cout << "======================================================================== " << endl;
  cout << "======================================================================== " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  //PetscFunctionReturn(0);
}
