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

#include "MORE_THUENTE.h"
#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MORE_THUENTE::MORE_THUENTE(void (*functionGradient)(const VECTOR&, Real&, VECTOR&)) :
  LINE_SEARCH(functionGradient)
{
  _bracket = 0;
  _infoc = 1;
  _initstep = 1;
}

////////////////////////////////////////////////////////////////////////////////////////
/*
     The subroutine mcstep is taken from the work of Jorge Nocedal.
     this is a variant of More' and Thuente's routine.

     subroutine mcstep

     the purpose of mcstep is to compute a safeguarded step for
     a linesearch and to update an interval of uncertainty for
     a minimizer of the function.

     the parameter stx contains the step with the least function
     value. the parameter stp contains the current step. it is
     assumed that the derivative at stx is negative in the
     direction of the step. if bracket is set true then a
     minimizer has been bracketed in an interval of uncertainty
     with endpoints stx and sty.

     the subroutine statement is

     subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,bracket,
                       stpmin,stpmax,info)

     where

       stx, fx, and dx are variables which specify the step,
         the function, and the derivative at the best step obtained
         so far. The derivative must be negative in the direction
         of the step, that is, dx and stp-stx must have opposite
         signs. On output these parameters are updated appropriately.

       sty, fy, and dy are variables which specify the step,
         the function, and the derivative at the other endpoint of
         the interval of uncertainty. On output these parameters are
         updated appropriately.

       stp, fp, and dp are variables which specify the step,
         the function, and the derivative at the current step.
         If bracket is set true then on input stp must be
         between stx and sty. On output stp is set to the new step.

       bracket is a logical variable which specifies if a minimizer
         has been bracketed.  If the minimizer has not been bracketed
         then on input bracket must be set false.  If the minimizer
         is bracketed then on output bracket is set true.

       stpmin and stpmax are input variables which specify lower
         and upper bounds for the step.

       info is an integer output variable set as follows:
         if info = 1,2,3,4,5, then the step has been computed
         according to one of the five cases below. otherwise
         info = 0, and this indicates improper input parameters.

     subprograms called

       fortran-supplied ... abs,max,min,sqrt

     argonne national laboratory. minpack project. june 1983
     jorge j. more', david j. thuente

*/
////////////////////////////////////////////////////////////////////////////////////////
void MORE_THUENTE::mcstep(Real& stx, Real& fx, Real& dx,
                          Real& sty, Real& fy, Real& dy,
                          Real& stp, Real& fp, Real& dp)
{
  //TAOLINESEARCH_MT_CTX *mtP = (TAOLINESEARCH_MT_CTX *) ls->data;
  //PetscReal gamma1, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  Real gamma1, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  //PetscInt bound;
  int bound;

  //PetscFunctionBegin;

  // Check the input parameters for errors
  _infoc = 0;
  //if (_bracket && (*stp <= PetscMin(*stx,*sty) || (*stp >= PetscMax(*stx,*sty)))) 
  if (_bracket && (stp <= std::min(stx,sty) || (stp >= std::max(stx,sty))))
  { 
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " bad stp in bracket " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }
  //if (*dx * (*stp-*stx) >= 0.0) 
  if (dx * (stp-stx) >= 0.0)
  { 
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "dx * (stp-stx) >= 0.0" << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }
  //if (ls->stepmax < ls->stepmin) 
  if (_stepmax < _stepmin)
  { 
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "stepmax > stepmin" << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  // Determine if the derivatives have opposite sign
  //sgnd = *dp * (*dx / PetscAbsReal(*dx));
  sgnd = dp * (dx / fabs(dx));

  //if (*fp > *fx) 
  if (fp > fx) 
  {
    // Case 1: a higher function value.
    // The minimum is bracketed. If the cubic step is closer
    // to stx than the quadratic step, the cubic step is taken,
    // else the average of the cubic and quadratic steps is taken.

    //mtP->infoc = 1;
    _infoc = 1;
    bound = 1;
    //theta = 3 * (*fx - *fp) / (*stp - *stx) + *dx + *dp;
    theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    //s = PetscMax(PetscAbsReal(theta),PetscAbsReal(*dx));
    s = std::max(fabs(theta),fabs(dx));
    //s = PetscMax(s,PetscAbsReal(*dp));
    s = std::max(s,fabs(dp));
    //gamma1 = s*PetscSqrtScalar(PetscPowScalar(theta/s,2.0) - (*dx/s)*(*dp/s));
    gamma1 = s*sqrt(pow(theta/s,(Real)2.0) - (dx/s)*(dp/s));
    //if (*stp < *stx) gamma1 = -gamma1;
    if (stp < stx) gamma1 = -gamma1;
    // Can p be 0?  Check
    //p = (gamma1 - *dx) + theta;
    p = (gamma1 - dx) + theta;
    //q = ((gamma1 - *dx) + gamma1) + *dp;
    q = ((gamma1 - dx) + gamma1) + dp;
    r = p/q;
    //stpc = *stx + r*(*stp - *stx);
    stpc = stx + r*(stp - stx);
    //stpq = *stx + ((*dx/((*fx-*fp)/(*stp-*stx)+*dx))*0.5) * (*stp - *stx);
    stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))*0.5) * (stp - stx);

    //if (PetscAbsReal(stpc-*stx) < PetscAbsReal(stpq-*stx)) 
    if (fabs(stpc-stx) < fabs(stpq-stx)) 
    {
      stpf = stpc;
    } 
    else {
      stpf = stpc + 0.5*(stpq - stpc);
    }
    //mtP->bracket = 1;
    _bracket = 1;
  }
  else if (sgnd < 0.0) 
  {
    // Case 2: A lower function value and derivatives of
    // opposite sign. The minimum is bracketed. If the cubic
    // step is closer to stx than the quadratic (secant) step,
    // the cubic step is taken, else the quadratic step is taken.

    //mtP->infoc = 2;
    _infoc = 2;
    bound = 0;
    //theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    theta = 3*(fx - fp)/(stp - stx) + dx + dp;
    //s = PetscMax(PetscAbsReal(theta),PetscAbsReal(*dx));
    s = std::max(fabs(theta),fabs(dx));
    //s = PetscMax(s,PetscAbsReal(*dp));
    s = std::max(s,fabs(dp));
    //gamma1 = s*PetscSqrtScalar(PetscPowScalar(theta/s,2.0) - (*dx/s)*(*dp/s));
    gamma1 = s*sqrt(pow(theta/s,(Real)2.0) - (dx/s)*(dp/s));
    //if (*stp > *stx) gamma1 = -gamma1;
    if (stp > stx) gamma1 = -gamma1;
    //p = (gamma1 - *dp) + theta;
    p = (gamma1 - dp) + theta;
    //q = ((gamma1 - *dp) + gamma1) + *dx;
    q = ((gamma1 - dp) + gamma1) + dx;
    r = p/q;
    //stpc = *stp + r*(*stx - *stp);
    stpc = stp + r*(stx - stp);
    //stpq = *stp + (*dp/(*dp-*dx))*(*stx - *stp);
    stpq = stp + (dp/(dp-dx))*(stx - stp);

    //if (PetscAbsReal(stpc-*stp) > PetscAbsReal(stpq-*stp)) 
    if (fabs(stpc-stp) > fabs(stpq-stp)) 
    {
      stpf = stpc;
    }
    else {
      stpf = stpq;
    }
    //mtP->bracket = 1;
    _bracket = 1;
  }
  //else if (PetscAbsReal(*dp) < PetscAbsReal(*dx)) 
  else if (fabs(dp) < fabs(dx)) 
  {
    // Case 3: A lower function value, derivatives of the
    // same sign, and the magnitude of the derivative decreases.
    // The cubic step is only used if the cubic tends to infinity
    // in the direction of the step or if the minimum of the cubic
    // is beyond stp. Otherwise the cubic step is defined to be
    // either stepmin or stepmax. The quadratic (secant) step is also
    // computed and if the minimum is bracketed then the the step
    // closest to stx is taken, else the step farthest away is taken.

    //mtP->infoc = 3;
    _infoc = 3;
    bound = 1;
    //theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    theta = 3*(fx - fp)/(stp - stx) + dx + dp;
    //s = PetscMax(PetscAbsReal(theta),PetscAbsReal(*dx));
    s = std::max(fabs(theta),fabs(dx));
    //s = PetscMax(s,PetscAbsReal(*dp));
    s = std::max(s,fabs(dp));

    // The case gamma1 = 0 only arises if the cubic does not tend
    //   to infinity in the direction of the step.
    //gamma1 = s*PetscSqrtScalar(PetscMax(0.0,PetscPowScalar(theta/s,2.0) - (*dx/s)*(*dp/s)));
    gamma1 = s*sqrt(std::max((Real)0.0,pow(theta/s,(Real)2.0) - (dx/s)*(dp/s)));
    //if (*stp > *stx) gamma1 = -gamma1;
    if (stp > stx) gamma1 = -gamma1;
    //p = (gamma1 - *dp) + theta;
    p = (gamma1 - dp) + theta;
    //q = (gamma1 + (*dx - *dp)) + gamma1;
    q = (gamma1 + (dx - dp)) + gamma1;
    r = p/q;
    //if (r < 0.0 && gamma1 != 0.0) stpc = *stp + r*(*stx - *stp);
    if (r < 0.0 && gamma1 != 0.0) stpc = stp + r*(stx - stp);
    //else if (*stp > *stx)        stpc = ls->stepmax;
    else if (stp > stx)          stpc = _stepmax;
    else                         stpc = _stepmin;
    //stpq = *stp + (*dp/(*dp-*dx)) * (*stx - *stp);
    stpq = stp + (dp/(dp-dx)) * (stx - stp);

    //if (mtP->bracket) 
    if (_bracket) 
    {
      //if (PetscAbsReal(*stp-stpc) < PetscAbsReal(*stp-stpq)) 
      if (fabs(stp-stpc) < fabs(stp-stpq)) 
      {
	      stpf = stpc;
      } 
      else 
      {
      	stpf = stpq;
      }
    }
    else {
      //if (PetscAbsReal(*stp-stpc) > PetscAbsReal(*stp-stpq)) 
      if (fabs(stp-stpc) > fabs(stp-stpq)) 
      {
      	stpf = stpc;
      }
      else 
      {
      	stpf = stpq;
      }
    }
  }
  else 
  {
    // Case 4: A lower function value, derivatives of the
    //   same sign, and the magnitude of the derivative does
    //   not decrease. If the minimum is not bracketed, the step
    //   is either stpmin or stpmax, else the cubic step is taken.

    //mtP->infoc = 4;
    _infoc = 4;
    bound = 0;
    //if (mtP->bracket) 
    if (_bracket) 
    {
      //theta = 3*(*fp - *fy)/(*sty - *stp) + *dy + *dp;
      theta = 3*(fp - fy)/(sty - stp) + dy + dp;
      //s = PetscMax(PetscAbsReal(theta),PetscAbsReal(*dy));
      s = std::max(fabs(theta),fabs(dy));
      //s = PetscMax(s,PetscAbsReal(*dp));
      s = std::max(s,fabs(dp));
      //gamma1 = s*PetscSqrtScalar(PetscPowScalar(theta/s,2.0) - (*dy/s)*(*dp/s));
      gamma1 = s*sqrt(pow(theta/s,(Real)2.0) - (dy/s)*(dp/s));
      //if (*stp > *sty) gamma1 = -gamma1;
      if (stp > sty) gamma1 = -gamma1;
      //p = (gamma1 - *dp) + theta;
      p = (gamma1 - dp) + theta;
      //q = ((gamma1 - *dp) + gamma1) + *dy;
      q = ((gamma1 - dp) + gamma1) + dy;
      r = p/q;
      //stpc = *stp + r*(*sty - *stp);
      stpc = stp + r*(sty - stp);
      //stpq = *stp + (*dp/(*dp-*dx)) * (*stx - *stp);
      stpq = stp + (dp/(dp-dx)) * (stx - stp);

      stpf = stpc;
    } 
    //else if (*stp > *stx) 
    else if (stp > stx) 
    {
      //stpf = ls->stepmax;
      stpf = _stepmax;
    } 
    else 
    {
      //stpf = ls->stepmin;
      stpf = _stepmin;
    }
  }
  
  // Update the interval of uncertainty.  This update does not
  //   depend on the new step or the case analysis above.
  //if (*fp > *fx) 
  if (fp > fx) 
  {
    //*sty = *stp;
    sty = stp;
    //*fy = *fp;
    fy = fp;
    //*dy = *dp;
    dy = dp;
  } 
  else 
  {
    if (sgnd < 0.0) 
    {
      //*sty = *stx;
      sty = stx;
      //*fy = *fx;
      fy = fx;
      //*dy = *dx;
      dy = dx;
    }
    //*stx = *stp;
    stx = stp;
    //*fx = *fp;
    fx = fp;
    //*dx = *dp;
    dx = dp;
  }
  
  // Compute the new step and safeguard it.
  //stpf = PetscMin(ls->stepmax,stpf);
  stpf = std::min(_stepmax,stpf);
  //stpf = PetscMax(ls->stepmin,stpf);
  stpf = std::max(_stepmin,stpf);
  //*stp = stpf;
  stp = stpf;
  //if (mtP->bracket && bound) 
  if (_bracket && bound) 
  {
    //if (*sty > *stx) 
    if (sty > stx) 
    {
      //*stp = PetscMin(*stx+0.66*(*sty-*stx),*stp);
      stp = std::min(stx+(Real)0.66*(sty-stx),stp);
    }
    else 
    {
      //*stp = PetscMax(*stx+0.66*(*sty-*stx),*stp);
      stp = std::max(stx+(Real)0.66*(sty-stx),stp);
    }
  }
  //PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////////
// The More-Thuente line search
////////////////////////////////////////////////////////////////////////////////////////
void MORE_THUENTE::LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s)
{
  const Real INF = 1.0e100;

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " s: " << s << endl;
  //TAOLINESEARCH_MT_CTX *mt;
  
  //PetscReal    xtrapf = 4.0;
  Real    xtrapf = 4.0;
  
  //PetscReal   finit, width, width1, dginit, fm, fxm, fym, dgm, dgxm, dgym;
  Real   finit, width, width1, dginit, fm, fxm, fym, dgm, dgxm, dgym;
  
  //PetscReal    dgx, dgy, dg, dg2, fx, fy, stx, sty, dgtest;
  Real    dgx, dgy, dg, dg2, fx, fy, stx, sty, dgtest;
  
  //PetscReal ftest1=0.0, ftest2=0.0;
  Real ftest1=0.0, ftest2=0.0;
  
  //PetscInt  i, stage1,n1,n2,nn1,nn2;
  int i, stage1;

  // TK: not used?
  //int n1,n2,nn1,nn2;
  
  //PetscReal bstepmin1, bstepmin2, bstepmax;
  Real bstepmin1 = 0;
  Real bstepmin2 = 0;
  Real bstepmax = 0;
  
  //PetscBool g_computed=PETSC_FALSE;
  bool g_computed=false;

  //PetscFunctionBegin;
  //PetscValidHeaderSpecific(ls,TAOLINESEARCH_CLASSID,1);
  //PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  //PetscValidScalarPointer(f,3);
  //PetscValidHeaderSpecific(g,VEC_CLASSID,4);
  //PetscValidHeaderSpecific(s,VEC_CLASSID,5);

  // comm,type,size checks are done in interface TaoLineSearchApply
  //mt = (TAOLINESEARCH_MT_CTX*)(ls->data);

  _reason = LINESEARCH_CONTINUE_ITERATING;

  // Check work vector
  //if (!mt->work) 
  if (_work.size() == 0)
  {
    //VecDuplicate(x,&mt->work);
    _work = x;
    //mt->x = x;
    _x = x;
    //PetscObjectReference((PetscObject)mt->x);
  }
  // If x has changed, then recreate work
  //else if (x != mt->x) 
  //else if ((_x - x).norm2() > 1e-6)
  else if ((_x - x).norm() > 1e-6)
  { 
    //VecDestroy(&mt->work);
    //VecDuplicate(x,&mt->work);
    _work = x;
    //PetscObjectDereference((PetscObject)mt->x);
    //mt->x = x;
    _x = x;
    //PetscObjectReference((PetscObject)mt->x);
  }
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " s: " << s << endl;
  
  if (_bounded) 
  {
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout << " ENTERED BOUNDED BLOCK " << endl;
    // Compute step length needed to make all variables equal a bound
    // Compute the smallest steplength that will make one nonbinding variable
    //  equal the bound
    //VecGetLocalSize(ls->upper,&n1);
    //VecGetLocalSize(mt->x, &n2);
    //VecGetSize(ls->upper,&nn1);
    //VecGetSize(mt->x,&nn2);
    //if (n1 != n2 || nn1 != nn2) {
    //    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Variable vector not compatible with bounds vector");
    //}

    //VecScale(s,-1.0);
    s *= -1.0;

    //VecBoundGradientProjection(s,x,ls->lower,ls->upper,s);
    //LMVM_VecBoundGradientProjection(s,x,_lower,_upper,s);
    VecBoundGradientProjection(s,x,_lower,_upper,s);
    //cout << " after projection s: " << s << endl;
    
    //VecScale(s,-1.0);
    s *= -1.0;

    //VecStepBoundInfo(x,ls->lower,ls->upper,s,&bstepmin1,&bstepmin2,&bstepmax);
    VecStepBoundInfo(x,_lower,_upper,s,bstepmin1,bstepmin2,bstepmax);
    //ls->stepmax = PetscMin(bstepmax,1.0e15);
    _stepmax = std::min(bstepmax,(Real)1.0e15);
  }

  //VecDot(g,s,&dginit);
  //dginit = g * s;
  dginit = g.dot(s);
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " dginit: " << dginit << endl;
  //cout << " g: " << g << endl;
  //cout << " s: " << s << endl;
  
  //if (PetscIsInfOrNanReal(dginit)) 
  if (isnan(dginit) || isinf(dginit)) 
  {
    //PetscInfo1(ls,"Initial Line Search step * g is Inf or Nan (%G)\n",dginit);
    _reason = LINESEARCH_FAILED_INFORNAN;
    //PetscFunctionReturn(0);
    return;
  }
  if (dginit >= 0.0) 
  {
    //PetscInfo1(ls,"Initial Line Search step * g is not descent direction (%G)\n",dginit);
    //ls->reason = TAOLINESEARCH_FAILED_ASCENT;
    _reason = LINESEARCH_FAILED_ASCENT;
    //PetscFunctionReturn(0);
    return;
  }

  // Initialization
  //mt->bracket = 0;
  _bracket = 0;
  stage1 = 1;
  //finit = *f;
  finit = f;
  //dgtest = ls->ftol * dginit;
  dgtest = _ftol * dginit;
  //width = ls->stepmax - ls->stepmin;
  width = _stepmax - _stepmin;
  width1 = width * 2.0;
  //VecCopy(x,mt->work);
  _work = x;
  // Variable dictionary:  
  // stx, fx, dgx - the step, function, and derivative at the best step
  // sty, fy, dgy - the step, function, and derivative at the other endpoint 
  //               of the interval of uncertainty
  // step, f, dg - the step, function, and derivative at the current step
  stx = 0.0;
  fx  = finit;
  dgx = dginit;
  sty = 0.0;
  fy  = finit;
  dgy = dginit;

  //ls->step=ls->initstep;
  _step=_initstep;
  //for (i=0; i< ls->max_funcs; i++) 
  for (i=0; i< _max_funcs; i++) 
  {
    // Set min and max steps to correspond to the interval of uncertainty
    //if (mt->bracket) 
    if (_bracket) 
    {
      //ls->stepmin = PetscMin(stx,sty); 
      _stepmin = std::min(stx,sty); 
      //ls->stepmax = PetscMax(stx,sty); 
      _stepmax = std::max(stx,sty); 
    } 
    else 
    {
      //ls->stepmin = stx;
      _stepmin = stx;
      //ls->stepmax = ls->step + xtrapf * (ls->step - stx);
      _stepmax = _step + xtrapf * (_step - stx);
    }

    // Force the step to be within the bounds
    //ls->step = PetscMax(ls->step,ls->stepmin);
    _step = std::max(_step,_stepmin);
    //ls->step = PetscMin(ls->step,ls->stepmax);
    _step = std::min(_step,_stepmax);
  
    // If an unusual termination is to occur, then let step be the lowest
    //   point obtained thus far
    //if ((stx!=0) && (((mt->bracket) && (ls->step <= ls->stepmin || ls->step >= ls->stepmax)) ||
    //     ((mt->bracket) && (ls->stepmax - ls->stepmin <= ls->rtol * ls->stepmax)) ||
    //     ((ls->nfeval+ls->nfgeval) >= ls->max_funcs - 1) || (mt->infoc == 0))) 
    if ((stx!=0) && (((_bracket) && (_step <= _stepmin || _step >= _stepmax)) ||
         ((_bracket) && (_stepmax - _stepmin <= _rtol * _stepmax)) ||
         ((_nfeval+_nfgeval) >= _max_funcs - 1) || (_infoc == 0))) 
    {
      //ls->step = stx;
      _step = stx;
    }

    //VecCopy(x,mt->work);
    _work = x;
    //VecAXPY(mt->work,ls->step,s);   // W = X + step*S
    _work += _step * s; 

    //if (ls->bounded) 
    if (_bounded) 
    {
      //VecMedian(ls->lower, mt->work, ls->upper, mt->work);
      _work = median(_lower, _work, _upper);
    }
    //if (ls->usegts) 
    if (_usegts) 
    {
      //TaoLineSearchComputeObjectiveAndGTS(ls,mt->work,f,&dg); 
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " GTS COMPUTATION NOT SUPPORTED " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      exit(0);
      g_computed = false;
    } 
    else 
    {
      //TaoLineSearchComputeObjectiveAndGradient(ls,mt->work,f,g);
      ComputeObjectiveAndGradient(_work,f,g);

      //g_computed=PETSC_TRUE;
      g_computed = true;
      //if (ls->bounded) 
      if (_bounded)
      {
        //VecDot(g,x,&dg);
        //dg = g * x;
        dg = g.dot(x);
        //VecDot(g,mt->work,&dg2);
        //dg2 = g * _work;
        dg2 = g.dot(_work);
        //dg = (dg2 - dg)/ls->step;
        dg = (dg2 - dg) / _step;
      } 
      else 
      {
        //VecDot(g,s,&dg);
        //dg = g * s;
        dg = g.dot(s);
      }
    }

    if (0 == i) 
    {
      //ls->f_fullstep=*f;
      _f_fullstep = f;
    }

    //if (PetscIsInfOrNanReal(*f) || PetscIsInfOrNanReal(dg)) 
    if (isinf(f) || isnan(f) || isinf(dg) || isnan(dg)) 
    {
      // User provided compute function generated Not-a-Number, assume 
      // domain violation and set function value and directional
      // derivative to infinity.
      //*f = TAO_INFINITY;
      f = INF;
      //dg = TAO_INFINITY;
      dg = INF;
    }

    //ftest1 = finit + ls->step * dgtest;
    ftest1 = finit + _step * dgtest;
    //if (ls->bounded) 
    if (_bounded) 
    {
      //ftest2 = finit + ls->step * dgtest * ls->ftol;
      ftest2 = finit + _step * dgtest * _ftol;
    }
    // Convergence testing
    //if (((*f - ftest1 <= 1.0e-10 * PetscAbsReal(finit)) &&  (PetscAbsReal(dg) + ls->gtol*dginit <= 0.0))) 
    if (((f - ftest1 <= 1.0e-10 * fabs(finit)) &&  (fabs(dg) + _gtol*dginit <= 0.0))) 
    {
      //PetscInfo(ls, "Line search success: Sufficient decrease and directional deriv conditions hold\n");
      _reason = LINESEARCH_SUCCESS;
      break;
    }

    // Check Armijo if beyond the first breakpoint
    //if (ls->bounded && (*f <= ftest2) && (ls->step >= bstepmin2)) 
    if (_bounded && (f <= ftest2) && (_step >= bstepmin2)) 
    {
      //PetscInfo(ls,"Line search success: Sufficient decrease.\n");
      break;
    }

    // Checks for bad cases
    //if (((mt->bracket) && (ls->step <= ls->stepmin||ls->step >= ls->stepmax)) || (!mt->infoc)) 
    if (((_bracket) && (_step <= _stepmin || _step >= _stepmax)) || (!_infoc)) 
    {
      cout << " bracket: " << _bracket << endl;
      //PetscInfo(ls,"Rounding errors may prevent further progress.  May not be a step satisfying\n");
      //PetscInfo(ls,"sufficient decrease and curvature conditions. Tolerances may be too small.\n");
      _reason = LINESEARCH_HALTED_OTHER;
      break;
    }
    //if ((ls->step == ls->stepmax) && (*f <= ftest1) && (dg <= dgtest)) 
    if ((_step == _stepmax) && (f <= ftest1) && (dg <= dgtest)) 
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " step: " << _step << endl;
      cout << " stepmax: " << _stepmax << endl;
      cout << " f: " << f << endl;
      cout << " ftest1: " << ftest1 << endl;
      cout << " dg: " << dg << endl;
      cout << " dgtest: " << dgtest << endl;
      //PetscInfo1(ls,"Step is at the upper bound, stepmax (%G)\n",ls->stepmax);
      //ls->reason = TAOLINESEARCH_HALTED_UPPERBOUND;
      _reason = LINESEARCH_HALTED_UPPERBOUND;
      break;
    }
    //if ((ls->step == ls->stepmin) && (*f >= ftest1) && (dg >= dgtest)) 
    if ((_step == _stepmin) && (f >= ftest1) && (dg >= dgtest)) 
    {
      //PetscInfo1(ls,"Step is at the lower bound, stepmin (%G)\n",ls->stepmin);
      //ls->reason = TAOLINESEARCH_HALTED_LOWERBOUND;
      _reason = LINESEARCH_HALTED_LOWERBOUND;
      break;
    }
    //if ((mt->bracket) && (ls->stepmax - ls->stepmin <= ls->rtol*ls->stepmax))
    if ((_bracket) && (_stepmax - _stepmin <= _rtol*_stepmax))
    {
      //PetscInfo1(ls,"Relative width of interval of uncertainty is at most rtol (%G)\n",ls->rtol);
      //ls->reason = TAOLINESEARCH_HALTED_RTOL;
      _reason = LINESEARCH_HALTED_RTOL;
      break;
    }

    // In the first stage, we seek a step for which the modified function
    //  has a nonpositive value and nonnegative derivative
    //if ((stage1) && (*f <= ftest1) && (dg >= dginit * PetscMin(ls->ftol, ls->gtol))) 
    if ((stage1) && (f <= ftest1) && (dg >= dginit * std::min(_ftol, _gtol))) 
    {
      stage1 = 0;
    }

    // A modified function is used to predict the step only if we
    //   have not obtained a step for which the modified function has a 
    //   nonpositive function value and nonnegative derivative, and if a
    //   lower function value has been obtained but the decrease is not
    //   sufficient

   //if ((stage1) && (*f <= fx) && (*f > ftest1)) 
   if ((stage1) && (f <= fx) && (f > ftest1)) 
   {
      //fm   = *f - ls->step * dgtest;	// Define modified function 
      fm   = f - _step * dgtest;	// Define modified function 
      fxm  = fx - stx * dgtest;	      // and derivatives
      fym  = fy - sty * dgtest;
      dgm  = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      // if (dgxm * (ls->step - stx) >= 0.0)
      // Update the interval of uncertainty and compute the new step
      //Tao_mcstep(ls,&stx,&fxm,&dgxm,&sty,&fym,&dgym,&ls->step,&fm,&dgm);
      mcstep(stx,fxm,dgxm,sty,fym,dgym,_step,fm,dgm);
          
      fx  = fxm + stx * dgtest;	// Reset the function and 
      fy  = fym + sty * dgtest;	// gradient values 
      dgx = dgxm + dgtest; 
      dgy = dgym + dgtest; 
    } 
    else {
      // Update the interval of uncertainty and compute the new step
      //Tao_mcstep(ls,&stx,&fx,&dgx,&sty,&fy,&dgy,&ls->step,f,&dg);
      mcstep(stx,fx,dgx,sty,fy,dgy,_step,f,dg);
    }
  
    // Force a sufficient decrease in the interval of uncertainty
    if (_bracket) {
      //if (PetscAbsReal(sty - stx) >= 0.66 * width1) ls->step = stx + 0.5*(sty - stx);
      if (fabs(sty - stx) >= 0.66 * width1) _step = stx + 0.5*(sty - stx);
      width1 = width;
      //width = PetscAbsReal(sty - stx);
      width = fabs(sty - stx);
    }
  }
  //if ((ls->nfeval+ls->nfgeval) > ls->max_funcs) 
  if ((_nfeval+_nfgeval) > _max_funcs) 
  {
    //PetscInfo2(ls,"Number of line search function evals (%D) > maximum (%D)\n",(ls->nfeval+ls->nfgeval),ls->max_funcs);
    //ls->reason = TAOLINESEARCH_HALTED_MAXFCN;
    _reason = LINESEARCH_HALTED_MAXFCN;
  }

  // Finish computations
  //PetscInfo2(ls,"%D function evals in line search, step = %G\n",(ls->nfeval+ls->nfgeval),ls->step);
  
  // Set new solution vector and compute gradient if needed
  //VecCopy(mt->work,x);
  x = _work;

  // this should never fire - GTS is not supported
  if (!g_computed) 
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " UNSUPPORTED CODE PATH " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    exit(0);
    //TaoLineSearchComputeGradient(ls,mt->work,g);
  }

  //PetscFunctionReturn(0);
}
