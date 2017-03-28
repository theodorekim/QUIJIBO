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

#include "LINE_SEARCH.h"
#include <vector>
#include <cstdio>

using namespace std;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
LINE_SEARCH::LINE_SEARCH(void (*functionGradient)(const VECTOR&, Real&, VECTOR&))
{
  _bounded = 0;
  _max_funcs=30;
  _ftol = 0.0001;
  _gtol = 0.9;
  _rtol = 1.0e-10;
  _stepmin=1.0e-20;
  _stepmax=1.0e+20;
  _step=1.0;
  _nfeval=0;
  _ngeval=0;
  _nfgeval=0;

  // this appears to never be initialized?
  _usegts = false;

  //_ops->computeobjective=0;
  //_ops->computegradient=0;
  //_ops->computeobjectiveandgradient=0;
  //_ops->computeobjectiveandgts=0;
  //_ops->setup=0;
  //_ops->apply=0;
  //_ops->view=0;
  //_ops->setfromoptions=0;
  //_ops->reset=0;
  //_ops->destroy=0;
  _setupcalled=false;
  _usetaoroutines=false;

  _functionGradient = functionGradient;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void LINE_SEARCH::VecStepBoundInfo(VECTOR& x, VECTOR& xl, VECTOR& xu, VECTOR& dx, Real& boundmin, Real& wolfemin, Real& boundmax)
{
  //PetscErrorCode ierr;
  //PetscInt n,i;
  int n,i;
  //PetscReal *x,*xl,*xu,*dx;
  //PetscReal t;
  Real t;
  //PetscReal localmin=1.0e300,localwolfemin=1.0e300,localmax=0;
  Real localmin=1.0e300,localwolfemin=1.0e300,localmax=0;
 
  //PetscFunctionBegin;
  //PetscValidHeaderSpecific(X,VEC_CLASSID,1);
  //PetscValidHeaderSpecific(XL,VEC_CLASSID,2);
  //PetscValidHeaderSpecific(XU,VEC_CLASSID,3);
  //PetscValidHeaderSpecific(DX,VEC_CLASSID,4);

  //ierr=VecGetArray(X,&x);CHKERRQ(ierr);
  //ierr=VecGetArray(XL,&xl);CHKERRQ(ierr);
  //ierr=VecGetArray(XU,&xu);CHKERRQ(ierr);
  //ierr=VecGetArray(DX,&dx);CHKERRQ(ierr);
  //ierr = VecGetLocalSize(X,&n);CHKERRQ(ierr);
  n = x.size();
  for (i=0;i<n;i++){
    if (dx[i]>0){
      t=(xu[i]-x[i])/dx[i];
      //localmin = PetscMin(t,localmin);
      localmin = std::min(t,localmin);
      if (localmin>0){
	      //localwolfemin = PetscMin(t,localwolfemin);
	      localwolfemin = std::min(t,localwolfemin);
      }
      //localmax = PetscMax(t,localmax);
      localmax = std::max(t,localmax);
    } else if (dx[i]<0){
      t=(xl[i]-x[i])/dx[i];
      //localmin = PetscMin(t,localmin);
      localmin = std::min(t,localmin);
      if (localmin>0){
      	//localwolfemin = PetscMin(t,localwolfemin);
      	localwolfemin = std::min(t,localwolfemin);
      }
      //localmax = PetscMax(t,localmax);
      localmax = std::max(t,localmax);
    }
  }
  /*
  ierr=VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr=VecRestoreArray(XL,&xl);CHKERRQ(ierr);
  ierr=VecRestoreArray(XU,&xu);CHKERRQ(ierr);
  ierr=VecRestoreArray(DX,&dx);CHKERRQ(ierr);
  ierr=PetscObjectGetComm((PetscObject)X,&comm);CHKERRQ(ierr);
  
  if (boundmin){ ierr = MPI_Allreduce(&localmin,boundmin,1,MPIU_REAL,MPIU_MIN,comm);CHKERRQ(ierr);}
  if (wolfemin){ ierr = MPI_Allreduce(&localwolfemin,wolfemin,1,MPIU_REAL,MPIU_MIN,comm);CHKERRQ(ierr);}
  if (boundmax) { ierr = MPI_Allreduce(&localmax,boundmax,1,MPIU_REAL,MPIU_MAX,comm);CHKERRQ(ierr);}

  ierr = PetscInfo3(X,"Step Bound Info: Closest Bound: %G, Wolfe: %G, Max: %G \n",*boundmin,*wolfemin,*boundmax); CHKERRQ(ierr);
  PetscFunctionReturn(0);  
  */
  boundmin = localmin;
  wolfemin = localwolfemin;
  boundmax = localmax;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//PetscErrorCode TaoLineSearchComputeObjectiveAndGradient(TaoLineSearch ls, Vec x, PetscReal *f, Vec g) 
void LINE_SEARCH::ComputeObjectiveAndGradient(VECTOR& x, Real& f, VECTOR& g)
{
  //PetscErrorCode ierr;
  //PetscFunctionBegin;
  //PetscValidHeaderSpecific(ls,TAOLINESEARCH_CLASSID,1);
  //PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  //PetscValidPointer(f,3);
  //PetscValidHeaderSpecific(g,VEC_CLASSID,4);
  //PetscCheckSameComm(ls,1,x,2);
  //PetscCheckSameComm(ls,1,g,4);

  if (_usetaoroutines) 
  {
    // this does not appear to fire in the cases I care about
    //TaoComputeObjectiveAndGradient(ls->taosolver,x,f,g); 
  } 
  else 
  {
    //ierr = PetscLogEventBegin(TaoLineSearch_EvalEvent,ls,0,0,0); CHKERRQ(ierr);
    //if (!ls->ops->computeobjective && !ls->ops->computeobjectiveandgradient) {
    //  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Line Search does not have objective function set");
    //}
    //if (!ls->ops->computegradient && !ls->ops->computeobjectiveandgradient) {
    //  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Line Search does not have gradient function set");
    //}

    //PetscStackPush("TaoLineSearch user objective/gradient routine"); 
    //CHKMEMQ;
    //if (ls->ops->computeobjectiveandgradient) 
    //{
    //  (*ls->ops->computeobjectiveandgradient)(ls,x,f,g,ls->userctx_funcgrad);
    //  this one always fires

    // TODO: need to make this generic
    //FormScaledOnlyFunctionGradient(x, optimize3D, f, g);
    //FormScaledFunctionGradient(x, optimize3D, f, g);
    //_functionGradient(x, optimize3D, f, g);
    _functionGradient(x, f, g);

    //} 
    //else 
    //{
    //  (*ls->ops->computeobjective)(ls,x,f,ls->userctx_func);
    //  (*ls->ops->computegradient)(ls,x,g,ls->userctx_grad);
    //}
    //CHKMEMQ;
    //PetscStackPop;
    //ierr = PetscLogEventEnd(TaoLineSearch_EvalEvent,ls,0,0,0); CHKERRQ(ierr);
  }
  //ierr = PetscInfo1(ls,"TaoLineSearch Function evaluation: %14.12e\n",*f);CHKERRQ(ierr);    
  //ls->nfgeval++;
  //PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void LINE_SEARCH::VecBoundGradientProjection(const VECTOR& G, const VECTOR& X, const VECTOR& XL, const VECTOR& XU, VECTOR& GP)
{
  int n,i;
  Real xval,gpval;

  // Project variables at the lower and upper bound
  n = X.size();
  assert(G.size() == n);
  assert(XU.size() == n);
  assert(GP.size() == n);

  for (i=0; i < n; ++i){
    gpval = G[i]; xval = X[i]; 

    if (gpval>0 && xval<=XL[i]){
      gpval = 0;
    } else if (gpval<0 && xval>=XU[i]){
      gpval = 0;
    }
    GP[i] = gpval;
  }
}
/*
{
  int n,i;
  const Real *xptr,*xlptr,*xuptr,*gptr;
  Real* gpptr;
  Real xval,gpval;

  // Project variables at the lower and upper bound
  n = X.size();
  assert(G.size() == n);
  assert(XU.size() == n);
  assert(GP.size() == n);

  xptr = X.dataConst();
  xlptr = XL.dataConst();
  xuptr = XU.dataConst();
  gptr = G.dataConst();
  gpptr = GP.data();

  for (i=0; i < n; ++i){
    gpval = gptr[i]; xval = xptr[i]; 

    if (gpval>0 && xval<=xlptr[i]){
      gpval = 0;
    } else if (gpval<0 && xval>=xuptr[i]){
      gpval = 0;
    }
    gpptr[i] = gpval;
  }
}
*/

////////////////////////////////////////////////////////////////////////////////////////
// return a vector with the median values from all three vectors
////////////////////////////////////////////////////////////////////////////////////////
VECTOR LINE_SEARCH::median(const VECTOR& v0, const VECTOR& v1, const VECTOR& v2)
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
// Get everything setup for the call to LineSearchImplementation
////////////////////////////////////////////////////////////////////////////////////////
void LINE_SEARCH::LineSearchApply(VECTOR& x, Real& f, VECTOR& g, VECTOR& s)
{
  //PetscErrorCode ierr;
  //PetscViewer viewer;
  //PetscInt low1,low2,low3,high1,high2,high3;
  //int low1,low2,low3,high1,high2,high3;
  //PetscBool flg;
  //bool flg;
  //char filename[PETSC_MAX_PATH_LEN];

  //PetscFunctionBegin;
  //*reason = TAOLINESEARCH_CONTINUE_ITERATING;
  _reason = LINESEARCH_CONTINUE_ITERATING;
  //PetscValidHeaderSpecific(ls,TAOLINESEARCH_CLASSID,1);
  //PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  //PetscValidScalarPointer(f,3);
  //PetscValidHeaderSpecific(g,VEC_CLASSID,4);
  //PetscValidHeaderSpecific(s,VEC_CLASSID,5);
  //PetscValidPointer(reason,7);
  //PetscCheckSameComm(ls,1,x,2);
  //PetscCheckSameTypeAndComm(x,2,g,4);
  //PetscCheckSameTypeAndComm(x,2,s,5);
  //ierr = VecGetOwnershipRange(x, &low1, &high1); CHKERRQ(ierr);
  //ierr = VecGetOwnershipRange(g, &low2, &high2); CHKERRQ(ierr);
  //ierr = VecGetOwnershipRange(s, &low3, &high3); CHKERRQ(ierr);
  //if ( low1!= low2 || low1!= low3 || high1!= high2 || high1!= high3) {
  // SETERRQ(PETSC_COMM_SELF,1,"InCompatible vector local lengths");
  //}

  //if (ls->stepdirection) {
  // ierr = PetscObjectDereference((PetscObject)(s)); CHKERRQ(ierr);
  //}
  //ls->stepdirection = s;
  _stepdirection = s;
  //ierr = PetscObjectReference((PetscObject)s); CHKERRQ(ierr);

  // this checks if the function pointer has been defined, which is 
  // guaranteed by the constructor
  //ierr = TaoLineSearchSetUp(ls); CHKERRQ(ierr);
  //if (!ls->ops->apply) {
  //  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Line Search Object does not have 'apply' routine");
  //}
  //ls->nfeval=0;
  _nfeval = 0;
  //ls->ngeval=0;
  _ngeval = 0;
  //ls->nfgeval=0;
  _nfgeval = 0;

  // Check parameter values
  //if (ls->ftol < 0.0) 
  if (_ftol < 0.0) 
  {
    //ierr = PetscInfo1(ls,"Bad Line Search Parameter: ftol (%G) < 0\n",ls->ftol); CHKERRQ(ierr);
    printf("Bad Line Search Parameter: ftol (%G) < 0\n",(double)_ftol);
    //*reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }
  //if (ls->rtol < 0.0) 
  if (_rtol < 0.0) 
  {
    //ierr = PetscInfo1(ls,"Bad Line Search Parameter: rtol (%G) < 0\n",ls->rtol); CHKERRQ(ierr);
    printf("Bad Line Search Parameter: rtol (%G) < 0\n",(double)_rtol);
    //*reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }

  //if (ls->gtol < 0.0) 
  if (_gtol < 0.0) 
  {
    //ierr = PetscInfo1(ls,"Bad Line Search Parameter: gtol (%G) < 0\n",ls->gtol); CHKERRQ(ierr);
    printf("Bad Line Search Parameter: gtol (%G) < 0\n",(double)_gtol);
    //*reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }
  //if (ls->stepmin < 0.0) 
  if (_stepmin < 0.0) 
  {
    //ierr = PetscInfo1(ls,"Bad Line Search Parameter: stepmin (%G) < 0\n",ls->stepmin); CHKERRQ(ierr);
    printf("Bad Line Search Parameter: stepmin (%G) < 0\n",(double)_stepmin);
    //*reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }
  //if (ls->stepmax < ls->stepmin) 
  if (_stepmax < _stepmin) 
  {
    //ierr = PetscInfo2(ls,"Bad Line Search Parameter: stepmin (%G) > stepmax (%G)\n",ls->stepmin,ls->stepmax); CHKERRQ(ierr);
    printf("Bad Line Search Parameter: stepmin (%G) > stepmax (%G)\n",(double)_stepmin,(double)_stepmax);
    //*reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }
  //if (ls->max_funcs < 0) 
  if (_max_funcs < 0) 
  {
    //ierr = PetscInfo1(ls,"Bad Line Search Parameter: max_funcs (%D) < 0\n",ls->max_funcs); CHKERRQ(ierr);
    printf("Bad Line Search Parameter: max_funcs (%i) < 0\n",_max_funcs);
    //*reason=TAOLINESEARCH_FAILED_BADPARAMETER;
    _reason = LINESEARCH_FAILED_BADPARAMETER;
  }
  //if (PetscIsInfOrNanReal(*f)) 
  if (isinf(f) || isnan(f)) 
  {
    //ierr = PetscInfo1(ls,"Initial Line Search Function Value is Inf or Nan (%G)\n",*f); CHKERRQ(ierr);
    printf("Initial Line Search Function Value is Inf or Nan (%G)\n",(double)f);
    //*reason=TAOLINESEARCH_FAILED_INFORNAN;
    _reason = LINESEARCH_FAILED_INFORNAN;
  }

  // TK: Just do the copy every time. If this becomes a bottleneck later, 
  // we can revisit.
  //if (x != ls->start_x) 
  //if (((x - _start_x).norm2() > 1e-7))
  {
    //ierr = PetscObjectReference((PetscObject)x);
    //if (ls->start_x)
    //  ierr = VecDestroy(&ls->start_x); CHKERRQ(ierr);
    _start_x = x;
  }

  //ierr = PetscLogEventBegin(TaoLineSearch_ApplyEvent,ls,0,0,0); CHKERRQ(ierr);
  // Make the implementation-specific linesearch call (finally)
  //ierr = (*ls->ops->apply)(ls,x,f,g,s); CHKERRQ(ierr);
  LineSearchImplementation(x,f,g,s);
  //ierr = PetscLogEventEnd(TaoLineSearch_ApplyEvent, ls, 0,0,0); CHKERRQ(ierr);
  //*reason=ls->reason;
  //ls->new_f = *f;
  _new_f = f;

  // TK: pass back the step length
  // don't need to do this -- it's stored locally anyway if we need it later
  //steplength = _step;

  //ierr = PetscOptionsGetString(((PetscObject)ls)->prefix,"-tao_ls_view",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  //if (ls->viewls && !PetscPreLoadingOn) {
  //  ierr = PetscViewerASCIIOpen(((PetscObject)ls)->comm,filename,&viewer); CHKERRQ(ierr);
  //  ierr = TaoLineSearchView(ls,viewer); CHKERRQ(ierr);
  //  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  //}
  //PetscFunctionReturn(0);
}
