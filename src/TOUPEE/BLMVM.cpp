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

#include "BLMVM.h"
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
BLMVM::BLMVM()
{
  _allocated = false;
  _lm = 5;
  _S.resize(0);
  _Y.resize(0);
  _delta_min = 1e-7;
  _delta_max = 100;

  _eps=0.0;
  _limitType = MatLMVM_Limit_None;
  _scaleType = MatLMVM_Scale_Broyden;
  _rScaleType = MatLMVM_Rescale_Scalar;
  _s_alpha = 1.0;
  _r_alpha = 1.0;
  _r_beta = 0.5;
  _mu = 1.0;
  _nu = 100.0;

  _phi = 0.125;		

  _scalar_history = 1;
  _rescale_history = 1;

  //  Complete configuration
  _rescale_history = min(_rescale_history, _lm);

  //PetscMalloc((ctx->lm+1)*sizeof(PetscReal),(void**)&ctx->rho); 
  _rho  = VECTOR(_lm + 1);
  //PetscMalloc((ctx->lm+1)*sizeof(PetscReal),(void**)&ctx->beta); 
  _beta = VECTOR(_lm + 1);

  int nhistory = max(_scalar_history,1);
  //PetscMalloc(nhistory*sizeof(PetscReal),(void**)&ctx->yy_history); 
  _yy_history = VECTOR(nhistory);
  //PetscMalloc(nhistory*sizeof(PetscReal),(void**)&ctx->ys_history);
  _ys_history = VECTOR(nhistory);
  //PetscMalloc(nhistory*sizeof(PetscReal),(void**)&ctx->ss_history); 
  _ss_history = VECTOR(nhistory);

  nhistory = max(_rescale_history,1);
  //PetscMalloc(nhistory*sizeof(PetscReal),(void**)&ctx->yy_rhistory);
  _yy_rhistory = VECTOR(nhistory);
  //PetscMalloc(nhistory*sizeof(PetscReal),(void**)&ctx->ys_rhistory);
  _ys_rhistory = VECTOR(nhistory);
  //PetscMalloc(nhistory*sizeof(PetscReal),(void**)&ctx->ss_rhistory); 
  _ss_rhistory = VECTOR(nhistory);

  //  Finish initializations
  _lmnow = 0;
  _iter = 0;
  _nupdates = 0;
  _nrejects = 0;
  _delta = 1.0;

  _Gprev.setZero();
  _Xprev.setZero();

  _scale.setZero(); 
  _useScale = false;

  //_H0 = 0;
  _useDefaultH0 = true;
  
  //MatCreateShell(comm, n, n, N, N, ctx, A); CHKERRQ(ierr);
  //MatShellSetOperation(*A,MATOP_DESTROY,(void(*)(void))MatDestroy_LMVM);
  //MatShellSetOperation(*A,MATOP_VIEW,(void(*)(void))MatView_LMVM);
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void BLMVM::MatReset(MATRIX& M)
{
  int i;

  // no need to do this
  //if (_Gprev) {
  //  PetscObjectDereference((PetscObject)_Gprev);
  //}
  //if (_Xprev) {
  //  PetscObjectDereference((PetscObject)_Xprev);
  //}
  _Gprev = _Y[_lm];
  assert(_lm < (int)_S.size());
  _Xprev = _S[_lm];

  // no need to do this
  //PetscObjectReference((PetscObject)_Gprev);
  //PetscObjectReference((PetscObject)_Xprev);
  //
  for (i=0; i<_lm; ++i) {
    _rho[i] = 0.0;
  }
  _rho[0] = 1.0;
  
  //  Set the scaling and diagonal scaling matrix
  switch(_scaleType) {
    case MatLMVM_Scale_None:
      _sigma = 1.0;
      break;
    case MatLMVM_Scale_Scalar:
      _sigma = _delta;
      break;
    case MatLMVM_Scale_Broyden:
      _D.setConstant(_delta);
      break;
  }

  _iter=0;
  _nupdates=0;
  _lmnow=0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void BLMVM::MatAllocateVectors(MATRIX& m, VECTOR& v)
{
  //  Perform allocations
  _S.resize(_lm + 1);
  _Y.resize(_lm + 1);
  for (int x = 0; x < _lm + 1; x++)
  {
    _S[x] = v;
    _Y[x] = v;
  }
  //VecDuplicateVecs(v,ctx->lm+1,&ctx->S);
  //VecDuplicateVecs(v,ctx->lm+1,&ctx->Y);

  _D = v;
  _U = v;
  _V = v;
  _W = v;
  _P = v;
  _Q = v;

  _allocated = true;
}

////////////////////////////////////////////////////////////////////////////////////////
// x = solution
// g = gradient
////////////////////////////////////////////////////////////////////////////////////////
void BLMVM::Update(MATRIX& M, VECTOR& x, VECTOR& g)
{
  //MatLMVMCtx *ctx;
  Real rhotemp, rhotol;
  Real y0temp, s0temp;
  Real yDy, yDs, sDs;
  Real sigmanew, denom;
  int i;
  //bool same;  // not used?
  Real yy_sum=0.0, ys_sum=0.0, ss_sum=0.0;

  //PetscFunctionBegin;

  // make sure these are vectors
  //PetscValidHeaderSpecific(x,VEC_CLASSID,2); 
  //PetscValidHeaderSpecific(g,VEC_CLASSID,3);
  
  // make sure the matrix can support LMVM
  //PetscObjectTypeCompare((PetscObject)M,MATSHELL,&same);
  
  //if (!same) {SETERRQ(PETSC_COMM_SELF,1,"Matrix M is not type MatLMVM");}

  // allocate the LMVM matrix
  //MatShellGetContext(M,(void**)&ctx);
  //if (!ctx->allocated) {
  //    MatLMVMAllocateVectors(M, x); 
  //}
  if (!_allocated)
  {
    cout << " NOT ALLOCATED " << endl;
    MatAllocateVectors(M, x); 
  }

  if (_iter == 0) 
  {
    cout << " RESETTING " << endl;
    MatReset(M);
  } 
  else 
  {
    //VecAYPX(ctx->Gprev,-1.0,g);
    _Gprev = g - _Gprev;
    //VecAYPX(ctx->Xprev,-1.0,x);
    _Xprev = x - _Xprev;

    //VecDot(ctx->Gprev,ctx->Xprev,&rhotemp);
    rhotemp = _Gprev.dot(_Xprev);
    //VecDot(ctx->Gprev,ctx->Gprev,&y0temp);
    y0temp = _Gprev.dot(_Gprev);

    //rhotol = ctx->eps * y0temp;
    rhotol = _eps * y0temp;
    
    if (rhotemp > rhotol) 
    {
      //++ctx->nupdates;
      _nupdates++;

      //ctx->lmnow = PetscMin(ctx->lmnow+1, ctx->lm);
      _lmnow = min(_lmnow+1, _lm);
      //for (i = ctx->lm-1; i >= 0; --i) 
      for (i = _lm-1; i >= 0; --i) {
        //ctx->S[i+1] = ctx->S[i];
        assert(i+1 < (int)_S.size());
        assert(i+1 < (int)_Y.size());
        assert(i >= 0);
        _S[i+1] = _S[i];
        //ctx->Y[i+1] = ctx->Y[i];
        _Y[i+1] = _Y[i];
        //ctx->rho[i+1] = ctx->rho[i];
        _rho[i+1] = _rho[i];
      }
      //ctx->S[0] = ctx->Xprev;
      _S[0] = _Xprev;
      //ctx->Y[0] = ctx->Gprev;
      _Y[0] = _Gprev;

      //PetscObjectReference((PetscObject)ctx->S[0]);
      //PetscObjectReference((PetscObject)ctx->Y[0]);
      //ctx->rho[0] = 1.0 / rhotemp;
      _rho[0] = 1.0 / rhotemp;

      //  Compute the scaling
      switch(_scaleType) {
      case MatLMVM_Scale_None:
        break;

      case MatLMVM_Scale_Scalar:
        //  Compute s^T s
        //VecDot(ctx->Xprev,ctx->Xprev,&s0temp);
        s0temp = _Xprev.dot(_Xprev);

        //  Scalar is positive; safeguards are not required.

        //  Save information for scalar scaling
        //ctx->yy_history[(ctx->nupdates - 1) % ctx->scalar_history] = y0temp;
        _yy_history[(_nupdates - 1) % _scalar_history] = y0temp;
        //ctx->ys_history[(ctx->nupdates - 1) % ctx->scalar_history] = rhotemp;
        _ys_history[(_nupdates - 1) % _scalar_history] = rhotemp;
        //ctx->ss_history[(ctx->nupdates - 1) % ctx->scalar_history] = s0temp;
        _ss_history[(_nupdates - 1) % _scalar_history] = s0temp;

        //  Compute summations for scalar scaling
        yy_sum = 0; //  No safeguard required; y^T y > 0
        ys_sum = 0; //  No safeguard required; y^T s > 0
        ss_sum = 0; //  No safeguard required; s^T s > 0
        //for (i = 0; i < PetscMin(ctx->nupdates, ctx->scalar_history); ++i)
        for (i = 0; i < min(_nupdates, _scalar_history); ++i) {
          //yy_sum += ctx->yy_history[i];
          yy_sum += _yy_history[i];
          //ys_sum += ctx->ys_history[i];
          ys_sum += _ys_history[i];
          //ss_sum += ctx->ss_history[i];
          ss_sum += _ss_history[i];
        }

        //if (0.0 == ctx->s_alpha) 
        if (0.0 == _s_alpha) 
        {
          //  Safeguard ys_sum
          if (0.0 == ys_sum) {
            ys_sum = ZERO_SAFEGUARD;
          }

          sigmanew = ss_sum / ys_sum;
        }
        //else if (1.0 == ctx->s_alpha) 
        else if (1.0 == _s_alpha) 
        {
          //  Safeguard yy_sum
          if (0.0 == yy_sum) {
            yy_sum = ZERO_SAFEGUARD;
          }

          sigmanew = ys_sum / yy_sum;
        }
        else {
          //denom = 2*ctx->s_alpha*yy_sum;
          denom = 2*_s_alpha*yy_sum;

          //  Safeguard denom
          if (0.0 == denom) {
            denom = ZERO_SAFEGUARD;
          }

          //sigmanew = ((2*ctx->s_alpha-1)*ys_sum + 
          //          PetscSqrtScalar((2*ctx->s_alpha-1)*(2*ctx->s_alpha-1)*ys_sum*ys_sum - 
          //                           4*(ctx->s_alpha)*(ctx->s_alpha-1)*yy_sum*ss_sum)) / denom;
          sigmanew = ((2*_s_alpha-1)*ys_sum + 
                                sqrt((2*_s_alpha-1)*(2*_s_alpha-1)*ys_sum*ys_sum - 
                                     4*(_s_alpha)*(_s_alpha-1)*yy_sum*ss_sum)) / denom;
        }

        //switch(ctx->limitType) 
        switch(_limitType) 
        {
          case MatLMVM_Limit_Average:
            //if (1.0 == ctx->mu) 
            if (1.0 == _mu) 
            {
              //ctx->sigma = sigmanew;
              _sigma = sigmanew;
            }
            //else if (ctx->mu) 
            else if (_mu) 
            {
              //ctx->sigma = ctx->mu * sigmanew + (1.0 - ctx->mu) * ctx->sigma;
              _sigma = _mu * sigmanew + (1.0 - _mu) * _sigma;
            }
          break;
          case MatLMVM_Limit_Relative:
            //if (ctx->mu) 
            if (_mu) 
            {
              //ctx->sigma = TaoMid((1.0 - ctx->mu) * ctx->sigma, sigmanew, (1.0 + ctx->mu) * ctx->sigma);
              _sigma = mid((1.0 - _mu) * _sigma, sigmanew, (1.0 + _mu) * _sigma);
            }
            break;

          case MatLMVM_Limit_Absolute:
            //if (ctx->nu) 
            if (_nu) 
            {
              //ctx->sigma = TaoMid(ctx->sigma - ctx->nu, sigmanew, ctx->sigma + ctx->nu);
              _sigma = mid(_sigma - _nu, sigmanew, _sigma + _nu);
            }
            break;

          default:
            //ctx->sigma = sigmanew;
            _sigma = sigmanew;
            break;
        }
        break;

      case MatLMVM_Scale_Broyden:
        //  Original version
        //  Combine DFP and BFGS

        //  This code appears to be numerically unstable.  We use the
        //  original version because this was used to generate all of
        //  the data and because it may be the least unstable of the
        //  bunch.

        //  P = Q = inv(D);
        //VecCopy(ctx->D,ctx->P);
        _P = _D;
        //VecReciprocal(ctx->P);
        //_P = reciprocal(_P);
        _P = VECTOR::Ones(_P.size()).array() / _P.array();
        //VecCopy(ctx->P,ctx->Q);
        _Q = _P;

        //  V = y*y
        //VecPointwiseMult(ctx->V,ctx->Gprev,ctx->Gprev);
        //_V = _Gprev.pointwise(_Gprev);
        _V = _Gprev.array() * _Gprev.array();

        //  W = inv(D)*s 
        //VecPointwiseMult(ctx->W,ctx->Xprev,ctx->P);
        //_W = _Xprev.pointwise(_P);
        _W = _Xprev.array() * _P.array();

        //VecDot(ctx->W,ctx->Xprev,&sDs);
        //sDs = _W * _Xprev;
        sDs = _W.dot( _Xprev);

        //  Safeguard rhotemp and sDs 
        if (0.0 == rhotemp) {
          rhotemp = ZERO_SAFEGUARD;
        }

        if (0.0 == sDs) {
          sDs = ZERO_SAFEGUARD;
        }

        //if (1.0 != ctx->phi) {
        if (1.0 != _phi) 
        {
          //  BFGS portion of the update
          //  U = (inv(D)*s)*(inv(D)*s)
          //VecPointwiseMult(ctx->U,ctx->W,ctx->W);
          //_U = _W.pointwise(_W);
          _U = _W.array() * _W.array();

          //  Assemble
          //VecAXPY(ctx->P,1.0/rhotemp,ctx->V);
          _P += 1.0/rhotemp * _V;
          //VecAXPY(ctx->P,-1.0/sDs,ctx->U);
          _P += -1.0/sDs * _U;
        }

        //if (0.0 != ctx->phi) 
        if (0.0 != _phi) 
        {
          //  DFP portion of the update
          //  U = inv(D)*s*y
          //VecPointwiseMult(ctx->U, ctx->W, ctx->Gprev);
          //_U = _W.pointwise(_Gprev);
          _U = _W.array() * _Gprev.array();

          //  Assemble
          //VecAXPY(ctx->Q,1.0/rhotemp + sDs/(rhotemp*rhotemp), ctx->V);
          _Q += (1.0/rhotemp + sDs/(rhotemp*rhotemp)) * _V;
          //VecAXPY(ctx->Q,-2.0/rhotemp,ctx->U);
          _Q += -2.0/rhotemp * _U;
        }

        //if (0.0 == ctx->phi) 
        if (0.0 == _phi) 
        {
          //VecCopy(ctx->P,ctx->U);
          _U = _P;
        }
        //else if (1.0 == ctx->phi) 
        else if (1.0 == _phi) 
        {
          //VecCopy(ctx->Q,ctx->U);
          _U = _Q;
        }
        else 
        {
          //  Broyden update U=(1-phi)*P + phi*Q 
          //VecCopy(ctx->Q,ctx->U);
          _U = _Q;
          //VecAXPBY(ctx->U,1.0-ctx->phi, ctx->phi, ctx->P);
          _U = (1.0-_phi)* _P + _phi * _U;
        }

        //  Obtain inverse and ensure positive definite
        //VecReciprocal(ctx->U);
        _U = VECTOR::Ones(_U.size()).array() / _U.array();
        //VecAbs(ctx->U);
        //_U.fabs();;
        _U = _U.array().abs();

        //switch(ctx->rScaleType) 
        switch(_rScaleType) 
        {
          case MatLMVM_Rescale_None:
            break;

          case MatLMVM_Rescale_Scalar:
          case MatLMVM_Rescale_GL:
            //if (ctx->rScaleType == MatLMVM_Rescale_GL) 
            if (_rScaleType == MatLMVM_Rescale_GL) 
            {
              //  Gilbert and Lemarachal use the old diagonal
              //VecCopy(ctx->D,ctx->P);
              _P = _D;
            }
            else {
              //  The default version uses the current diagonal
              //VecCopy(ctx->U,ctx->P);
              _P = _U;
            }

            //  Compute s^T s
            //VecDot(ctx->Xprev,ctx->Xprev,&s0temp);
            //s0temp = _Xprev * _Xprev;
            s0temp = _Xprev.dot(_Xprev);

            //  Save information for special cases of scalar rescaling
            //ctx->yy_rhistory[(ctx->nupdates - 1) % ctx->rescale_history] = y0temp;
            _yy_rhistory[(_nupdates - 1) % _rescale_history] = y0temp;
            //ctx->ys_rhistory[(ctx->nupdates - 1) % ctx->rescale_history] = rhotemp;
            _ys_rhistory[(_nupdates - 1) % _rescale_history] = rhotemp;
            //ctx->ss_rhistory[(ctx->nupdates - 1) % ctx->rescale_history] = s0temp;
            _ss_rhistory[(_nupdates - 1) % _rescale_history] = s0temp;

            //if (0.5 == ctx->r_beta) 
            if (0.5 == _r_beta) 
            {
              //if (1 == PetscMin(ctx->nupdates, ctx->rescale_history)) 
              if (1 == min(_nupdates, _rescale_history)) 
              {
                //VecPointwiseMult(ctx->V,ctx->Y[0],ctx->P);
                //_V = _Y[0].pointwise(_P);
                _V = _Y[0].array() * _P.array();
                //VecDot(ctx->V,ctx->Y[0],&yy_sum);
                //yy_sum = _V * _Y[0];
                yy_sum = _V.dot(_Y[0]);
              
                //VecPointwiseDivide(ctx->W,ctx->S[0],ctx->P);
                //_W = _S[0].pointwiseDivide(_P);
                _W = _S[0].array() / _P.array();
                //VecDot(ctx->W,ctx->S[0],&ss_sum);
                //ss_sum = _W * _S[0];
                ss_sum = _W.dot(_S[0]);

                //ys_sum = ctx->ys_rhistory[0];
                ys_sum = _ys_rhistory[0];
              }
              else 
              {
                //VecCopy(ctx->P,ctx->Q);
                _Q = _P;
                //VecReciprocal(ctx->Q);
                //_Q = _Q.reciprocal();
                _Q = VECTOR::Ones(_Q.size()).array() / _Q.array();

                //  Compute summations for scalar scaling
                yy_sum = 0; //  No safeguard required
                ys_sum = 0; //  No safeguard required
                ss_sum = 0; //  No safeguard required
                //for (i = 0; i < PetscMin(ctx->nupdates, ctx->rescale_history); ++i) 
                for (i = 0; i < min(_nupdates, _rescale_history); ++i) 
                {
                  //VecPointwiseMult(ctx->V,ctx->Y[i],ctx->P);
                  //_V = _Y[i].pointwise(_P);
                  _V = _Y[i].array() * _P.array();
                  //VecDot(ctx->V,ctx->Y[i],&yDy);
                  //yDy = _V * _Y[i];
                  yDy = _V.dot(_Y[i]);
                  //yy_sum += yDy;
                  yy_sum += yDy;

                  //VecPointwiseMult(ctx->W,ctx->S[i],ctx->Q);
                  //_W = _S[i].pointwise(_Q);
                  //_W = pointwise(_S[i],_Q);
                  _W = _S[i].array() * _Q.array();
                  //VecDot(ctx->W,ctx->S[i],&sDs);
                  //sDs = _W * _S[i];
                  sDs = _W.dot(_S[i]);
                  //ss_sum += sDs;
                  ss_sum += sDs;
                  //ys_sum += ctx->ys_rhistory[i];
                  ys_sum += _ys_rhistory[i];
                }
              }
            }
            //else if (0.0 == ctx->r_beta) 
            else if (0.0 == _r_beta) 
            {
              //if (1 == PetscMin(ctx->nupdates, ctx->rescale_history)) 
              if (1 == min(_nupdates, _rescale_history)) 
              {
                //  Compute summations for scalar scaling
                //VecPointwiseDivide(ctx->W,ctx->S[0],ctx->P);
                //_W = _S[0].pointwiseDivide(_P);
                _W = _S[0].array() / _P.array();

                //VecDot(ctx->W, ctx->Y[0], &ys_sum);
                //ys_sum = _W * _Y[0];
                ys_sum = _W.dot(_Y[0]);
                //VecDot(ctx->W, ctx->W, &ss_sum);
                //ss_sum = _W * _W;
                ss_sum = _W.dot(_W);
                //yy_sum += ctx->yy_rhistory[0];
                yy_sum += _yy_rhistory[0];
              }
              else 
              {
                //VecCopy(ctx->Q, ctx->P);
                _P = _Q;
                //VecReciprocal(ctx->Q);
                //_Q = _Q.reciprocal();
                _Q = VECTOR::Ones(_Q.size()).array() / _Q.array();

                //  Compute summations for scalar scaling
                yy_sum = 0; //  No safeguard required
                ys_sum = 0; //  No safeguard required
                ss_sum = 0; //  No safeguard required
                //for (i = 0; i < PetscMin(ctx->nupdates, ctx->rescale_history); ++i) 
                for (i = 0; i < min(_nupdates, _rescale_history); ++i) 
                {
                  //VecPointwiseMult(ctx->W, ctx->S[i], ctx->Q);
                  //_W = _S[i].pointwise(_Q);
                  _W = _S[i].array() * _Q.array();
                  //VecDot(ctx->W, ctx->Y[i], &yDs);
                  //yDs = _W * _Y[i];
                  yDs = _W.dot(_Y[i]);
                  ys_sum += yDs;

                  //VecDot(ctx->W, ctx->W, &sDs);
                  //sDs = _W * _W;
                  sDs = _W.dot(_W);
                  ss_sum += sDs;
        
                  //yy_sum += ctx->yy_rhistory[i];
                  yy_sum += _yy_rhistory[i];
                }
              }
            }
            //else if (1.0 == ctx->r_beta) 
            else if (1.0 == _r_beta) 
            {
              //  Compute summations for scalar scaling
              yy_sum = 0; //  No safeguard required
              ys_sum = 0; //  No safeguard required
              ss_sum = 0; //  No safeguard required
              //for (i = 0; i < PetscMin(ctx->nupdates, ctx->rescale_history); ++i) 
              for (i = 0; i < min(_nupdates, _rescale_history); ++i) 
              {
                //VecPointwiseMult(ctx->V, ctx->Y[i], ctx->P);
                //_V = _Y[i].pointwise(_P);
                _V = _Y[i].array() * _P.array();
                //VecDot(ctx->V, ctx->S[i], &yDs);
                //yDs = _V * _S[i];
                yDs = _V.dot(_S[i]);
                ys_sum += yDs;

                //VecDot(ctx->V, ctx->V, &yDy);
                //yDy = _V * _V;
                yDy = _V.dot(_V);
                yy_sum += yDy;

                //ss_sum += ctx->ss_rhistory[i];
                ss_sum += _ss_rhistory[i];
              }
            }
            else 
            {
              //VecCopy(ctx->Q, ctx->P);
              _P = _Q;

              //VecPow(ctx->P, ctx->r_beta);
              //_P.pow(_r_beta);
              _P = _P.array().pow(_r_beta);
              //VecPointwiseDivide(ctx->Q, ctx->P, ctx->Q);
              //_Q = _P.pointwiseDivide(_Q);
              _Q = _P.array() / _Q.array();

              //  Compute summations for scalar scaling
              yy_sum = 0; //  No safeguard required
              ys_sum = 0; //  No safeguard required
              ss_sum = 0; //  No safeguard required
              //for (i = 0; i < PetscMin(ctx->nupdates, ctx->rescale_history); ++i) 
              for (i = 0; i < min(_nupdates, _rescale_history); ++i) 
              {
                //VecPointwiseMult(ctx->V, ctx->P, ctx->Y[i]);
                //_V = _P.pointwise(_Y[i]);
                _V = _P.array() * _Y[i].array();
                //VecPointwiseMult(ctx->W, ctx->Q, ctx->S[i]);
                //_W = _Q.pointwise(_S[i]);
                _W = _Q.array(), _S[i].array();

                //VecDot(ctx->V, ctx->V, &yDy);
                yDy = _V.dot(_V);
                //VecDot(ctx->V, ctx->W, &yDs);
                yDs = _V.dot(_W);
                //VecDot(ctx->W, ctx->W, &sDs);
                sDs = _W.dot(_W);

                yy_sum += yDy;
                ys_sum += yDs;
                ss_sum += sDs;
              }
            }

            //if (0.0 == ctx->r_alpha) 
            if (0.0 == _r_alpha) 
            {
              //  Safeguard ys_sum
              if (0.0 == ys_sum) 
              {
                ys_sum = ZERO_SAFEGUARD;
              }
              sigmanew = ss_sum / ys_sum;
            }
            //else if (1.0 == ctx->r_alpha) 
            else if (1.0 == _r_alpha) 
            {
              //  Safeguard yy_sum  
              if (0.0 == yy_sum) 
              {
                ys_sum = ZERO_SAFEGUARD;
              }
              sigmanew = ys_sum / yy_sum;
            }
            else 
            {
              //denom = 2*ctx->r_alpha*yy_sum;
              denom = 2*_r_alpha*yy_sum;

              //  Safeguard denom
              if (0.0 == denom) 
              {
                denom = ZERO_SAFEGUARD;
              }

              //sigmanew = ((2*ctx->r_alpha-1)*ys_sum +
              //             PetscSqrtScalar((2*ctx->r_alpha-1)*(2*ctx->r_alpha-1)*ys_sum*ys_sum -
              //             4*ctx->r_alpha*(ctx->r_alpha-1)*yy_sum*ss_sum)) / denom;
              sigmanew = ((2*_r_alpha-1)*ys_sum +
                           sqrt((2*_r_alpha-1)*(2*_r_alpha-1)*ys_sum*ys_sum -
                           4*_r_alpha*(_r_alpha-1)*yy_sum*ss_sum)) / denom;
            }

            //  If Q has small values, then Q^(r_beta - 1)
            //  can have very large values.  Hence, ys_sum
            //  and ss_sum can be infinity.  In this case,
            //  sigmanew can either be not-a-number or infinity.

            if (isnan(sigmanew) || isinf(sigmanew)) {
              //  sigmanew is not-a-number; skip rescaling
            }
            else if (!sigmanew) {
              //  sigmanew is zero; this is a bad case; skip rescaling
            }
            else {
              //  sigmanew is positive
              //VecScale(ctx->U, sigmanew);
              _U *= sigmanew;
            }
          break;
        }

        //  Modify for previous information
        //switch(ctx->limitType) 
        switch(_limitType) 
        {
          case MatLMVM_Limit_Average:
            //if (1.0 == ctx->mu) 
            if (1.0 == _mu) 
            {
              //VecCopy(ctx->D, ctx->U);
              _U = _D;
            }
            //else if (ctx->mu) 
            else if (_mu) 
            {
              // VecAXPBY(y, alpha, beta, x)
              // y = alpha * x + beta * y
              //VecAXPBY(ctx->D,ctx->mu, 1.0-ctx->mu,ctx->U);
              _D = _mu * _U + (1.0-_mu) * _D;
            }
            break;
 
          case MatLMVM_Limit_Relative:
            //if (ctx->mu) 
            if (_mu) 
            {
              //  P = (1-mu) * D
              //VecAXPBY(ctx->P, 1.0-ctx->mu, 0.0, ctx->D);
              _P = (1.0-_mu) * _D;

              //  Q = (1+mu) * D
              //VecAXPBY(ctx->Q, 1.0+ctx->mu, 0.0, ctx->D);
              _Q = (1.0+_mu) * _D;
              //VecMedian(ctx->P, ctx->U, ctx->Q, ctx->D);
              _D = median(_P, _U, _Q);
            }
            break;

          case MatLMVM_Limit_Absolute:
            //if (ctx->nu) 
            if (_nu) 
            {
              //VecCopy(ctx->P, ctx->D);
              _D = _P; 
              
              //VecShift(ctx->P, -ctx->nu);
              //_P += (-_nu);
              _P = _P.array() - _nu;
              
              //VecCopy(ctx->D, ctx->Q);
              _Q = _D;
              
              //VecShift(ctx->Q, ctx->nu);
              //_Q += _nu;
              _Q = _Q.array() + _nu;
              
              //VecMedian(ctx->P, ctx->U, ctx->Q, ctx->P);
              _P = median(_U, _Q, _P);
            }
            break;

          default:
            //VecCopy(ctx->U, ctx->D);
            _D = _U;
            break;
        } 
        break;
      }
      //PetscObjectDereference((PetscObject)ctx->Xprev);
      //PetscObjectDereference((PetscObject)ctx->Gprev);
      assert(_lm < (int)_S.size());
      _Xprev = _S[_lm]; 
      _Gprev = _Y[_lm];
      //PetscObjectReference((PetscObject)ctx->S[ctx->lm]);
      //iPetscObjectReference((PetscObject)ctx->Y[ctx->lm]);

    } 
    else { 
      _nrejects++;
    }
  }
  
  _iter++;
  _Xprev = x;
  _Gprev = g;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void BLMVM::MatSolve(MATRIX& A, VECTOR& b, VECTOR& x) 
{
  Real sq, yq, dd;
  int ll;
  bool scaled;
  //MatLMVMCtx     *shell;
  //PetscErrorCode ierr;

  //PetscFunctionBegin;
  //PetscValidHeaderSpecific(A,MAT_CLASSID,1);
  //PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  //PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  //ierr = MatShellGetContext(A,(void**)&shell); CHKERRQ(ierr);
  //if (shell->lmnow < 1)
  if (_lmnow < 1) 
  {
    //shell->rho[0] = 1.0;
    _rho[0] = 1.0;
  }

  //VecCopy(b,x);
  x = b;
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " x 1:" << x << endl;
  //for (ll = 0; ll < shell->lmnow; ++ll) 
  for (ll = 0; ll < _lmnow; ++ll) 
  {
    //VecDot(x,shell->S[ll],&sq);
    assert(ll < (int)_S.size());
    //sq = x * _S[ll];
    sq = x.dot(_S[ll]);
    //shell->beta[ll] = sq * shell->rho[ll];
    _beta[ll] = sq * _rho[ll];
    //VecAXPY(x,-shell->beta[ll],shell->Y[ll]);
    x += -_beta[ll] * _Y[ll];
  }
  //cout << " x 2:" << x << endl;

  scaled = false;
  //if (!scaled && !shell->useDefaultH0 && shell->H0) 
  if (!scaled && !(_useDefaultH0) && H0.rows() > 0) 
  {
    //MatSolve(shell->H0,x,shell->U);
    VECTOR backup = x;
    //H0.solve(x);
    x = H0.colPivHouseholderQr().solve(x);
    _U = x;
    x = backup;
   
    //VecDot(x,shell->U,&dd);
    //dd = x * _U;
    dd = x.dot(_U);
    if ((dd > 0.0) && !(isnan(dd) || isinf(dd))) 
    {
        //  Accept Hessian solve
        //VecCopy(shell->U,x);
        x = _U;
        //scaled = PETSC_TRUE;
        scaled = true;
    }
  }
  //cout << " x 3:" << x << endl;

  //if (!scaled && shell->useScale) 
  if (!scaled && _useScale) 
  {
    //VecPointwiseMult(shell->U,x,shell->scale);
    //_U = x.pointwise(_scale);
    _U = x.array() * _scale.array();
    //VecDot(x,shell->U,&dd);
    //dd = x * _U;
    dd = x.dot(_U);
    //if ((dd > 0.0) && !PetscIsInfOrNanReal(dd)) 
    if ((dd > 0.0) && !(isnan(dd) || isinf(dd))) 
    {
        //  Accept scaling
        //VecCopy(shell->U,x);
        x = _U;
        //scaled = PETSC_TRUE;
        scaled = true;
    }
  }
  //cout << " x 4:" << x << endl;
  //cout << " scaled:"  << scaled << endl;
  //cout << " scale type: " << _scaleType << endl;
  //cout << " _D: " << _D << endl;

  if (!scaled) {
    //switch(shell->scaleType) 
    switch(_scaleType) 
    {
      case MatLMVM_Scale_None:
      break;

      case MatLMVM_Scale_Scalar:
      //VecScale(x,shell->sigma);
      x *= _sigma;
      break;
    
      case MatLMVM_Scale_Broyden:
      //VecPointwiseMult(x,x,shell->D);
      //x = x.pointwise(_D);
      x = x.array() * _D.array();
      //cout << " D: " << _D << endl;
      //cout << " x pointwise:" << x << endl;
      break;
    }
  } 
  //cout << " x 5:" << x << endl;

  //for (ll = shell->lmnow-1; ll >= 0; --ll) 
  for (ll = _lmnow-1; ll >= 0; --ll) 
  {
    //VecDot(x,shell->Y[ll],&yq);
    //yq = x * _Y[ll];
    yq = x.dot(_Y[ll]);
    //VecAXPY(x,shell->beta[ll]-yq*shell->rho[ll],shell->S[ll]);
    assert(ll < (int)_S.size());
    x += (_beta[ll]-yq*_rho[ll]) * _S[ll];
  }
  //cout << " x 6:" << x << endl;
  //PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void BLMVM::VecBoundGradientProjection(const VECTOR& G, const VECTOR& X, const VECTOR& XL, const VECTOR& XU, VECTOR& GP)
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
////////////////////////////////////////////////////////////////////////////////////////
void BLMVM::SetDelta(Real& delta)
{
  _delta = fabs(delta);
  _delta = (_delta > _delta_min) ? _delta : _delta_min;
  _delta = (_delta < _delta_max) ? _delta : _delta_max;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
double BLMVM::mid(const double& a, const double& b, const double& c)
{
  vector<Real> sortable(3);
  sortable[0] = a;
  sortable[1] = b;
  sortable[2] = c;
  sort(sortable.begin(), sortable.end());

  return sortable[1];
}

////////////////////////////////////////////////////////////////////////////////////////
// return a vector with the median values from all three vectors
////////////////////////////////////////////////////////////////////////////////////////
VECTOR BLMVM::median(const VECTOR& v0, const VECTOR& v1, const VECTOR& v2)
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
