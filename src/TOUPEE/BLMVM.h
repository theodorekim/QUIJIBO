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

#ifndef _BLMVM_H_
#define _BLMVM_H_

#include "SETTINGS.h"

#include <vector>
using namespace std;

#ifndef ZERO_SAFEGUARD
#define ZERO_SAFEGUARD	1e-8
#endif

#ifndef INF_SAFEGUARD
#define INF_SAFEGUARD	1e+8
#endif

////////////////////////////////////////////////////////////////////////////////////////
// port of the TAO BLMVM implementation
////////////////////////////////////////////////////////////////////////////////////////
class BLMVM 
{
public:
  BLMVM();
  void MatReset(MATRIX& M);
  void MatAllocateVectors(MATRIX& m, VECTOR& v);
  void Update(MATRIX& M, VECTOR& x, VECTOR& g);
  void MatSolve(MATRIX& A, VECTOR& b, VECTOR& x); 
  void VecBoundGradientProjection(const VECTOR& G, const VECTOR& X, const VECTOR& XL, const VECTOR& XU, VECTOR& GP);
  void SetDelta(Real& delta);
  double mid(const double& a, const double& b, const double& c);
  VECTOR median(const VECTOR& v0, const VECTOR& v1, const VECTOR& v2);

public:
  enum SCALE_TYPE {MatLMVM_Scale_None, MatLMVM_Scale_Scalar, MatLMVM_Scale_Broyden, MatLMVM_Scale_Types};
  enum RESCALE_TYPE {MatLMVM_Rescale_None, MatLMVM_Rescale_Scalar, MatLMVM_Rescale_GL, MatLMVM_Rescale_Types};
  enum LIMIT_TYPE {MatLMVM_Limit_None, MatLMVM_Limit_Average,	MatLMVM_Limit_Relative, MatLMVM_Limit_Absolute, MatLMVM_Limit_Types};

  bool _allocated;
  int _lm;
  vector<VECTOR> _S;
  vector<VECTOR> _Y;
  VECTOR _D;
  VECTOR _U;
  VECTOR _V;
  VECTOR _W;
  VECTOR _P;
  VECTOR _Q;
  Real _delta_min;
  Real _delta_max;

  Real _eps;
  LIMIT_TYPE _limitType;
  SCALE_TYPE _scaleType;
  RESCALE_TYPE _rScaleType;
    
  Real _s_alpha;	/*  Factor for scalar scaling */
  Real _r_alpha;	/*  Factor on scalar for rescaling diagonal matrix */
  Real _r_beta;	/*  Factor on diagonal for rescaling diagonal matrix */
  Real _mu;		/*  Factor for using historical information */
  Real _nu;		/*  Factor for using historical information */
  Real _phi;		/*  Factor for Broyden scaling */

  int _scalar_history;	/*  Amount of history to keep for scalar scaling */
  VECTOR _yy_history;	/*  Past information for scalar scaling */
  VECTOR _ys_history;	/*  Past information for scalar scaling */
  VECTOR _ss_history;	/*  Past information for scalar scaling */

  int _rescale_history;  /*  Amount of history to keep for rescaling diagonal */
  VECTOR _yy_rhistory;	/*  Past information for scalar rescaling */
  VECTOR _ys_rhistory;	/*  Past information for scalar rescaling */
  VECTOR _ss_rhistory;	/*  Past information for scalar rescaling */

  int _lmnow;
  int _iter;
  int _nupdates;
  int _nrejects;

  VECTOR _Gprev;
  VECTOR _Xprev;

  Real _delta;
  Real _sigma;
  VECTOR _rho;
  VECTOR _beta;

  bool _useDefaultH0;
  MATRIX H0;

  bool _useScale;
  VECTOR _scale;
};
#endif
