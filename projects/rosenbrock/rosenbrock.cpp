//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
#include "TAO_SOLVER.h"
#include "BLMVM.h"
#include "MORE_THUENTE.h"

#include <float.h>
#include <iostream>

using namespace std;

// pointer to the current energy / gradient computation function
void (*functionGradient)(const VECTOR& state, Real& function, VECTOR& gradient);

//////////////////////////////////////////////////////////////////////////////
// The Hessian for the Rosenbrock function
//
// f = (a - x^2) + b(y - x^2)^2
//////////////////////////////////////////////////////////////////////////////
void rosenbrockHessian(const VECTOR& state, MATRIX& hessian)
{
  Real x = state[0];
  Real y = state[1];
  //const Real a = 1.0;
  const Real b = 100.0;

  Real am = 12.0 * b * x * x - 4.0 * b * y + 2.0;
  Real bm = -4.0 * b * x;
  Real cm = -4.0 * b * x;
  Real dm = 2.0 * b;

  Real det = 1.0 / (am * dm - bm * cm); 

  //hessian(0,0) = 12.0 * b * x * x - 4.0 * b * y + 2.0;
  //hessian(0,1) = -4.0 * b * x;
  //hessian(1,0) = -4.0 * b * x;
  //hessian(1,1) = 2.0 * b; 
  hessian(0,0) = dm * det;
  hessian(0,1) = -bm * det;
  hessian(1,0) = -cm * det;
  hessian(1,1) = am * det;
}

//////////////////////////////////////////////////////////////////////////////
// The function gradient for the Rosenbrock function
//
// f = (a - x^2) + b(y - x^2)^2
//////////////////////////////////////////////////////////////////////////////
void rosenbrockFunctionGradient(const VECTOR& state, Real& function, VECTOR& gradient)
{
  cout << ". " << flush;
  Real x = state[0];
  Real y = state[1];
  const Real a = 1.0;
  const Real b = 100.0;

  function = pow(a - x, (Real)2.0) + b * pow((y - x * x), (Real)2.0);

  gradient[0] = -2.0 * a + 4.0 * b * x * x * x - 4.0 * b * x * y + 2.0 * x;
  gradient[1] = 2.0 * b * (y - x * x);
  //cout << " state: " << state << endl;
  cout << function << endl;
  //cout << " grad: " << gradient[0] << " " << gradient[1] << endl;
  //cout << " state: " << x << " " << y << endl;
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  TAO_SOLVER solver;

  VECTOR guess(2);
  //guess[0] = -3;
  //guess[1] = -4;
  
  // for the gradient-only methods
  //guess[0] = -1000;
  //guess[1] = 1000;

  // for the Hessian methods
  // 
  // this is actually really badly conditioned
  //guess[0] = -100;
  //guess[1] = 100;
  guess[0] = -10;
  guess[1] = 10;

  VECTOR upper(2);
  VECTOR lower(2);
  upper = VECTOR::Constant(2, FLT_MAX);
  lower = VECTOR::Constant(2, -FLT_MAX);

  solver.useMoreThuente();
  //solver.useUnitLineSearch();
  //solver.useArmijo();
  //solver.useGPCG();
  //solver.optimizeBLMVM(guess, upper, lower, rosenbrockFunctionGradient);
  //solver.optimizeLMVM(guess, rosenbrockFunctionGradient);
  solver.optimizeNLS(guess, rosenbrockFunctionGradient, rosenbrockHessian);

  return 0;
}
