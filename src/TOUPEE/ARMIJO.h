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

#ifndef _ARMIJO_H_
#define _ARMIJO_H_

#include "LINE_SEARCH.h"

/* Context for an Armijo (nonmonotone) linesearch for unconstrained 
   minimization.
   
   Given a function f, the current iterate x, and a descent direction d:
   Find the smallest i in 0, 1, 2, ..., such that:
  
      f(x + (beta**i)d) <= f(x) + (sigma*beta**i)<grad f(x),d>
  
   The nonmonotone modification of this linesearch replaces the f(x) term
   with a reference value, R, and seeks to find the smallest i such that:
  
      f(x + (beta**i)d) <= R + (sigma*beta**i)<grad f(x),d>
  
   This modification does effect neither the convergence nor rate of 
   convergence of an algorithm when R is chosen appropriately.  Essentially,
   R must decrease on average in some sense.  The benefit of a nonmonotone
   linesearch is that local minimizers can be avoided (by allowing increase
   in function value), and typically, fewer iterations are performed in
   the main code.
  
   The reference value is chosen based upon some historical information
   consisting of function values for previous iterates.  The amount of
   historical information used is determined by the memory size where the
   memory is used to store the previous function values.  The memory is
   initialized to alpha*f(x^0) for some alpha >= 1, with alpha=1 signifying
   that we always force decrease from the initial point.
  
   The reference value can be the maximum value in the memory or can be
   chosen to provide some mean descent.  Elements are removed from the
   memory with a replacement policy that either removes the oldest
   value in the memory (FIFO), or the largest value in the memory (MRU).
  
   Additionally, we can add a watchdog strategy to the search, which
   essentially accepts small directions and only checks the nonmonotonic
   descent criteria every m-steps.  This strategy is NOT implemented in
   the code.
  
   Finally, care must be taken when steepest descent directions are used.
   For example, when the Newton direction is not not satisfy a sufficient
   descent criteria.  The code will apply the same test regardless of
   the direction.  This type of search may not be appropriate for all
   algorithms.  For example, when a gradient direction is used, we may 
   want to revert to the best point found and reset the memory so that
   we stay in an appropriate level set after using a gradient steps.
   This type of search is currently NOT supported by the code.
  
   References:
    Armijo, "Minimization of Functions Having Lipschitz Continuous 
      First-Partial Derivatives," Pacific Journal of Mathematics, volume 16,
      pages 1-3, 1966.
    Ferris and Lucidi, "Nonmonotone Stabilization Methods for Nonlinear
      Equations," Journal of Optimization Theory and Applications, volume 81,
      pages 53-71, 1994.
    Grippo, Lampariello, and Lucidi, "A Nonmonotone Line Search Technique
      for Newton's Method," SIAM Journal on Numerical Analysis, volume 23,
      pages 707-716, 1986.
    Grippo, Lampariello, and Lucidi, "A Class of Nonmonotone Stabilization
      Methods in Unconstrained Optimization," Numerische Mathematik, volume 59,
      pages 779-805, 1991. */
class ARMIJO : public LINE_SEARCH
{
public:
  ARMIJO(void (*functionGradient)(const VECTOR&, Real&, VECTOR&)); 

  ////////////////////////////////////////////////////////////////////////////////////////
  // Originally:
  //    static PetscErrorCode TaoLineSearchApply_Armijo(TaoLineSearch ls,Vec x,PetscReal *f,
  //                                                  Vec g,Vec step_direction)
  // In file: src/linesearch/impls/armijo/armijo.c
  ////////////////////////////////////////////////////////////////////////////////////////
  virtual void LineSearchImplementation(VECTOR& x, Real& f, VECTOR& g, VECTOR& s);

private:
  VECTOR _memory;

  Real _alpha;			    // Initial reference factor >= 1
  Real _beta;			      // Steplength determination < 1
  Real _beta_inf;		    // Steplength determination < 1
  Real _sigma;			    // Acceptance criteria < 1)
  Real _minimumStep;		// Minimum step size
  Real _lastReference;	// Reference value of last iteration

  int _memorySize;		    // Number of functions kept in memory
  int _current;			      // Current element for FIFO
  int _referencePolicy;		// Integer for reference calculation rule
  int _replacementPolicy;	// Policy for replacing values in memory
  
  bool _nondescending;     
  bool _memorySetup;

  VECTOR _x;        // Maintain reference to variable vector to check for changes
  VECTOR _work; 
};

#endif
