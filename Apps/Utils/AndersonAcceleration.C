// $Id$
//==============================================================================
//!
//! \file AndersonAcceleration.C
//!
//! \date Sep 18 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Anderson acceleration of a fixed point iteration.
//!
//==============================================================================

#include "AndersonAcceleration.h"

#include "SAM.h"

#include <iostream>


void AndersonAcceleration::clear ()
{
  dF.clear();
  dG.clear();
  prevF.clear();
  prevG.clear();
  bestF = -1.0;
}


bool AndersonAcceleration::accelerate (const SAM& sam,
                                       const Vector& x, Vector& gx)
{
  if (depth == 0)
    return true;

  // Residual of the step taken
  Vector f(gx);
  f -= x;

  const double normF = sam.norm2(f,'A');
  if (bestF >= 0.0 && normF > maxGrowth*bestF)
    this->clear(); // The combination is not helping, start over from here

  // Not an else: the history restarts from this step, so this step is what
  // the next one is measured against. Leaving the baseline unset until the
  // step after would let that one grow unchecked
  if (bestF < 0.0 || normF < bestF)
    bestF = normF;

  if (!prevF.empty()) {
    Vector df(f);
    df -= prevF;
    Vector dg(gx);
    dg -= prevG;
    dF.push_back(std::move(df));
    dG.push_back(std::move(dg));
    // A loop rather than a test, so that a depth lowered after the history
    // has been built is reached rather than only approached
    while (dF.size() > depth) {
      dF.erase(dF.begin());
      dG.erase(dG.begin());
    }
  }

  prevF = f;
  prevG = gx;

  const size_t m = dF.size();
  if (m == 0)
    return true; // Nothing to combine with yet

  // Normal equations of the least squares problem, which are of the size of
  // the history and so are formed and solved directly
  Matrix A(m,m);
  RealArray b(m);
  double trace = 0.0;
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = i; j < m; ++j)
      A(1+i,1+j) = A(1+j,1+i) = sam.dot(dF[i],dF[j],'A');
    b[i] = sam.dot(dF[i],f,'A');
    trace += A(1+i,1+i);
  }

  if (trace <= 0.0)
    return true; // The steps have stopped differing

  for (size_t i = 0; i < m; ++i)
    A(1+i,1+i) += eps*trace;

  if (!utl::solve(A,b)) {
    std::cerr <<"  ** AndersonAcceleration: The least squares problem over "
              << m <<" steps could not be solved, taking the step as it is."
              << std::endl;
    this->clear();
    return false;
  }

  // Combine the arrival points
  for (size_t i = 0; i < m; ++i)
    gx.add(dG[i],-b[i]);

  return true;
}
