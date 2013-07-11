// $Id$
//==============================================================================
//!
//! \file GenAlphaMats.C
//!
//! \date Jul 04 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//!
//==============================================================================

#include "GenAlphaMats.h"


GenAlphaMats::GenAlphaMats (double alpha, double a, double b) : NewmarkMats(a,b)
{
  beta = 0.25*(1.0-alpha)*(1.0-alpha);
  gamma = 0.5 - alpha;
}


const Matrix& GenAlphaMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A.front());

  N = A[2];
  if (alpha2 > 0.0)
    N.multiply(1.0 + alpha2*gamma/(beta*h));

  N.add(A[1],(gamma*alpha1 + 1.0/h)/(beta*h));

  return A.front();
}


const Vector& GenAlphaMats::getRHSVector () const
{
  int ia = vec.size() - 1; // index to element acceleration vector (a)
  int iv = vec.size() - 2; // index to element velocity vector (v)

  if (A.size() > 2 && b.size() > 1 && vec.size() > 2)
  {
    Vector& Fi = const_cast<Vector&>(b[1]);

    // Find the actual inertia force from the dynamic equilibrium equation
    // (Fi = Fext - Fs - Fd) and store it in the second RHS vector Fi=b[1]
    Fi = b.front(); // Fi = Fext - Fs
    if (alpha1 > 0.0 && iv > 0)
      Fi.add(A[1]*vec[iv],-alpha1); // Fi -= alpha1*M*v
    if (alpha2 > 0.0 && iv > 0)
      Fi.add(A[2]*vec[iv],-alpha2); // Ri -= alpha2*K*v
  }

  if (A.size() > 2 && !b.empty() && vec.size() > 2)
  {
    Vector& RHS = const_cast<Vector&>(b.front());

    // Find the right-hand-side vector
    double alphaPlus1 = 1.5 - gamma;
    RHS *= alphaPlus1;
    if (alpha1 > 0.0)
      RHS.add(A[1]*vec[iv], alpha1*(isPredictor ? -alphaPlus1 : alphaPlus1));
    if (alpha2 > 0.0)
      RHS.add(A[2]*vec[iv], alpha2*(isPredictor ? -alphaPlus1 : alphaPlus1));

    RHS.add(A[1]*vec[ia], isPredictor ? 1.0 : -1.0);
  }

  return b.front();
}
