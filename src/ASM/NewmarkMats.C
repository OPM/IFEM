// $Id$
//==============================================================================
//!
//! \file NewmarkMats.C
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//! \details See eq. (6.50) and (6.52) in Cottrell et. al. (2009).
//!
//==============================================================================

#include "NewmarkMats.h"


NewmarkMats::NewmarkMats (double a1, double a2, double b, double c)
{
  alpha1 = a1;
  alpha2 = a2;
  beta   = b;
  gamma  = c;

  isPredictor = true;
  h = 0.0;
}


const Matrix& NewmarkMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A.front());

  N = A[1];
  if (alpha1 > 0.0)
    N.multiply(1.0 + alpha1*gamma*h);    // [N]  = (1+alpha1*gamma*h)*[M]

  N.add(A[2],(alpha2*gamma + beta*h)*h); // [N] += (alpha2*gamma+beta*h)*h*[K]

  return A.front();
}


const Vector& NewmarkMats::getRHSVector () const
{
  if (!A.empty() && vec.size() > 2)
  {
    Vector& dF = const_cast<Vector&>(b.front());

    int ia = vec.size() - 1; // index to element acceleration vector (a)
    int iv = vec.size() - 2; // index to element velocity vector (v)

    dF.add(A[1]*vec[ia],-1.0);      // dF = Fext - M*a

    if (alpha1 > 0.0)
      dF.add(A[1]*vec[iv],-alpha1); // dF -= alpha1*M*v

    if (alpha2 > 0.0)
      dF.add(A[2]*vec[iv],-alpha2); // dF -= alpha2*K*v

    dF.add(A[2]*vec.front(),-1.0);  // dF -= K*d
  }

  return b.front();
}
