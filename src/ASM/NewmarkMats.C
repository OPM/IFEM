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

#if SP_DEBUG > 2
  std::cout <<"\nElement mass matrix"<< A[1];
  std::cout <<"Element stiffness matrix"<< A[2];
  std::cout <<"Resulting Newton matrix"<< A[0];
#endif

  return A.front();
}


const Vector& NewmarkMats::getRHSVector () const
{
  if (!A.empty() && vec.size() > 2)
  {
    Vector& dF = const_cast<Vector&>(b.front());

    int ia = vec.size() - 1; // index to element acceleration vector (a)
    int iv = vec.size() - 2; // index to element velocity vector (v)
#if SP_DEBUG > 2
    std::cout <<"\nf_ext"<< dF;
    std::cout <<"f_i = M*a"<< A[1]*vec[ia];
    if (alpha1 > 0.0)
      std::cout <<"f_d1/alpha1 = M*v (alpha1="<< alpha1 <<")"<< A[1]*vec[iv];
    if (alpha2 > 0.0)
      std::cout <<"f_d2/alpha2 = K*v (alpha2="<< alpha2 <<")"<< A[2]*vec[iv];
    std::cout <<"f_s = K*d"<< A[2]*vec.front();
#endif

    dF.add(A[1]*vec[ia],-1.0);      // dF = Fext - M*a

    if (alpha1 > 0.0)
      dF.add(A[1]*vec[iv],-alpha1); // dF -= alpha1*M*v

    if (alpha2 > 0.0)
      dF.add(A[2]*vec[iv],-alpha2); // dF -= alpha2*K*v

    dF.add(A[2]*vec.front(),-1.0);  // dF -= K*d
  }

#if SP_DEBUG > 2
  std::cout <<"\nElement right-hand-side vector"<< b.front();
#endif

  return b.front();
}
