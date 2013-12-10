//==============================================================================
//!
//! \file BDFMats.C
//!
//! \date Nov 2 2013
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem
//! using backward difference formulae (BDF)
//!
//==============================================================================

#include "BDFMats.h"

const Matrix& BDFMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A.front());

  double dtcoeff = h;
  if (bdf.getDegree() == 2)
    dtcoeff *= h;

  N = A[1];
  N.multiply(bdf[0]/dtcoeff);    // [N]  = (bdf[0]/dt*dt)*[M]

  N.add(A[2]); // [N] += [K]

#if SP_DEBUG > 2
  std::cout <<"\nElement mass matrix"<< A[1];
  std::cout <<"Element stiffness matrix"<< A[2];
  std::cout <<"Resulting Newton matrix"<< A[0];
#endif

  return A.front();
}


const Vector& BDFMats::getRHSVector () const
{
  if (!A.empty() && vec.size() > 3)
  {
    Vector& db = const_cast<Vector&>(b.front());

    double dtcoeff = h;
    if (bdf.getDegree() == 2)
      dtcoeff *= h;

    const size_t ncoeff = bdf.getCoefs().size();
    for (size_t i = 0;i < ncoeff;i++) 
      db.add(A[1]*vec[i],-bdf[i]/dtcoeff);
  }

#if SP_DEBUG > 2
  std::cout <<"\nElement right-hand-side vector"<< b.front();
#endif

  return b.front();
}
