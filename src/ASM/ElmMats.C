// $Id$
//==============================================================================
//!
//! \file ElmMats.C
//!
//! \date Aug 23 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a FEM problem.
//!
//==============================================================================

#include "ElmMats.h"


void ElmMats::resize (size_t nA, size_t nB, size_t nC)
{
  A.resize(nA);
  b.resize(nB);
  c.resize(nC);
}


void ElmMats::redim (size_t ndim)
{
  for (Matrix& Amat : A) Amat.resize(ndim,ndim);
  for (Vector& bvec : b) bvec.resize(ndim);
}


const Matrix& ElmMats::getNewtonMatrix () const
{
  if (A.empty())
  {
#ifdef SP_DEBUG
    std::cout <<"  ** ElMats::getNewtonMatrix: No element matrices"<< std::endl;
#endif
    static Matrix empty;
    return empty;
  }

#if SP_DEBUG > 2
  std::cout <<"\nElement coefficient matrix"<< A.front();
#endif
  return A.front();
}


const Vector& ElmMats::getRHSVector () const
{
  if (b.empty())
  {
#ifdef SP_DEBUG
    std::cout <<"  ** ElMats::getRHSVector: No element vectors"<< std::endl;
#endif
    static Vector empty;
    return empty;
  }

#if SP_DEBUG > 2
  std::cout <<"\nElement right-hand-side vector"<< b.front();
#endif
  return b.front();
}
