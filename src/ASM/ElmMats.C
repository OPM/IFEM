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
  for (std::vector<Matrix>::iterator ait = A.begin(); ait != A.end(); ++ait)
    ait->resize(ndim,ndim);

  for (std::vector<Vector>::iterator bit = b.begin(); bit != b.end(); ++bit)
    bit->resize(ndim);
}


const Matrix& ElmMats::getNewtonMatrix () const
{
  if (A.empty())
  {
    std::cerr <<" *** ElMats::getNewtonMatrix: No element matrices"<< std::endl;
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
    std::cerr <<" *** ElMats::getRHSVector: No element vectors"<< std::endl;
    static Vector empty;
    return empty;
  }
#if SP_DEBUG > 2
  std::cout <<"\nElement right-hand-side vector"<< b.front();
#endif
  return b.front();
}
