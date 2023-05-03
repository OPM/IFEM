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
  this->printMat(std::cout);
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
  this->printVec(std::cout);
#endif
  return b.front();
}


void ElmMats::printMat (std::ostream& os, size_t idx) const
{
  if (idx > A.size()) return;

  os <<"\nElement coefficient matrix";
  if (idx > 0) os <<" "<< idx;

  if (idx < Aname.size() && Aname[idx])
    os <<" ("<< Aname[idx] <<")";

  os << A[idx];
}


void ElmMats::printVec (std::ostream& os, size_t idx) const
{
  if (idx > b.size()) return;

  os <<"\nElement right-hand-side vector";
  if (idx > 0) os <<" "<< idx;

  if (idx < Bname.size() && Bname[idx])
    os <<" ("<< Bname[idx] <<")";

  os << b[idx];
}


void ElmMats::printScl (std::ostream& os, size_t idx) const
{
  if (idx > c.size()) return;

  os <<"Element scalar";
  if (idx > 0) os <<" "<< idx;

  if (idx < Cname.size() && Cname[idx])
    os <<" ("<< Cname[idx] <<"): ";

  os << c[idx] << std::endl;
}
