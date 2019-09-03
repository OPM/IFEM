// $Id$
//==============================================================================
//!
//! \file DiagMatrix.C
//!
//! \date Aug 29 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Diagonal system matrix representation.
//!
//==============================================================================

#include "DiagMatrix.h"
#include "SAM.h"


DiagMatrix::DiagMatrix (const RealArray& data, size_t nrows)
{
  if (nrows == 0) nrows = data.size();

  myMat.resize(nrows);
  memcpy(myMat.ptr(),&data.front(),nrows*sizeof(Real));
}


void DiagMatrix::initAssembly (const SAM& sam, bool)
{
  myMat.resize(sam.neq,true);
}


void DiagMatrix::dump (std::ostream& os, char format, const char* label)
{
  switch (format)
    {
    case 'M':
    case 'm':
      utl::writeMatlab(label,myMat,os);
      break;

    default:
      if (label) os << label <<" =";
      os << myMat;
    }
}


bool DiagMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  if (myMat.size() != (size_t)sam.neq)
    return false;

  std::vector<int> meen;
  if (!sam.getElmEqns(meen,e))
    return false;
  else if (meen.size() != 1 || eM.size() != 1)
  {
    std::cerr <<" *** DiagMatrix::assemble: Only for one-dof elements, nedof = "
              << meen.size() <<" size(eM) = " << eM.size() << std::endl;
    return false;
  }

  int ieq = meen.front();
  if (ieq < 1 || ieq > sam.neq)
  {
    std::cerr <<" *** DiagMatrix::assemble: ieq="<< ieq <<" is out or range [1,"
              << sam.neq <<"]"<< std::endl;
    return false;
  }

  myMat(ieq) += eM(1,1);
  return true;
}


bool DiagMatrix::assemble (const Matrix& eM, const SAM& sam,
                           SystemVector&, int e)
{
  return this->assemble(eM,sam,e);
}


bool DiagMatrix::redim (size_t r)
{
  if (r == myMat.size())
    return false;

  myMat.std::vector<Real>::resize(r);
  return true;
}


bool DiagMatrix::add (const SystemMatrix& B, Real alpha)
{
  const DiagMatrix* Bptr = dynamic_cast<const DiagMatrix*>(&B);
  if (!Bptr || myMat.size() != Bptr->myMat.size())
    return false;

  myMat.add(Bptr->myMat,alpha);
  return true;
}


bool DiagMatrix::add (Real sigma)
{
  for (Real& v : myMat)
    v += sigma;
  return true;
}


bool DiagMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const Real* a = myMat.ptr();
  const Real* b = B.getRef();
  Real*       c = C.getPtr();

  for (size_t i = 0; i < C.size(); i++)
    c[i] = i < B.size() && i < myMat.size() ? a[i]*b[i] : Real(0);

  return true;
}


bool DiagMatrix::solve (SystemVector& B, bool, Real*)
{
  if (myMat.empty()) return true; // Nothing to solve

  int nzero = 0;
  Real* b = B.getPtr();
  for (Real pivot : myMat)
    if (fabs(pivot) < 1.0e-16)
      ++nzero;
    else
      *(b++) /= pivot;

  if (nzero == 0) return true;
  std::cerr <<" *** DiagMatrix::solve: Singular matrix, "
            << nzero <<" pivot elements are zero."<< std::endl;
  return false;
}
