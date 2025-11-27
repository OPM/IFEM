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
  memcpy(myMat.ptr(),data.data(),nrows*sizeof(Real));
}


void DiagMatrix::initAssembly (const SAM& sam, char)
{
  myMat.resize(sam.neq,true);
}


void DiagMatrix::dump (std::ostream& os, LinAlg::StorageFormat format,
                       const char* label)
{
  switch (format)
  {
    case LinAlg::MATLAB:
      utl::writeMatlab(label,myMat,os);
      break;

    case LinAlg::MATRIX_MARKET:
      os <<"%%MatrixMarket matrix coordinate real general";
      if (label)
        os <<"\n% label = "<< label;
      os <<'\n'<< myMat.size() <<' '<< myMat.size() <<' '<< myMat.size();
      for (size_t row = 1; row <= myMat.size(); row++)
        os <<'\n'<< row <<' '<< row <<' '<< myMat(row);
      os << std::endl;
      break;

    case LinAlg::FLAT:
      if (label) os << label <<" =";
      os << myMat;
      break;

    default:
      break;
  }
}


bool DiagMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  if (static_cast<int>(myMat.size()) != sam.neq)
    return false;

  StdVector B;
  IntVec meen;
  return sam.getElmEqns(meen,e,eM.rows()) && this->assemble(eM,sam,B,meen);
}


bool DiagMatrix::assemble (const Matrix& eM, const SAM& sam,
                           SystemVector&, const IntVec& meq)
{
  const size_t nedof = meq.size();
  if (eM.rows() < nedof || eM.cols() < nedof)
  {
    std::cerr <<" *** DiagMatrix::assemble: Invalid element matrix, nedof = "
              << nedof <<" size(eM) = "<< eM.rows() <<","<< eM.cols()
              << std::endl;
    return false;
  }

  // Add elements corresponding to free dofs in eM,
  // and (appropriately weighted) elements corresponding to constrained dofs,
  // into the diagonal system matrix
  for (size_t j = 1; j <= nedof; j++)
    if (int jeq = meq[j-1]; jeq > 0)
      myMat(jeq) += eM(j,j);
    else if (int jceq = -meq[j-1]; jceq > 0)
      for (int jp = sam.mpmceq[jceq-1]; jp < sam.mpmceq[jceq]-1; jp++)
        if (int jeq = sam.mmceq[jp] > 0 ? sam.meqn[sam.mmceq[jp]-1]:0; jeq > 0)
          for (size_t i = 1; i <= nedof; i++)
            if (meq[i-1] == jeq)
              myMat(jeq) += sam.ttcc[jp]*(eM(i,j) + eM(j,i));
            else if (int iceq = -meq[i-1]; iceq > 0)
              for (int ip = sam.mpmceq[iceq-1]; ip < sam.mpmceq[iceq]-1; ip++)
                if (sam.mmceq[ip] > 0 && sam.meqn[sam.mmceq[ip]-1] == jeq)
                  myMat(jeq) += sam.ttcc[ip]*sam.ttcc[jp]*eM(i,j);

  return this->flagNonZeroEqs(meq);
}


/*!
  This method is doing the same as the above assemble() method,
  but where all element matrices are assumed identity matrices.
  The purpose is to verify that each equation only gets a single contribution
  within each threading group.
*/

bool DiagMatrix::assembleStruct (int val, const SAM& sam, const IntVec& meq)
{
  // Lambda function marking matrix contributions due to linear couplings.
  auto&& addEq = [&diag=myMat,val](int ieq)
  {
    int current = diag(ieq)/1000.0;
    if (current != val)
      diag(ieq) += 1000.0*val;
  };

  const size_t nedof = meq.size();
  for (size_t j = 1; j <= nedof; j++)
    if (int jeq = meq[j-1]; jeq > 0)
      myMat(jeq) += 1.0;
    else if (int jceq = -meq[j-1]; jceq > 0)
      for (int jp = sam.mpmceq[jceq-1]; jp < sam.mpmceq[jceq]-1; jp++)
        if (int jeq = sam.mmceq[jp] > 0 ? sam.meqn[sam.mmceq[jp]-1]:0; jeq > 0)
          for (size_t i = 1; i <= nedof; i++)
            if (meq[i-1] == jeq)
              addEq(jeq);
            else if (int iceq = -meq[i-1]; iceq > 0)
              for (int ip = sam.mpmceq[iceq-1]; ip < sam.mpmceq[iceq]-1; ip++)
                if (sam.mmceq[ip] > 0 && sam.meqn[sam.mmceq[ip]-1] == jeq)
                  addEq(jeq);

  return this->flagNonZeroEqs(meq);
}


bool DiagMatrix::assemble (const Matrix& eM, const SAM& sam,
                           SystemVector&, int e)
{
  return this->assemble(eM,sam,e);
}


bool DiagMatrix::assemble (const Matrix& eM, const IntVec& meq)
{
  for (size_t i = 0; i < meq.size(); ++i)
    (*this)(meq[i]+1) += eM(i+1, i+1);

  return true;
}


bool DiagMatrix::add (const SystemMatrix& B, Real alpha)
{
  const DiagMatrix* Bptr = dynamic_cast<const DiagMatrix*>(&B);
  if (!Bptr || myMat.size() != Bptr->myMat.size())
    return false;

  myMat.add(Bptr->myMat,alpha);

  return this->flagNonZeroEqs(B);
}


bool DiagMatrix::add (Real sigma, int ieq)
{
  if (ieq > static_cast<int>(myMat.size()))
    return false;
  else if (ieq > 0)
    myMat[ieq-1] += sigma;
  else for (Real& v : myMat)
    v += sigma;

  if (ieq == 0)
    return this->flagNonZeroEqs();
  else
    return this->flagNonZeroEqs({ieq});
}


bool DiagMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const Real* a = myMat.ptr();
  const Real* b = B.getRef();
  Real*       c = C.getPtr();

  for (size_t i = 0; i < C.dim(); i++)
    c[i] = i < B.dim() && i < myMat.size() ? a[i]*b[i] : Real(0);

  return true;
}


bool DiagMatrix::solve (SystemVector& B, Real*)
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
