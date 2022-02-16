// $Id$
//==============================================================================
//!
//! \file BlockElmMats.C
//!
//! \date Jul 10 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the block element matrices for a FEM problem.
//!
//==============================================================================

#include "BlockElmMats.h"


bool BlockElmMats::redim (size_t blkIndex, size_t nen, size_t ncmp, char basis)
{
  bool ok = true;
  if (basis < 0) // Zero diagonal block
    basis = -basis;
  else if (blkIndex <= blockInfo.size() && blkIndex < A.size())
    A[blkIndex].resize(nen*ncmp,nen*ncmp);
  else
    ok = false;

  if (blkIndex <= blockInfo.size() && blkIndex < b.size())
    b[blkIndex].resize(nen*ncmp);
  else
    ok = false;

  if (blkIndex == 0)
    return ok;
  else if (blkIndex <= blockInfo.size())
    blockInfo[blkIndex-1] = Block(ncmp,basis);
  else
    ok = false;

  if (basis > 0 && (size_t)basis <= basisInfo.size())
  {
    basisInfo[basis-1].nen = nen;
    basisInfo[basis-1].ncmp += ncmp;
  }
  else
    ok = false;

  return ok;
}


bool BlockElmMats::redimOffDiag (size_t blkIndex, char symmetric)
{
  size_t nDiagB = blockInfo.size();
  if (blkIndex <= nDiagB || blkIndex > A.size())
    return false; // Not an off-diagonal sub-matrix

  size_t kb = nDiagB;
  bool symm = A.size()-1 < nDiagB*(nDiagB+1)/2;
  if (symmFlag.empty())
    symmFlag.resize(nDiagB*(nDiagB-1) / (symm ? 2 : 1), 0);
  for (size_t ib = 0; ib < nDiagB; ib++)
    for (size_t jb = symm ? ib+1 : 0; jb < nDiagB; jb++)
      if (jb != ib && ++kb == blkIndex)
      {
        A[kb].resize(basisInfo[blockInfo[ib].basis-1].nen*blockInfo[ib].ncmp,
                     basisInfo[blockInfo[jb].basis-1].nen*blockInfo[jb].ncmp);
        kb -= nDiagB;
        symmFlag[kb-1] = symmetric;
        return true;
      }

  return false;
}


bool BlockElmMats::finalize ()
{
  // Calculate the total number of element dofs
  neldof = 0;
  size_t nDiagB = blockInfo.size();
  for (size_t i = 0; i < nDiagB; i++)
    neldof += blockInfo[i].ncmp*basisInfo[blockInfo[i].basis-1].nen;

  // Calculate the offset of each block sub-matrix
  size_t idof = 1;
  for (size_t j = 0; j < basisInfo.size(); j++)
  {
    size_t ndof = 0;
    size_t jdof = idof;
    for (size_t i = 0; i < nDiagB; i++)
      if (blockInfo[i].basis == (int)j+1)
      {
        blockInfo[i].idof = jdof;
        jdof += blockInfo[i].ncmp;
        ndof += blockInfo[i].ncmp*basisInfo[j].nen;
      }
    idof += ndof;
  }

  return idof-1 == neldof;
}


const Matrix& BlockElmMats::getNewtonMatrix () const
{
  if (!A.front().empty())
    return A.front();

  Matrix& N = const_cast<Matrix&>(A.front());
  // Set the element Newton matrix dimension
  N.resize(neldof,neldof);

  size_t ib, jb, kb;
  size_t nDiagB = blockInfo.size();
  bool symmetry = A.size()-1 <= nDiagB*(nDiagB+1)/2;
  for (ib = 0, kb = nDiagB-1; ib < nDiagB; ib++)
  {
    size_t ibas = blockInfo[ib].basis-1;
    size_t neni = basisInfo[ibas].nen;
    size_t ndi  = basisInfo[ibas].ncmp;
    size_t nci  = blockInfo[ib].ncmp;
    size_t idof = blockInfo[ib].idof;

    // Insert the diagonal block matrix
    const Matrix& Aii = A[1+ib];
    if (!Aii.empty())
      for (size_t in = 0; in < neni; in++)
        for (size_t jn = 0; jn < neni; jn++)
          for (size_t i = 0; i < nci; i++)
            for (size_t jd = 0; jd < nci; jd++)
              N(idof+ndi*in+i,idof+ndi*jn+jd) = Aii(1+nci*in+i,1+nci*jn+jd);

    for (jb = symmetry ? ib+1 : 0; jb < nDiagB && 2+kb < A.size(); jb++)
    {
      if (jb == ib)
        continue;
      else
        ++kb;

      size_t jbas = blockInfo[jb].basis-1;
      size_t nenj = basisInfo[jbas].nen;
      size_t ndj  = basisInfo[jbas].ncmp;
      size_t ncj  = blockInfo[jb].ncmp;
      size_t jdof = blockInfo[jb].idof;

      // Insert the off-diagonal block matrix
      const Matrix& Aij = A[1+kb];
      if (!Aij.empty())
        for (size_t in = 0; in < neni; in++)
          for (size_t jn = 0; jn < nenj; jn++)
            for (size_t i = 0; i < nci; i++)
              for (size_t j = 0; j < ncj; j++)
              {
                N(idof+ndi*in+i,jdof+ndj*jn+j) = Aij(1+nci*in+i,1+ncj*jn+j);
                if (symmFlag[kb-nDiagB] > 0)
                  N(jdof+ndj*jn+j,idof+ndi*in+i) = Aij(1+nci*in+i,1+ncj*jn+j);
                else if (symmFlag[kb-nDiagB] < 0)
                  N(jdof+ndj*jn+j,idof+ndi*in+i) = -Aij(1+nci*in+i,1+ncj*jn+j);
              }
    }
  }

#if SP_DEBUG > 2
  std::cout << std::endl;
  for (ib = 1; ib < A.size(); ib++)
    std::cout <<"Block matrix "<< ib << A[ib];
  std::cout <<"Resulting Newton matrix"<< A.front();
#endif

  return A.front();
}


const Vector& BlockElmMats::getRHSVector () const
{
  if (!b.front().empty())
    return b.front();

  Vector& R = const_cast<Vector&>(b.front());
  R.resize(neldof);

  for (size_t ib = 0; ib < blockInfo.size(); ib++)
  {
    size_t ibas = blockInfo[ib].basis-1;
    size_t neni = basisInfo[ibas].nen;
    size_t ndi  = basisInfo[ibas].ncmp;
    size_t nci  = blockInfo[ib].ncmp;
    size_t idof = blockInfo[ib].idof;

    const Vector& bi = b[1+ib];
    for (size_t in = 0; in < neni; in++)
      for (size_t i = 0; i < nci; i++)
        R(idof+ndi*in+i) = bi(1+nci*in+i);
  }

#if SP_DEBUG > 2
  std::cout << std::endl;
  for (size_t i = 1; i < b.size(); i++)
    std::cout <<"Block vector "<< i << b[i];
  std::cout <<"Resulting right-hand-side vector"<< b.front();
#endif

  return b.front();
}
