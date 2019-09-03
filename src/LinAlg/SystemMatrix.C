// $Id$
//==============================================================================
//!
//! \file SystemMatrix.C
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General representation of the system matrix on different forms.
//!
//==============================================================================

#include "SystemMatrix.h"
#include "DenseMatrix.h"
#include "SPRMatrix.h"
#include "SparseMatrix.h"
#include "DiagMatrix.h"
#ifdef HAS_PETSC
#include "PETScMatrix.h"
#endif
#ifdef HAS_ISTL
#include "ISTLMatrix.h"
#endif
#include "LinSolParams.h"


SystemVector* SystemVector::create (const ProcessAdm* adm,
                                    LinAlg::MatrixType vType)
{
  switch (vType)
    {
    case LinAlg::ISTL:
#ifdef HAS_ISTL
      if (adm) return new ISTLVector(*adm);
#endif
      break;

    case LinAlg::PETSC:
#ifdef HAS_PETSC
      if (adm) return new PETScVector(*adm);
#endif
      break;

    default:
      return new StdVector();
    }

  std::cerr <<"SystemVector::create: Unsupported vector type "
            << vType << std::endl;
  return nullptr;
}


SystemVector& SystemVector::copy (const SystemVector& x)
{
  this->redim(x.size());
  Real* vec = this->getPtr();
  memcpy(vec,x.getRef(),x.dim()*sizeof(Real));
  this->restore(vec);

  return *this;
}


void StdVector::dump (std::ostream& os, char format, const char* label)
{
  switch (format)
    {
    case 'M':
    case 'm':
      utl::writeMatlab(label,*this,os);
      break;

    default:
      if (label) os << label <<" =";
      os << *this;
    }
}


SystemMatrix* SystemMatrix::create (const ProcessAdm* adm,
                                    LinAlg::MatrixType mType,
                                    const LinSolParams& spar)
{
#ifdef HAS_PETSC
  if (mType == LinAlg::PETSC && adm)
    return new PETScMatrix(*adm,spar);
#endif
#ifdef HAS_ISTL
  if (mType == LinAlg::ISTL && adm)
    return new ISTLMatrix(*adm,spar);
#endif

  return SystemMatrix::create(adm,mType);
}


SystemMatrix* SystemMatrix::create (const ProcessAdm* adm,
                                    LinAlg::MatrixType mType,
                                    int num_thread_SLU)
{
#if defined(HAS_PETSC) || defined(HAS_ISTL)
  // Use default settings when no parameters are provided by user
  static LinSolParams defaultPar;
#endif

  switch (mType)
    {
    case LinAlg::DENSE:
      return new DenseMatrix();

    case LinAlg::SPR:
      return new SPRMatrix();

    case LinAlg::SPARSE:
      return new SparseMatrix(SparseMatrix::SUPERLU,num_thread_SLU);

    case LinAlg::SAMG:
      return new SparseMatrix(SparseMatrix::S_A_M_G);

    case LinAlg::PETSC:
#ifdef HAS_PETSC
      if (adm) return new PETScMatrix(*adm,defaultPar);
#else
      std::cerr <<"SystemMatrix::create: PETSc not compiled in."<< std::endl;
#endif
      break;

    case LinAlg::ISTL:
#ifdef HAS_ISTL
      if (adm) return new ISTLMatrix(*adm,defaultPar);
#else
      std::cerr <<"SystemMatrix::create: ISTL not compiled in."<< std::endl;
#endif
      break;

    case LinAlg::UMFPACK:
      return new SparseMatrix(SparseMatrix::UMFPACK);

    case LinAlg::DIAG:
      return new DiagMatrix();

    default:
      break;
    }

  std::cerr <<"SystemMatrix::create: Unsupported matrix type "
            << mType << std::endl;
  return nullptr;
}


bool SystemMatrix::assemble (const Matrix&, const SAM&,
                             SystemVector&, const std::vector<int>&)
{
  std::cerr <<"SystemMatrix::assemble(const Matrix&,const SAM&,"
            <<"SystemVector&,const std::vector<int>&): "
            <<"Not implemented for the chosen matrix type."<< std::endl;
  return false;
}


StdVector SystemMatrix::operator* (const SystemVector& b) const
{
  StdVector results;
  this->multiply(b,results);
  return results;
}


StdVector SystemMatrix::operator/ (const SystemVector& b)
{
  StdVector results;
  this->solve(b,results);
  return results;
}
