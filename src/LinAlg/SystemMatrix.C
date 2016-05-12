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
#ifdef HAS_PETSC
#include "PETScMatrix.h"
#include "PETScBlockMatrix.h"
#include "LinSolParams.h"
#endif


SystemVector* SystemVector::create (const ProcessAdm& adm, Type vectorType)
{
  switch (vectorType)
    {
    case STD   : return new StdVector();
#ifdef HAS_PETSC
    case PETSC : return new PETScVector(adm);
#endif
    default:
      std::cerr <<"SystemVector::create: Unsupported vector type "
                << vectorType << std::endl;
    }

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


SystemMatrix* SystemMatrix::create (const ProcessAdm& padm, Type matrixType,
                                    const LinSolParams& spar,
                                    LinAlg::LinearSystemType ltype)
{
#ifdef HAS_PETSC
  if (matrixType == PETSC) {
    if (spar.getNoBlocks() > 1)
      return new PETScBlockMatrix(padm,spar);
    else
      return new PETScMatrix(padm,spar,ltype);
  }
  else if (matrixType == PETSCBLOCK)
    return new PETScBlockMatrix(padm,spar);
#endif

  return SystemMatrix::create(padm,matrixType,ltype);
}


SystemMatrix* SystemMatrix::create (const ProcessAdm& padm, Type matrixType,
                                    LinAlg::LinearSystemType ltype,
                                    int num_thread_SLU)
{
#ifndef HAS_PETSC
  if (matrixType == PETSC || matrixType == PETSCBLOCK) {
    std::cerr <<"SystemMatrix::create: PETSc not compiled in, bailing out..."
              << std::endl;
    exit(1);
  }
#else
  // Use default PETSc settings when no parameters are provided by user
  static LinSolParams defaultPar(3);
#endif

  switch (matrixType)
    {
    case DENSE : return new DenseMatrix();
    case SPR   : return new SPRMatrix();
    case SPARSE: return new SparseMatrix(SparseMatrix::SUPERLU,num_thread_SLU);
    case SAMG  : return new SparseMatrix(SparseMatrix::S_A_M_G);
#ifdef HAS_PETSC
    case PETSC :      return new PETScMatrix(padm,defaultPar,ltype);
    case PETSCBLOCK : return new PETScBlockMatrix(padm,defaultPar);
#endif
    default:
      std::cerr <<"SystemMatrix::create: Unsupported matrix type "
                << matrixType << std::endl;
    }

  return 0;
}


bool SystemMatrix::assemble (const Matrix&, const SAM&,
                             SystemVector&, const std::vector<int>&)
{
  std::cerr <<"SystemMatrix::assemble(const Matrix&,const SAM&,"
            <<"SystemVector&,const std::vector<int>&): "
            <<"Not implemented for the chosen matrix type."<< std::endl;
  return false;
}

//! \brief Matrix-vector product
StdVector SystemMatrix::operator*(const StdVector& b) const
{
  StdVector results;
  multiply(b, results);
  return results;
}

//! \brief Solve linear system
StdVector SystemMatrix::operator/(const StdVector& b)
{
  StdVector results;
  solve(b, results);
  return results;
}
