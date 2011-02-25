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
#endif


SystemVector* SystemVector::create (Type vectorType)
{
  switch (vectorType)
    {
    case STD   : return new StdVector();
#ifdef HAS_PETSC
    case PETSC : return new PETScVector();
#endif
    default:
      std::cerr <<"SystemVector::create: Unsupported vector type "
		<< vectorType << std::endl;
    }

  return 0;
}


SystemMatrix* SystemMatrix::create (Type matrixType, const LinSolParams& spar)
{
#ifdef HAS_PETSC
  if (matrixType == PETSC) return new PETScMatrix(spar);
#endif
  return SystemMatrix::create(matrixType);
}


SystemMatrix* SystemMatrix::create (Type matrixType, int num_thread_SLU)
{
  switch (matrixType)
    {
    case DENSE : return new DenseMatrix();
    case SPR   : return new SPRMatrix();
    case SPARSE: return new SparseMatrix(SparseMatrix::SUPERLU,num_thread_SLU);
    case SAMG  : return new SparseMatrix(SparseMatrix::S_A_M_G);
#ifdef HAS_PETSC
    case PETSC :
      // Use default PETSc settings when no parameters are provided by user
      static LinSolParams defaultPar;
      return new PETScMatrix(defaultPar);
#endif
    default:
      std::cerr <<"SystemMatrix::create: Unsupported matrix type "
		<< matrixType << std::endl;
    }

  return 0;
}
