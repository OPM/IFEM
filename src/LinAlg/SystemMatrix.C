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
#include "LinSolParams.h"


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


SystemMatrix* SystemMatrix::create (Type matrixType, const LinSolParams& spar)
{
#ifdef HAS_PETSC
  if (matrixType == PETSC)
    return new PETScMatrix(spar);
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
    case PETSC :
#ifdef HAS_PETSC
      // Use default PETSc settings when no parameters are provided by user
      static LinSolParams defaultPar;
      return new PETScMatrix(defaultPar);
#else
      std::cerr <<"SystemMatrix::create: PETSc not compiled in, bailing out..."
		<< std::endl;
      exit(1);
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
