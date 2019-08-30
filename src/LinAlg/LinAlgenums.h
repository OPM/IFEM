// $Id$
//==============================================================================
//!
//! \file LinAlgenums.h
//!
//! \date Oct 23 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various enums for linear algebra scope.
//!
//==============================================================================

#ifndef _LINALG_ENUMS_H
#define _LINALG_ENUMS_H


namespace LinAlg //! Linear algebra scope
{
  //! \brief The available system matrix formats and associated solvers.
  enum MatrixType
  {
    DENSE   = 0, //!< Dense matrices / LAPack solver
    SPR     = 1, //!< Sparse matrices / SPR solver
    SPARSE  = 2, //!< Sparse matrices / SuperLU solver
    SAMG    = 3, //!< Sparse matrices / SAMG solver
    PETSC   = 4, //!< Sparse matrices / PETSc solver
    ISTL    = 5, //!< Sparse matrices / Dune solver
    UMFPACK = 6, //!< Sparse matrices / UmfPack solver
    DIAG    = 7  //!< Diagonal matrices / Trivial solver
  };

  //! \brief Enum defining linear system properties.
  enum LinearSystemType
  {
    GENERAL_MATRIX = 0, //!< General matrix
    SYMMETRIC      = 1, //!< Symmetric matrix (value and structure)
    SPD            = 2  //!< Symmetric, positive definite matrix
  };
}

#endif
