// $Id$
//==============================================================================
//!
//! \file EigSolver.h
//!
//! \date Apr 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Interface to LAPack, ARPack and SLEPc eigenvalue solvers.
//!
//==============================================================================

#ifndef _EIG_SOLVER_H
#define _EIG_SOLVER_H

#include "MatVec.h"

class SystemMatrix;


namespace eig //! Top-level functions for invoking eigenproblem solvers.
{
  //! \brief Solves the eigenvalue problem using LAPACK::DSYGVX.
  //! \details For dense matrices only and no shift.
  //! \param A The system stiffness matrix
  //! \param B The system mass matrix
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors
  //! \param[in] nev Number of eigenvalues/vector to compute
  bool solve(SystemMatrix* A, SystemMatrix* B,
	     Vector& eigVal, Matrix& eigVec, int nev);

  //! \brief Solves the eigenvalue problem (A-lambda*B)*x = 0 using ARPACK.
  //! \param A The system stiffness matrix
  //! \param B The system mass matrix
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors
  //! \param[in] nev Number of eigenvalues/vectors (see ARPack documentation)
  //! \param[in] ncv Number of Arnoldi vectors (see ARPack documentation)
  //! \param[in] mode Eigensolver method (1,...5, see ARPack documentation)
  //! \param[in] shift Eigenvalue shift
  bool solve(SystemMatrix* A, SystemMatrix* B,
	     Vector& eigVal, Matrix& eigVec, int nev, int ncv,
	     int mode = 4, double shift = 0.0);
}

#endif
