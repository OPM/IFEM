// $Id$
//==============================================================================
//!
//! \file LinSolParams.h
//!
//! \date Jan 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//! \details Includes linear solver method, preconditioner
//! and convergence criteria.
//!
//==============================================================================

#ifndef _LINSOLPARAMS_MATRIX_H
#define _LINSOLPARAMS_MATRIX_H

#include <iostream>
#ifdef HAS_PETSC
#include <string>
#include "petscksp.h"
#endif


//! \brief Null-space for matrix operator
enum NullSpace { NONE, CONSTANT, RIGID_BODY };

/*!
  \brief Class for linear solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class LinSolParams
{
public:
  //! \brief Default constructor.
  LinSolParams() { this->setDefault(); }
  //! \brief Copy constructor.
  LinSolParams(const LinSolParams& spar) { this->copy(spar); }

  //! \brief Set default values.
  virtual void setDefault();

  //! \brief Copy linear solver parameters.
  virtual void copy(const LinSolParams& spar);

  //! \brief Read linear solver parameters from stream.
  virtual bool read(std::istream& is, int nparams = 10);

  //! \brief Read linear solver parameters from file.
  virtual bool read(const char* filename);

#ifdef HAS_PETSC
  //! \brief Set linear solver method
  virtual void setMethod(KSPType type) { method.assign(type); }

  //! \brief Set preconditioner
  virtual void setPreconditioner(PCType type) { prec.assign(type); }

  //! \brief Set linear solver package
  virtual void setPackage(MatSolverPackage stype) { package.assign(stype); }

  //! \brief Set absolute convergence tolerance
  virtual void setAbsTolerance(real eps) { atol = eps; }

  //! \brief Set relative convergence tolerance
  virtual void setRelTolerance(real eps) { rtol = eps; }

  //! \brief Set divergence tolerance
  virtual void setDivTolerance(real eps) { dtol = eps; }

  //! \brief Set maximum number of iterations
  virtual void setMaxIterations(int its) { maxIts = its; }

  //! \brief Set number of overlap
  virtual void setOverlap(int olap) { overlap = olap; }

  //! \brief Set null-space of matrix
  virtual void setNullSpace(NullSpace nspc) { nullspc = nspc; }

  //! \brief Get linear solver method
  virtual const char* getMethod() const { return method.c_str(); }

  //! \brief Get preconditioner
  virtual const char* getPreconditioner() const { return prec.c_str(); }

  //! \brief Get linear solver package
  virtual const char* getPackage() const { return package.c_str(); }

  //! \brief Get absolute convergence tolerance
  virtual real getAbsTolerance() const { return atol; }

  //! \brief Get relative convergence tolerance
  virtual real getRelTolerance() const { return rtol; }

  //! \brief Get divergence tolerance
  virtual real getDivTolerance() const { return dtol; }

  //! \brief Get maximum number of iterations
  virtual int getMaxIterations() const { return maxIts; }

  //! \brief Get number of overlaps
  virtual int getOverlap() const { return overlap; }

  //! \brief Get number of overlaps
  virtual NullSpace getNullSpace() const { return nullspc; }

  //! \brief Set linear solver parameters for KSP object
  virtual void setParams(KSP& ksp) const;

private:
  std::string method;   // Linear solver method
  std::string prec;     // Preconditioner
  std::string package;  // Linear software package (petsc, superlu_dist, ...)
  PetscReal atol;       // Absolute tolerance
  PetscReal rtol;       // Relative tolerance
  PetscReal dtol;       // Divergence tolerance
  PetscInt  levels;     // Number of levels of fill to use
  PetscInt  maxIts;     // Maximum number of iterations
  PetscInt  overlap;    // Number of overlaps
  NullSpace nullspc;    // Null-space for matrix
#endif // HAS_PETSC
};

#endif
