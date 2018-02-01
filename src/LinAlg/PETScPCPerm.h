// $Id
//==============================================================================
//!
//! \file PETScPCPerm.C
//!
//! \date Jan 15 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Permutation for PETSc preconditioners
//!
//==============================================================================

#ifndef PCPERM_H_IS_INCLUDED
#define PCPERM_H_IS_INCLUDED

#include "petscksp.h"


//! \brief Define context for user-defined preconditioner
typedef struct {
  IS  *order; //!< The index set describing the ordering
  Mat Aperm; //!< Permutated matrix
  PC  pc; //!< Preconditioner definition
  PetscBool identity; //!< \e true if the permutation is identity matrix
} PCPerm;

//! \brief Create the permutated preconditioner.
//! \param[out] pc Pointer to the preconditioner to create.
extern PetscErrorCode PCPermCreate(PCPerm** pc);

//! \brief Configure the permutated preconditioner.
//! \param[out] pc Preconditioner to create
//! \param perm The index set describing the ordering
//! \param A The matrix to permutate
//! \param type type for preconditioner to create
extern PetscErrorCode PCPermSetUp(PC pc, IS *perm, Mat A, const char* type);

//! \brief Apply the permutated preconditioner.
//! \param pc Preconditioner to apply
//! \param[in] x Vector to apply preconditioner to
//! \param[out] y The evaluate preconditioner value
extern PetscErrorCode PCPermApply(PC pc, Vec x, Vec y);

//! \brief Destroy a permutated preconditioner.
//! \param pc Preconditioner to destroy
extern PetscErrorCode PCPermDestroy(PC pc);

#endif
