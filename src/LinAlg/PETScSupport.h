// $Id$
//==============================================================================
//!
//! \file PETScSupport.h
//!
//! \date Sep 9 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief IFEM PETSc support
//!
//==============================================================================

#ifndef _PETSC_SUPPORT_H_
#define _PETSC_SUPPORT_H_

#include <vector>

#ifdef HAS_PETSC
#include "petscversion.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscvec.h"
#include "petscis.h"
#include "petscsys.h"

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 4
#include "petscpcmg.h"
#else
#include "petscpc.h"
#endif
#include "petscvec.h"
#define PETSCMANGLE(x) &x
#else
#include "petscmg.h"
#define PETSCMANGLE(x) x
#endif

#ifdef HAS_SLEPC
#include "slepceps.h"
#endif

#else
typedef int    PetscInt;  //!< To avoid compilation failures
typedef double PetscReal; //!< To avoid compilation failures
typedef int    IS;        //!< To avoid compilation failures
#endif

typedef std::vector<PetscInt>    PetscIntVec;  //!< PETSc integer vector
typedef std::vector<PetscReal>   PetscRealVec; //!< PETSc real vector
typedef std::vector<PetscIntVec> PetscIntMat;  //!< PETSc integer matrix

#endif
