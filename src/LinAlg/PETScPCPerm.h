#ifndef PCPERM_H_IS_INCLUDED
#define PCPERM_H_IS_INCLUDED

#include "petscksp.h"


/* Define context for user-defined preconditioner */
typedef struct {
  IS  *order;
  Mat Aperm;
  PC  pc;
  PetscBool identity;
} PCPerm;

/* Declare routines for user-defined preconditioner */
extern PetscErrorCode PCPermCreate(PCPerm** pc);
extern PetscErrorCode PCPermSetUp(PC pc, IS *perm, Mat A, const char* type);
extern PetscErrorCode PCPermApply(PC pc, Vec x, Vec y);
extern PetscErrorCode PCPermDestroy(PC pc);

#endif
