#ifndef PCScale_H_IS_INCLUDED
#define PCScale_H_IS_INCLUDED

#include "petscksp.h"

/* Define context for user-defined preconditioner */
typedef struct {
  Vec *scaling;
} PCScale;

/* Declare routines for user-defined preconditioner */
extern PetscErrorCode PCScaleCreate(PCScale** pc);
extern PetscErrorCode PCScaleSetUp(PC pc, Vec *s);
extern PetscErrorCode PCScaleApply(PC pc, Vec x, Vec y);
extern PetscErrorCode PCScaleDestroy(PC pc);

#endif
