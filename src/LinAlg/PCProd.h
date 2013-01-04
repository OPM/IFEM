#ifndef PCProd_H_IS_INCLUDED
#define PCProd_IS_INCLUDED

#include "petscksp.h"

/* Define context for user-defined preconditioner */
typedef struct {
  PC *pc1;
  PC *pc2;
} PCProd;

/* Declare routines for user-defined preconditioner */
extern PetscErrorCode PCProdCreate(PCProd** pc);
extern PetscErrorCode PCProdSetUp(PC pc, PC *prec1, PC *prec2);
extern PetscErrorCode PCProdApply(PC pc, Vec x, Vec y);
extern PetscErrorCode PCProdDestroy(PC pc);

#endif
