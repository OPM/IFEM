// $Id$
//==============================================================================
//!
//! \file PCProd.C
//!
//! \date Jan 15 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Tensor-product of PETSc preconditioners
//!
//==============================================================================


#include "PETScPCProd.h"
#include <iostream>


PetscErrorCode PCProdCreate(PCProd **pcprod)
{
  PCProd *newctx;

#if PETSC_VERSION_MINOR < 5
  PetscNew(PCProd,&newctx);
#else
  PetscNew(&newctx);
#endif
  newctx->pc1 = nullptr;
  newctx->pc2 = nullptr;
  newctx->pc3 = nullptr;
  *pcprod = newctx;

  return 0;
}


PetscErrorCode PCProdSetUp(PC pc, Vec *prec1, PC *prec2, PC *prec3)
{
  PCProd *shell;

  PCShellGetContext(pc,(void**) &shell);
  shell->pc1 = prec1;
  shell->pc2 = prec2;
  shell->pc3 = prec3;
  PCSetUp(*prec2);
  PCSetUp(*prec3);

  return 0;
}


PetscErrorCode PCProdApply(PC pc, Vec x, Vec y)
{
  PCProd *shell;

  PCShellGetContext(pc,(void**)&shell);
  if (shell->pc1) {
    PetscInt n;
    PetscScalar *svec, *xvec, *yvec;
    VecGetArray(*(shell->pc1),&svec);
    VecGetArray(x,&xvec);
    VecGetArray(y,&yvec);
    VecGetLocalSize(*(shell->pc1),&n);

    for (int i = 0;i < n;i++)
      yvec[i] = xvec[i]*svec[i];

    VecRestoreArray(*(shell->pc1),&svec);
    VecRestoreArray(x,&xvec);
    VecRestoreArray(y,&yvec);

    VecCopy(y,x);
  }
  if (shell->pc2) {
    PCApply(*(shell->pc2),x,y);
    VecCopy(y,x);
  }
  if (shell->pc3)
     PCApply(*(shell->pc3),x,y);

  return 0;
}


PetscErrorCode PCProdDestroy(PC pc)
{
  PCProd *shell;

  PCShellGetContext(pc,(void**)&shell);
  if (shell->pc1)
    VecDestroy(shell->pc1);
  if (shell->pc2)
    PCDestroy(shell->pc2);
 if (shell->pc3)
    PCDestroy(shell->pc3);
  PetscFree(shell);

  return 0;
}
