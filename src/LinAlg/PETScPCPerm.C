// $Id$
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


#include "PETScPCPerm.h"
#include <iostream>


PetscErrorCode PCPermCreate(PCPerm **pcperm)
{
  PCPerm *newctx;

#if PETSC_VERSION_MINOR < 5
  PetscNew(PCPerm,&newctx);
#else
  PetscNew(&newctx);
#endif
  newctx->order = nullptr;
  *pcperm = newctx;

  return 0;
}


PetscErrorCode PCPermSetUp(PC pc, IS *perm, Mat A, const char* type)
{
  PCPerm *shell;

  PCShellGetContext(pc,(void**) &shell);
  shell->order = perm;
  ISIdentity(*perm,&(shell->identity));
  if (!shell->identity)
    MatPermute(A,*perm,*perm,&(shell->Aperm));
  else
    shell->Aperm = A;
  PCCreate(PETSC_COMM_WORLD,&(shell->pc));
  PCSetType(shell->pc,type);
#if PETSC_VERSION_MINOR < 5
  PCSetOperators(shell->pc,shell->Aperm,shell->Aperm,SAME_PRECONDITIONER);
#else
  PCSetOperators(shell->pc,shell->Aperm,shell->Aperm);
#endif
  if (!shell->identity)
#if PETSC_VERSION_MINOR < 6
    PCFactorSetUseInPlace(shell->pc);
#else
    PCFactorSetUseInPlace(shell->pc,PETSC_TRUE);
#endif
  PCSetUp(shell->pc);

  return 0;
}


PetscErrorCode PCPermApply(PC pc, Vec x, Vec y)
{

  PCPerm *shell;
  PCShellGetContext(pc,(void**)&shell);
  if (!shell->identity)
    VecPermute(x,*(shell->order),PETSC_FALSE);
  PCApply(shell->pc,x,y);
  if (!shell->identity)
    VecPermute(y,*(shell->order),PETSC_TRUE);

  return 0;
}


PetscErrorCode PCPermDestroy(PC pc)
{
  PCPerm *shell;

  PCShellGetContext(pc,(void**)&shell);
  //if (shell->order)
  //  ISDestroy(shell->order);
  if (shell->pc)
    PCDestroy(&(shell->pc));
  PetscFree(shell);

  return 0;
}
