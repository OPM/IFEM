#ifdef HAS_PETSC
#include "PCPerm.h"
#include <iostream>

PetscErrorCode PCPermCreate(PCPerm **pcperm)
{
  PCPerm *newctx;

  PetscNew(PCPerm,&newctx);
  newctx->order = NULL;
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
  PCSetOperators(shell->pc,shell->Aperm,shell->Aperm,SAME_PRECONDITIONER);
  if (!shell->identity)
    PCFactorSetUseInPlace(shell->pc);
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
  if (shell->order)
    ISDestroy(shell->order);
  if (shell->pc)
    PCDestroy(&(shell->pc));
  PetscFree(shell);

  return 0;
}

#endif
