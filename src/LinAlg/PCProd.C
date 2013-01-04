#include "PCProd.h"
#include <iostream>

PetscErrorCode PCProdCreate(PCProd **pcprod)
{
  PCProd *newctx;

  PetscNew(PCProd,&newctx);
  newctx->pc1  = NULL;
  newctx->pc2 = NULL;
  *pcprod = newctx;

  return 0;
}


PetscErrorCode PCProdSetUp(PC pc,PC *prec1, PC *prec2)
{
  PCProd *shell;
  
  PCShellGetContext(pc,(void**) &shell);
  shell->pc1 = prec1;
  shell->pc2 = prec2;
  PCSetUp(*prec1);
  PCSetUp(*prec2);

  return 0;
}


PetscErrorCode PCProdApply(PC pc, Vec x, Vec y)
{
  PCProd *shell;

  PCShellGetContext(pc,(void**)&shell);
  PCApply(*(shell->pc1),x,y);
  VecCopy(y,x);
  PCApply(*(shell->pc2),x,y);
  
  return 0;
}


PetscErrorCode PCProdDestroy(PC pc)
{
  PCProd *shell;
  
  PCShellGetContext(pc,(void**)&shell);
  if (shell->pc1)
    PCDestroy(shell->pc1);
  if (shell->pc2)
    PCDestroy(shell->pc2);
  PetscFree(shell);

  return 0;
}
