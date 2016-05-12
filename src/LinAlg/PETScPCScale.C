// $Id$
//==============================================================================
//!
//! \file PCScale.C
//!
//! \date Jan 15 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Scaling for PETSc preconditioners
//!
//==============================================================================


#include "PETScPCScale.h"
#include <iostream>


PetscErrorCode PCScaleCreate(PCScale **pcscale)
{
  PCScale *newctx;

#if PETSC_VERSION_MINOR < 5
  PetscNew(PCScale,&newctx);
#else
  PetscNew(&newctx);
#endif
  newctx->scaling  = nullptr;
  *pcscale = newctx;

  return 0;
}


PetscErrorCode PCScaleSetUp(PC pc, Vec *s)
{
  PCScale *shell;

  PCShellGetContext(pc,(void**) &shell);
  shell->scaling = s;

  return 0;
}


PetscErrorCode PCScaleApply(PC pc, Vec x, Vec y)
{
  PCScale *shell;

  PCShellGetContext(pc,(void**)&shell);

  PetscInt n;
  PetscScalar *svec, *xvec, *yvec;

  VecGetLocalSize(*(shell->scaling),&n);
  VecGetArray(*(shell->scaling),&svec);
  VecGetArray(x,&xvec);
  VecGetArray(y,&yvec);

  for (int i = 0;i < n;i++)
    yvec[i] = xvec[i]*svec[i];

  VecRestoreArray(*(shell->scaling),&svec);
  VecRestoreArray(x,&xvec);
  VecRestoreArray(y,&yvec);

  return 0;
}


PetscErrorCode PCScaleDestroy(PC pc)
{
  PCScale *shell;

  PCShellGetContext(pc,(void**)&shell);
  if (shell->scaling)
    VecDestroy(shell->scaling);

  return 0;
}
