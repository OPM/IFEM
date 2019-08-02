// $Id$
//==============================================================================
//!
//! \file PETScSchurPC.C
//!
//! \date Jun 4 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Schur-complement preconditioner using PETSc.
//!
//==============================================================================

#include "PETScSchurPC.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#include <tinyxml.h>


PETScSchurPC::PETScSchurPC (PC& pc_init, const std::vector<Mat>& blocks,
                            const LinSolParams::BlockParams& params, const ProcessAdm& adm)
  : m_blocks(&blocks)
{
  PCSetType(pc_init, PCSHELL);
  PCShellSetContext(pc_init, this);
  PCShellSetName(pc_init, "Schur complement preconditioner");
  PCShellSetApply(pc_init, PETScSchurPC::Apply_Outer);
  PCShellSetDestroy(pc_init, PETScSchurPC::Destroy);

  KSPCreate(*adm.getCommunicator(), &inner_ksp);
  KSPSetType(inner_ksp, KSPPREONLY);
  KSPSetOperators(inner_ksp, blocks[0], blocks[0]);
  PC pc;
  KSPGetPC(inner_ksp, &pc);
  std::string schurpc = params.getStringValue("schurpc");
  if (schurpc.empty())
    schurpc = PCGAMG;
  PCSetType(pc, schurpc.c_str());
  KSPSetFromOptions(inner_ksp);
  KSPSetUp(inner_ksp);
  KSPView(inner_ksp, PETSC_VIEWER_STDOUT_WORLD);

  KSPCreate(*adm.getCommunicator(), &outer_ksp);
  KSPGetPC(outer_ksp, &pc);
  PCSetType(pc, PCNONE);
  int maxits = params.getIntValue("schur_maxits");
  if (maxits < 1)
    maxits = 1000;
  double atol = params.getDoubleValue("schur_atol");
  if (atol == 0.0)
    atol = 1e-16;
  double dtol = params.getDoubleValue("schur_dtol");
  if (dtol == 0.0)
    dtol = 1e100;
  double rtol = params.getDoubleValue("schur_rtol");
  if (rtol == 0.0)
    rtol = 1e-12;
  std::string type = params.getStringValue("schur_type");
  if (type.empty())
    type = KSPGMRES;
  KSPSetType(outer_ksp, type.c_str());
  KSPSetTolerances(outer_ksp, rtol, atol, dtol, maxits);

  MatCreate(*adm.getCommunicator(), &outer_mat);
  PetscInt r;
  MatGetSize(blocks[3], &r, &r);
  MatSetSizes(outer_mat, r, r, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetFromOptions(outer_mat);
  MatSetType(outer_mat, MATSHELL);
  MatShellSetContext(outer_mat, this);
  MatShellSetOperation(outer_mat, MATOP_MULT,
                       (void(*)(void))&PETScSchurPC::Apply_Schur);
  MatSetUp(outer_mat);

  KSPSetOperators(outer_ksp, outer_mat, outer_mat);
  KSPSetFromOptions(outer_ksp);
  KSPSetUp(outer_ksp);
  KSPView(outer_ksp, PETSC_VIEWER_STDOUT_WORLD);
  PCSetUp(pc_init);

  VecCreate(*adm.getCommunicator(), &tmp);
  VecSetFromOptions(tmp);
  MatGetSize(blocks[0], &r, &r);
  VecSetSizes(tmp, r, PETSC_DETERMINE);
  VecSetUp(tmp);
}


PETScSchurPC::~PETScSchurPC ()
{
  KSPDestroy(&inner_ksp);
  KSPDestroy(&outer_ksp);
  MatDestroy(&outer_mat);
  VecDestroy(&tmp);
}


PetscErrorCode PETScSchurPC::Apply_Schur (Mat A, Vec x, Vec y)
{
  void* p;
  MatShellGetContext(A, &p);
  PETScSchurPC* spc = static_cast<PETScSchurPC*>(p);
  MatMult(spc->m_blocks->at(1), x, spc->tmp);
  KSPSolve(spc->inner_ksp, spc->tmp, spc->tmp);
  MatMult(spc->m_blocks->at(2), spc->tmp, y);

  return 0;
}


PetscErrorCode PETScSchurPC::Apply_Outer (PC pc, Vec x, Vec y)
{
  void* p;
  PCShellGetContext(pc, &p);
  PETScSchurPC* spc = static_cast<PETScSchurPC*>(p);
  KSPSolve(spc->outer_ksp, x, y);

  return 0;
}


PetscErrorCode PETScSchurPC::Destroy (PC pc)
{
  PETScSchurPC *shell;

  PCShellGetContext(pc,(void**)&shell);
  delete shell;

  return 0;
}
