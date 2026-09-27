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

#include <cstring>
#include <iostream>


PETScSchurPC::PETScSchurPC (PC& pc_init, const std::vector<Mat>& blocks,
                            const LinSolParams::BlockParams& params,
                            const ProcessAdm& adm, int verbosity)
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
  if (verbosity > 0)
    KSPView(inner_ksp, PETSC_VIEWER_STDOUT_(*adm.getCommunicator()));

  KSPCreate(*adm.getCommunicator(), &outer_ksp);
  KSPGetPC(outer_ksp, &pc);

  // Applying the Schur complement means a Krylov solve against it, and what
  // that costs decides whether the preconditioner is worth having at all.
  // Unpreconditioned it takes hundreds of iterations, every one of them an
  // inner solve of the momentum operator, so it is cut short long before it
  // converges and the outer solver is handed something inconsistent.
  //
  // Whoever set this preconditioner up may have left an approximation of the
  // Schur complement on it. It is a matrix over the same equations, and
  // whatever it is worth as an approximation it is worth more as something
  // to precondition the solve against the real thing with.
  Mat sA = nullptr, sP = nullptr;
  MatType ptype = nullptr;
  PCGetOperators(pc_init, &sA, &sP);
  if (sP)
    MatGetType(sP, &ptype);

  // Only an assembled one is of any use: where the Schur complement is
  // preconditioned with itself, what comes back is the operator again and
  // there is nothing in it to factorize.
  const bool preconditioned = sP && ptype &&
                              strcmp(ptype, MATSHELL) &&
                              strcmp(ptype, MATSCHURCOMPLEMENT);

  // What is applied above is the negative of the Schur complement, that
  // being the positive definite one of the two, while the approximation
  // left here is of the Schur complement itself. Preconditioning the one
  // with the other takes a change of sign, which an operator with a sign
  // to speak of does not survive being without.
  if (preconditioned)
  {
    MatDuplicate(sP, MAT_COPY_VALUES, &prec);
    MatScale(prec, -1.0);
  }

  PCSetType(pc, preconditioned ? PCILU : PCNONE);
  if (preconditioned)
    // A constraint on the integrated pressure leaves its multiplier in these
    // equations, and the row it sits in has nothing on the diagonal for a
    // factorization to work with.
    PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
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
  PetscInt m, n;
  MatGetLocalSize(blocks[3], &m, &n);
  MatSetSizes(outer_mat, m, n, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetFromOptions(outer_mat);
  MatSetType(outer_mat, MATSHELL);
  MatShellSetContext(outer_mat, this);
  MatShellSetOperation(outer_mat, MATOP_MULT,
                       (void(*)(void))&PETScSchurPC::Apply_Schur);
  MatSetUp(outer_mat);

  KSPSetOperators(outer_ksp, outer_mat,
                  preconditioned ? prec : outer_mat);
  KSPSetFromOptions(outer_ksp);
  KSPSetUp(outer_ksp);
  if (verbosity > 0)
    KSPView(outer_ksp, PETSC_VIEWER_STDOUT_(*adm.getCommunicator()));
  PCSetUp(pc_init);

  // The inner solve works on the rows of the momentum operator
  MatCreateVecs(blocks[0], nullptr, &tmp);
}


PETScSchurPC::~PETScSchurPC ()
{
  KSPDestroy(&inner_ksp);
  KSPDestroy(&outer_ksp);
  MatDestroy(&outer_mat);
  if (prec) MatDestroy(&prec);
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
