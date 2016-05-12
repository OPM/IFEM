// $Id$
//==============================================================================
//!
//! \file SAMpatchPETSc.C
//!
//! \date Sep 6 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Evaluation of global norms for distributed models.
//!
//==============================================================================

#include "SAMpatchPETSc.h"
#include "LinAlgInit.h"
#include "ASMstruct.h"
#include "Vec3.h"
#include "PETScMatrix.h"
#include "ProcessAdm.h"


SAMpatchPETSc::SAMpatchPETSc(const std::map<int,int>& g2ln, const ProcessAdm& padm) : adm(padm)
{
  nProc = adm.getNoProcs();
  LinAlgInit::increfs();
}


SAMpatchPETSc::~SAMpatchPETSc()
{
  for (auto& it : dofIS) {
    if (it.second.scatterCreated)
      VecScatterDestroy(&it.second.ctx);
    ISDestroy(&it.second.local);
    ISDestroy(&it.second.global);
  }
  dofIS.clear();
  LinAlgInit::decrefs();
}


bool SAMpatchPETSc::init(const ASMVec& model, int numNod)
{
  // Get local 2 global node mapping for each patch
  patch = model;

  return SAMpatch::init(model,numNod);
}


Real SAMpatchPETSc::dot (const Vector& x, const Vector& y, char dofType) const
{
#ifdef HAVE_MPI
  if (adm.isParallel()) {
    if (dofIS.find(dofType) == dofIS.end())
      setupIS(dofType);

    Vec lx, ly;
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.size(), x.data(), &lx);
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, y.size(), y.data(), &ly);
    Vec gx, gy;
    VecCreate(*adm.getCommunicator(), &gx);
    VecCreate(*adm.getCommunicator(), &gy);
    VecSetSizes(gx, dofIS[dofType].nDofs, PETSC_DETERMINE);
    VecSetSizes(gy, dofIS[dofType].nDofs, PETSC_DETERMINE);
    VecSetFromOptions(gx);
    VecSetFromOptions(gy);
    if (!dofIS[dofType].scatterCreated) {
      VecScatterCreate(lx, dofIS[dofType].local, gx, dofIS[dofType].global, &dofIS[dofType].ctx);
      dofIS[dofType].scatterCreated = true;
    }

    VecScatterBegin(dofIS[dofType].ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(dofIS[dofType].ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterBegin(dofIS[dofType].ctx, ly, gy, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(dofIS[dofType].ctx, ly, gy, INSERT_VALUES, SCATTER_FORWARD);

    PetscReal d;
    VecDot(gx, gy, &d);
    VecDestroy(&lx);
    VecDestroy(&gx);
    VecDestroy(&ly);
    VecDestroy(&gy);

    return d;
  }
#endif

  return this->SAM::dot(x,y,dofType);
}


Real SAMpatchPETSc::normL2(const Vector& x, char dofType) const
{
#ifdef HAVE_MPI
  if (adm.isParallel()) {
    if (dofIS.find(dofType) == dofIS.end())
      setupIS(dofType);

    Vec lx;
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.size(), x.data(), &lx);
    Vec gx;
    VecCreate(*adm.getCommunicator(), &gx);
    VecSetSizes(gx, dofIS[dofType].nDofs, PETSC_DETERMINE);
    VecSetFromOptions(gx);
    PetscInt n;
    VecGetSize(gx, &n);

    if (!dofIS[dofType].scatterCreated) {
      VecScatterCreate(lx, dofIS[dofType].local, gx, dofIS[dofType].global, &dofIS[dofType].ctx);
      dofIS[dofType].scatterCreated = true;
    }

    VecScatterBegin(dofIS[dofType].ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(dofIS[dofType].ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    PetscReal d;
    VecNorm(gx, NORM_2, &d);
    VecDestroy(&lx);
    VecDestroy(&gx);

    return d / sqrt(double(n));
  }
#endif

  return this->SAM::normL2(x, dofType);
}


Real SAMpatchPETSc::normInf(const Vector& x, size_t& comp, char dofType) const
{
  Real max = this->SAM::normInf(x, comp, dofType);
#ifdef HAVE_MPI
  if (adm.isParallel())
    max = adm.allReduce(max, MPI_MAX);
#endif

  return max;
}


void SAMpatchPETSc::setupIS(char dofType) const
{
  PetscIntVec ldofs;
  PetscIntVec gdofs;
  int gdof = 0;
  if (adm.getProcId() > 0)
    adm.receive(gdof, adm.getProcId()-1);

  for (size_t i = 0; i < adm.dd.getMLGN().size(); ++i) {
    if ((dofType == 'A' || nodeType.empty() || this->SAM::getNodeType(i+1) == dofType) &&
        adm.dd.getMLGN()[i] >= adm.dd.getMinNode() &&
        adm.dd.getMLGN()[i] <= adm.dd.getMaxNode()) {
      std::pair<int, int> dofs = this->SAM::getNodeDOFs(i+1);
      for (int dof = dofs.first; dof <= dofs.second; ++dof) {
        ldofs.push_back(dof-1);
        gdofs.push_back(gdof++);
      }
    }
  }

  if (adm.getProcId() < adm.getNoProcs()-1)
    adm.send(gdof, adm.getProcId()+1);

  ISCreateGeneral(*adm.getCommunicator(), ldofs.size(), ldofs.data(), PETSC_COPY_VALUES, &dofIS[dofType].local);
  ISCreateGeneral(*adm.getCommunicator(), gdofs.size(), gdofs.data(), PETSC_COPY_VALUES, &dofIS[dofType].global);

  dofIS[dofType].nDofs = gdof - gdofs.front();
}


bool SAMpatchPETSc::expandSolution(const SystemVector& solVec,
                                   Vector& dofVec, Real scaleSD) const
{
#ifdef HAVE_MPI
  PETScVector* Bptr = const_cast<PETScVector*>(dynamic_cast<const PETScVector*>(&solVec));
  if (!Bptr)
    return false;

  if (adm.isParallel()) {
    if (!glob2LocEq) {
      std::vector<int> mlgeq(adm.dd.getMLGEQ());
      for (auto& it : mlgeq)
        --it;

      ISCreateGeneral(*adm.getCommunicator(),adm.dd.getMLGEQ().size(),
                      mlgeq.data(), PETSC_COPY_VALUES, &glob2LocEq);
    }

    Vec solution;
    VecCreateSeq(PETSC_COMM_SELF, Bptr->dim(), &solution);
    VecScatter ctx;
    VecScatterCreate(Bptr->getVector(), glob2LocEq, solution, nullptr, &ctx);
    VecScatterBegin(ctx, Bptr->getVector(), solution, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, Bptr->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
    PetscScalar* data;
    VecGetArray(solution, &data);
    std::copy(data, data + Bptr->dim(), Bptr->getPtr());
    VecDestroy(&solution);
  } else {
    PetscScalar* data;
    VecGetArray(Bptr->getVector(), &data);
    std::copy(data, data + Bptr->dim(), Bptr->getPtr());
    VecRestoreArray(Bptr->getVector(), &data);
  }
#endif

  return this->SAMpatch::expandSolution(solVec, dofVec, scaleSD);
}
