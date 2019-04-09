// $Id$
//==============================================================================
//!
//! \file SAMpatchPETSc.C
//!
//! \date Aug 31 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Evaluation of global norms for distributed models.
//!
//==============================================================================

#include "SAMpatchPETSc.h"
#include "LinAlgInit.h"
#include "PETScMatrix.h"
#include "ProcessAdm.h"


SAMpatchPETSc::SAMpatchPETSc(const std::map<int,int>& g2ln,
                             const ProcessAdm& padm) : adm(padm)
{
  nProc = adm.getNoProcs();
  LinAlgInit::increfs();
}


SAMpatchPETSc::~SAMpatchPETSc()
{
#ifdef HAVE_MPI
  for (auto& it : dofIS) {
    if (it.second.scatterCreated)
      VecScatterDestroy(&it.second.ctx);
    ISDestroy(&it.second.local);
    ISDestroy(&it.second.global);
  }
  dofIS.clear();
  if (glob2LocEq)
    ISDestroy(&glob2LocEq);
#endif
  LinAlgInit::decrefs();
}


bool SAMpatchPETSc::init (const std::vector<ASMbase*>& model, int numNod,
                          const std::vector<char>& dTypes)
{
  patch = model;
  return this->SAMpatch::init(model,numNod,dTypes);
}


Real SAMpatchPETSc::dot (const Vector& x, const Vector& y, char dofType) const
{
#ifdef HAVE_MPI
  if (adm.isParallel()) {
    if (adm.dd.isPartitioned())
      return this->SAMpatch::dot(x,y,dofType);

    DofIS& dis = this->getIS(dofType);

    Vec lx, ly;
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.size(), x.data(), &lx);
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, y.size(), y.data(), &ly);
    Vec gx, gy;
    VecCreate(*adm.getCommunicator(), &gx);
    VecCreate(*adm.getCommunicator(), &gy);
    VecSetSizes(gx, dis.nDofs, PETSC_DETERMINE);
    VecSetSizes(gy, dis.nDofs, PETSC_DETERMINE);
    VecSetFromOptions(gx);
    VecSetFromOptions(gy);
    if (!dis.scatterCreated) {
      VecScatterCreate(lx, dis.local, gx, dis.global, &dis.ctx);
      dis.scatterCreated = true;
    }

    VecScatterBegin(dis.ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(dis.ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterBegin(dis.ctx, ly, gy, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(dis.ctx, ly, gy, INSERT_VALUES, SCATTER_FORWARD);

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


Real SAMpatchPETSc::normL2 (const Vector& x, char dofType) const
{
#ifdef HAVE_MPI
  if (adm.isParallel() && !adm.dd.isPartitioned()) {
    DofIS& dis = this->getIS(dofType);

    Vec lx;
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.size(), x.data(), &lx);
    Vec gx;
    VecCreate(*adm.getCommunicator(), &gx);
    VecSetSizes(gx, dis.nDofs, PETSC_DETERMINE);
    VecSetFromOptions(gx);
    PetscInt n;
    VecGetSize(gx, &n);

    if (!dis.scatterCreated) {
      VecScatterCreate(lx, dis.local, gx, dis.global, &dis.ctx);
      dis.scatterCreated = true;
    }

    VecScatterBegin(dis.ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(dis.ctx, lx, gx, INSERT_VALUES, SCATTER_FORWARD);
    PetscReal d;
    VecNorm(gx, NORM_2, &d);
    VecDestroy(&lx);
    VecDestroy(&gx);

    return d / sqrt(double(n));
  }
#endif

  return this->SAM::normL2(x,dofType);
}


Real SAMpatchPETSc::normInf (const Vector& x, size_t& comp, char dofType) const
{
  Real max = this->SAM::normInf(x,comp,dofType);
#ifdef HAVE_MPI
  if (adm.isParallel() && !adm.dd.isPartitioned())
    max = adm.allReduce(max, MPI_MAX);
  // TODO: comp (node index of the max value) is not necessarily correct here
#endif

  return max;
}


#ifdef HAVE_MPI
SAMpatchPETSc::DofIS& SAMpatchPETSc::getIS (char dofType) const
{
  std::map<char,DofIS>::iterator it = dofIS.find(dofType);
  if (it != dofIS.end())
    return it->second;

  PetscIntVec ldofs;
  PetscIntVec gdofs;
  int gdof = 0;
  if (adm.getProcId() > 0)
    adm.receive(gdof, adm.getProcId()-1);

  size_t inod = 0;
  for (int n : adm.dd.getMLGN())
  {
    if (dofType == 'A' || inod >= nodeType.size() || nodeType[inod] == dofType)
      if (n >= adm.dd.getMinNode() && n <= adm.dd.getMaxNode())
        for (int dof = madof[inod]; dof < madof[inod+1]; dof++)
          if (dofType == 'A' || dof_type.empty() || dof_type[dof-1] == dofType)
          {
            ldofs.push_back(dof-1);
            gdofs.push_back(gdof++);
          }
    ++inod;
  }

  if (adm.getProcId() < adm.getNoProcs()-1)
    adm.send(gdof, adm.getProcId()+1);

  DofIS& newIS = dofIS[dofType];
  ISCreateGeneral(*adm.getCommunicator(), ldofs.size(), ldofs.data(), PETSC_COPY_VALUES, &newIS.local);
  ISCreateGeneral(*adm.getCommunicator(), gdofs.size(), gdofs.data(), PETSC_COPY_VALUES, &newIS.global);
  newIS.nDofs = gdof - gdofs.front();
  return newIS;
}
#endif


bool SAMpatchPETSc::expandSolution(const SystemVector& solVec,
                                   Vector& dofVec, Real scaleSD) const
{
  PETScVector* Bptr = const_cast<PETScVector*>(dynamic_cast<const PETScVector*>(&solVec));
  if (!Bptr)
    return false;

#ifdef HAVE_MPI
  if (adm.isParallel()) {
    if (!glob2LocEq && !adm.dd.isPartitioned()) {
      IntVec mlgeq(adm.dd.getMLGEQ());
      for (int& ieq : mlgeq) --ieq;
      ISCreateGeneral(*adm.getCommunicator(),adm.dd.getMLGEQ().size(),
                      mlgeq.data(), PETSC_COPY_VALUES, &glob2LocEq);
    }

   VecScatter ctx;
   Vec solution;
   if (adm.dd.isPartitioned())
      VecScatterCreateToAll(Bptr->getVector(), &ctx, &solution);
    else {
      VecCreateSeq(PETSC_COMM_SELF, Bptr->dim(), &solution);
      VecScatterCreate(Bptr->getVector(), glob2LocEq, solution, nullptr, &ctx);
   }

    VecScatterBegin(ctx, Bptr->getVector(), solution, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, Bptr->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
    PetscScalar* data;
    VecGetArray(solution, &data);
    if (adm.dd.isPartitioned()) {
      for (const auto& it : adm.dd.getG2LEQ(0))
        (*Bptr)(it.second) = data[it.first-1];
    } else
      std::copy(data, data + Bptr->dim(), Bptr->getPtr());

    VecDestroy(&solution);
  } else
#endif
  {
    PetscScalar* data;
    VecGetArray(Bptr->getVector(), &data);
    std::copy(data, data + Bptr->dim(), Bptr->getPtr());
    VecRestoreArray(Bptr->getVector(), &data);
  }

  return this->SAMpatch::expandSolution(solVec, dofVec, scaleSD);
}
