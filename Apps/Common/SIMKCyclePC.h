//==============================================================================
//!
//! \file SIMKCyclePC.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Use a simulator as a preconditioner for K-refinement cycles.
//!
//==============================================================================
#ifndef SIM_K_CYCLE_PC_H_
#define SIM_K_CYCLE_PC_H_

#ifdef HAS_PETSC

#include "AlgEqSystem.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "SAM.h"
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/trivariate/SplineVolume.h>


/*!
 \brief Template for K-cycle preconditioning.
 \details Currently only single patch models.
*/

template<class Sim>
class SIMKCyclePC : public PETScPC
{
public:
  //! \brief Default constructor.
  //! \param fromSim The model for the solution
  //! \param pcSim The model for the preconditioner
  SIMKCyclePC(Sim& fromSim, Sim& pcSim) :
    S1(pcSim), fSim(fromSim),
    B(SparseMatrix::SUPERLU), BT(SparseMatrix::SUPERLU)
  {
    setup();
  }

  //! \brief Setup the restriction/prolongation operators.
  //! \details Operators are constructed through Greville point projection.
  void setup()
  {
    if (fSim.getPatch(1)->getNoSpaceDim() == 2) {
      const ASMs2D* fpch = dynamic_cast<const ASMs2D*>(fSim.getPatch(1));
      const ASMs2D* pch = dynamic_cast<const ASMs2D*>(S1.getPatch(1));
      std::array<RealArray,2> gPrm;
      fpch->getGrevilleParameters(gPrm[0], 0);
      fpch->getGrevilleParameters(gPrm[1], 1);
      N.resize(gPrm[0].size()*gPrm[1].size(), S1.getNoDOFs());
      NT.resize(S1.getNoDOFs(), gPrm[0].size()*gPrm[1].size());
      std::vector<Go::BasisDerivsSf> spline;
      pch->getBasis(1)->computeBasisGrid(gPrm[0], gPrm[1], spline);
      int p1 = pch->getBasis(1)->order_u();
      int p2 = pch->getBasis(1)->order_v();
      int n1 = pch->getBasis(1)->numCoefs_u();
      int n2 = pch->getBasis(1)->numCoefs_v();
      for (size_t ip = 0; ip < gPrm[0].size()*gPrm[1].size(); ++ip) {
        IntVec idx;
        ASMs2D::scatterInd(n1,n2,p1,p2,spline[ip].left_idx,idx);
        for (size_t k = 0; k < spline[ip].basisValues.size(); ++k) {
          N(ip+1, idx[k]+1) += spline[ip].basisValues[k];
          NT(idx[k]+1, ip+1) += spline[ip].basisValues[k];
        }
      }
      fpch->getBasis(1)->computeBasisGrid(gPrm[0], gPrm[1], spline);
      p1 = fpch->getBasis(1)->order_u();
      p2 = fpch->getBasis(1)->order_v();
      n1 = fpch->getBasis(1)->numCoefs_u();
      n2 = fpch->getBasis(1)->numCoefs_v();
      B.resize(fSim.getNoDOFs(), fSim.getNoDOFs());
      BT.resize(fSim.getNoDOFs(), fSim.getNoDOFs());

      for (size_t ip = 0; ip < gPrm[0].size()*gPrm[1].size(); ++ip) {
        IntVec idx;
        ASMs2D::scatterInd(n1,n2,p1,p2,spline[ip].left_idx,idx);
        for (size_t k = 0; k < spline[ip].basisValues.size(); ++k) {
          B(ip+1, idx[k]+1) += spline[ip].basisValues[k];
          BT(idx[k]+1, ip+1) += spline[ip].basisValues[k];
        }
      }
    } else {
      const ASMs3D* fpch = dynamic_cast<const ASMs3D*>(fSim.getPatch(1));
      const ASMs3D* pch = dynamic_cast<const ASMs3D*>(S1.getPatch(1));
      std::array<RealArray,3> gPrm;
      fpch->getGrevilleParameters(gPrm[0], 0);
      fpch->getGrevilleParameters(gPrm[1], 1);
      fpch->getGrevilleParameters(gPrm[2], 2);
      N.resize(gPrm[0].size()*gPrm[1].size()*gPrm[2].size(), S1.getNoDOFs());
      NT.resize(S1.getNoDOFs(), gPrm[0].size()*gPrm[1].size()*gPrm[2].size());
      std::vector<Go::BasisDerivs> spline;
      pch->getBasis(1)->computeBasisGrid(gPrm[0], gPrm[1], gPrm[2], spline);
      int p1 = pch->getBasis(1)->order(0);
      int p2 = pch->getBasis(1)->order(1);
      int p3 = pch->getBasis(1)->order(2);
      int n1 = pch->getBasis(1)->numCoefs(0);
      int n2 = pch->getBasis(1)->numCoefs(1);
      int n3 = pch->getBasis(1)->numCoefs(2);
      for (size_t ip = 0; ip < gPrm[0].size()*gPrm[1].size()*gPrm[2].size(); ++ip) {
        IntVec idx;
        ASMs3D::scatterInd(n1,n2,n3,p1,p2,p3,spline[ip].left_idx,idx);
        for (size_t k = 0; k < spline[ip].basisValues.size(); ++k) {
          N(ip+1, idx[k]+1) += spline[ip].basisValues[k];
          NT(idx[k]+1, ip+1) += spline[ip].basisValues[k];
        }
      }
      fpch->getBasis(1)->computeBasisGrid(gPrm[0], gPrm[1], gPrm[2], spline);
      p1 = fpch->getBasis(1)->order(0);
      p2 = fpch->getBasis(1)->order(1);
      p3 = fpch->getBasis(1)->order(2);
      n1 = fpch->getBasis(1)->numCoefs(0);
      n2 = fpch->getBasis(1)->numCoefs(1);
      n3 = fpch->getBasis(1)->numCoefs(2);
      B.resize(fSim.getNoDOFs(), fSim.getNoDOFs());
      BT.resize(fSim.getNoDOFs(), fSim.getNoDOFs());

      for (size_t ip = 0; ip < gPrm[0].size()*gPrm[1].size()*gPrm[2].size(); ++ip) {
        IntVec idx;
        ASMs3D::scatterInd(n1,n2,n3,p1,p2,p3,spline[ip].left_idx,idx);
        for (size_t k = 0; k < spline[ip].basisValues.size(); ++k) {
          B(ip+1, idx[k]+1) += spline[ip].basisValues[k];
          BT(idx[k]+1, ip+1) += spline[ip].basisValues[k];
        }
      }
    }
  }

  //! \brief Evaluate preconditioner.
  bool eval(Vec& x, Vec& y) override
  {
    // copy vector into a standard vector for expansion
    int len;
    VecGetSize(x, &len);
    StdVector X(len);
    std::vector<PetscInt> idx(len);
    std::iota(idx.begin(), idx.end(), 0);
    VecGetValues(x, len, idx.data(), X.getPtr());

    // Expand solution vector x to DOF vector
    StdVector dofSol(fSim.getNoDOFs());
    fSim.getSAM()->expandSolution(X, dofSol);

    // Evaluate restriction - N^T * (B^T)^-1 * dofSol
    BT.solve(dofSol);
    StdVector restrict(NT.rows());
    NT.multiply(dofSol, restrict);

    // copy to solution vector
    StdVector* eqsVec = S1.getVector();
    eqsVec->init();
    vectorToSolVec(S1, *eqsVec, restrict);
    eqsVec->beginAssembly();
    eqsVec->endAssembly();

    // Perform inner solve, store in restrict
    double cond;
    int oldLev = this->S1.msgLevel;
    this->S1.msgLevel = 0;
    if (!this->S1.solveSystem(restrict, 0, &cond, "displacement", false))
      return false;
    this->S1.msgLevel = oldLev;

    // Evaluate prolonged solution - dofSol = B^-1 * N * restrict
    N.multiply(restrict, dofSol);
    B.solve(dofSol);

    // Inject into solution vector y
    const int* meqn = fSim.getSAM()->getMEQN();
    for (size_t i = 1; i <= fSim.getPatch(1)->getNoNodes(); ++i) {
      int eq = meqn[i-1];
      if (eq > 0)
        VecSetValue(y, eq-1, dofSol(i), INSERT_VALUES);
    }

    return true;
  }

  //! \brief Evaluate the right-hand-side system - the prolongiation of the coarse system solution.
  //! \param rhs The resulting vector
  void evalRHS(PETScVector& rhs)
  {
    // Solve coarse system, store in solCoarse
    StdVector solCoarse(S1.getNoDOFs());
    S1.solveSystem(solCoarse);

    // Evaluate prolongiation - rhs = B^-1 * N * solCoarse
    StdVector prolong(N.rows());
    N.multiply(solCoarse, prolong);
    B.solve(prolong);
    vectorToSolVec(fSim, rhs, prolong);

    rhs.beginAssembly();
    rhs.endAssembly();
  }

  //! \brief Returns name of this preconditioner.
  const char* getName() const override { return "K-cycle preconditioner"; }

protected:
  //! \brief Helper template for copying a equation vector to a DOF vector.
  //! \param sim The simulator to use
  //! \param V1 The resulting DOF vector
  //! \param V2 The initial equation vecto
  template <class V1, class V2>
  void vectorToSolVec(Sim& sim, V1& dof, const V2& eqn)
  {
    const int* meqn = sim.getSAM()->getMEQN();
    for (size_t i = 1; i <= sim.getPatch(1)->getNoNodes(); ++i) {
      int eq = meqn[i-1];
      if (eq > 0)
        dof(eq) = eqn(i);
    }
  }

  Sim& S1;   //!< Preconditioner simulator
  Sim& fSim; //!< Solution model
  SparseMatrix B, BT; //!< Greville mass matrix for fSim
  SparseMatrix N; //!< Basis functions of S1 in greville of fSim
  SparseMatrix NT; //!< Stupid no-transpose multiply
};


/*!
  \brief Driver class for K-refined NURBS-based FEM analysis of linear problems.
*/

template<class T1>
class SIMKRefWrap : public T1
{
public:
  //! \brief Default constructor.
  explicit SIMKRefWrap(const typename T1::SetupProps& props) : T1(props) {}

  //! \brief Clears out patches, FE model and linear system.
  void clearModel()
  {
    this->clearProperties();
    delete this->mySam;
    this->mySam=nullptr;
    for (ASMbase* patch : this->myModel)
      delete patch;
    this->myModel.clear();
    delete this->myEqSys;
    this->myEqSys = nullptr;
  }

  //! \brief Set a linear solver parameter.
  //! \param key Parameter to set
  //! \param value Value for parameter
  //! \param old If non-null, previous value is stored here
  void setLinSolParam(const std::string& key, const std::string& value, std::string* old = nullptr)
  {
    if (!this->mySolParams)
      this->mySolParams = new LinSolParams;

    if (old)
      *old = this->mySolParams->getStringValue(key);

    this->mySolParams->addValue(key, value);
  }

  //! \brief Remove shell matrix-vector product.
  //! \param oldSolver Solver used before it was overridden with PETSc for shell solver
  //! \param oldRtol Old relative tolerance before it was overridden with shell solver tolerance
  void clearMxV(int oldSolver, const std::string& oldRtol)
  {
    this->mySolParams->addValue("matrixfree", "0");
    this->mySolParams->addValue("rtol", oldRtol);
    this->getMatrix()->clearMxV();

    // We have to flag that we want to solve using sparse direct solver for preconditioner.
    // We have assembled into a PETSc matrix previously.
    if (oldSolver != SystemMatrix::PETSC)
      static_cast<PETScMatrix*>(this->myEqSys->getMatrix(0))->setSolveSparse(true);
  }

  //! \brief Obtain a pointer to current PETSc system matrix.
  PETScMatrix* getMatrix()
  {
    return dynamic_cast<PETScMatrix*>(this->myEqSys->getMatrix(0));
  }

  //! \brief Obtain a pointer to current system vector.
  StdVector* getVector()
  {
    return dynamic_cast<StdVector*>(this->myEqSys->getVector(0));
  }

  //! \brief Assemble system.
  bool assembleStep()
  {
    return this->assembleSystem();
  }
};

#endif

#endif
