// $Id$
//==============================================================================
//!
//! \file SIMSolverKRef.h
//!
//! \date Feb 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Stationary solver class template for K-refinement.
//!
//==============================================================================

#ifndef _SIM_SOLVER_KREF_H_
#define _SIM_SOLVER_KREF_H_

#include "SIMSolver.h"
#include "SIMenums.h"
#include "SIMMxV.h"
#include "SIMKCyclePC.h"
#include "Utilities.h"


/*!
  \brief Template class for stationary K-refinement simulator drivers.
  \details This template can be instanciated over any type implementing the
  ISolver interface. It provides a solver loop with data output.
*/

template<class T1> class SIMSolverKRef : public SIMSolverStat<SIMKRefWrap<T1>>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMSolverKRef(SIMKRefWrap<T1>& s1, SIMKRefWrap<T1>& s2) :
    SIMSolverStat<SIMKRefWrap<T1>>(s1), S2(s2)
  {
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverKRef() = default;

  //! \brief Read k-refinement settings.
  bool read(const char* infile) override { return this->SIMadmin::read(infile); }

  //! \brief Solves the problem (all k-refinement cycles).
  int solveProblem(char* infile,
                   const char* heading = nullptr, bool = false) override
  {
    int geoBlk=0, nBlock=0;

    this->printHeading(heading);

    this->S1.getProcessAdm().cout <<"\nK-refinement cycle 1" << std::endl;

    TimeStep tp;
    if (!this->initSystem(infile, this->S2, false))
      return 1;

    SIMKRefWrap<T1>* s1 = &this->S1;
    SIMKRefWrap<T1>* s2 = &this->S2;

    // First assemble first preconditioner system.
    s2->initDirichlet();
    s2->assembleStep();

    for (tp.iter = 1; tp.iter <= cycles; ++tp.iter) {
      if (tp.iter != 1)
        this->S1.getProcessAdm().cout <<"\nK-refinement cycle "<< tp.iter << std::endl;

      // Reload main simulator
      s1->increaseAdditionalRefinement();
      if (!this->initSystem(infile,*s1))
        return 3;

      SIMMxV<SIMKRefWrap<T1>> mxv(*s1);
      SIMKCyclePC<SIMKRefWrap<T1>> pc(*s1,*s2);
      mxv.setPC(&pc);

      // First RHS is solution of coarse system prolonged onto fine mesh.
      PETScVector rhs(s1->getProcessAdm(),
                      s1->getSAM()->getNoEquations());
      s2->initDirichlet();
      pc.evalRHS(rhs);
      s1->getMatrix()->setInitialGuess(&rhs);
      s2->updateDirichlet(0.0);

      s1->initDirichlet();
      if (!mxv.solveStep())
        return 1;

      if (!s1->saveModel(tp.iter == 1 ? infile : nullptr, geoBlk, nBlock))
        return 2;
      else if (!s1->saveStep(tp, nBlock))
        return 2;

      // Reset linear solver parameters
      s1->clearMxV(oldSolver, oldRtol);

      // Increase additional refinement for s2 (s1 in next step) and clear out model.
      s2->increaseAdditionalRefinement();
      s2->clearModel();
      s2->setVTF(s1->getVTF());

      // Main solver will be preconditioner in next cycle
      std::swap(s1, s2);
    }

    this->S2.setVTF(nullptr);

    return 0;
  }

  //! \brief Load and initialize a simulator.
  //! \param infile Input file to use
  //! \param sim The simulator to initialize
  //! \param setPetsc \e true to force a petsc solver
  int initSystem(const char* infile, SIMKRefWrap<T1>& sim, bool setPetsc=true)
  {
    if (!sim.read(infile))
      return 1;

    if (setPetsc)
      sim.opt.solver = SystemMatrix::PETSC;
    else
      oldSolver = sim.opt.solver;

    if (!sim.preprocess())
      return 2;

    sim.setQuadratureRule(sim.opt.nGauss[0],true);
    if (setPetsc) {
      sim.setLinSolParam("matrixfree", "1");
      std::stringstream str;
      str << rtol;
      sim.setLinSolParam("rtol", str.str(), &oldRtol);
    }

    sim.setMode(SIM::STATIC);
    return sim.initSystem(sim.opt.solver,1,1,0,true);
  }

protected:
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override
  {
    if (strcasecmp(elem->Value(),"krefinement"))
      return this->SIMSolverStat<SIMKRefWrap<T1>>::parse(elem);

    if (utl::getAttribute(elem,"cycles",cycles))
      IFEM::cout << "\tDoing " << cycles << " k-refinement cycles." << std::endl;
    if (utl::getAttribute(elem,"rtol",rtol))
      IFEM::cout << "\tK-solver tolerance " << rtol << std::endl;

    return true;
  }

  int cycles = 1;      //!< Number of k-refinement cycles
  double rtol = 1e-6;  //!< Relative tolerance in outer petsc solver
  int oldSolver = -1;  //!< Old matrix type before we overwrote with PETSc for k-refinement
  std::string oldRtol; //!< Old relative tolerance before we overwrite with k-refinement tolerance
  SIMKRefWrap<T1>& S2; //!< Reference to second simulator
};

#endif
