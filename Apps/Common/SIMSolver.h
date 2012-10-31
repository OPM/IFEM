//==============================================================================
//!
//! \file SIMSolver.h
//!
//! \date Oct 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief SIM solver class template
//==============================================================================
#ifndef SIM_SOLVER_H_
#define SIM_SOLVER_H_

#include "DataExporter.h"
#include "TimeStep.h"

/*!
  \brief Solver interface. This template can be instanced over any type
         implementing the ISolver interface. It provides a time stepping
         loop with data output.
*/
  template<class T1>
class SIMSolver
{
public:
  //! \brief Constructor 
  //! \param[in] s1 Pointer to model solver simulator
  SIMSolver(T1& s1) : S1(s1)
  {
  }

  //! \brief Destructor
  virtual ~SIMSolver()
  {
  }

  //! \brief Solves the problem up to the final time.
  virtual bool solveProblem(char* infile, DataExporter* exporter = NULL)
  {
    // Save initial step to VTF
    TimeStep tp = S1.getTimePrm();

    // Save FE model to VTF file for visualization
    int nBlock;
    if (!S1.saveModel(infile, nBlock))
      return 3;

    // Save initial step
    if (!S1.saveStep(tp, nBlock))
      return false;

    // Solve for each time step up to final time
    for (int iStep = 1; S1.advanceStep(tp); iStep++)
    {
      if (!S1.solveStep(tp))
        return false;
      if (!S1.saveStep(tp, nBlock))
        return false;
      if (exporter)
        exporter->dumpTimeLevel(&tp);
    }

    return true;
  }
protected:
  T1& S1; //!< Solver
};

#endif
