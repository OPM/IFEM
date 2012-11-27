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
class SIMSolver : public SIMinput
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

  bool advanceStep()
  {
    return S1.advanceStep(tp) && tp.increment();
  }

  void fastForward(int steps)
  {
    for (int i=0;i<steps;++i)
      advanceStep();
  }

  //! \brief Solves the problem up to the final time.
  virtual bool solveProblem(char* infile, DataExporter* exporter = NULL, int maxIter = 1)
  {
    // Save FE model to VTF file for visualization
    int nBlock;
    if (!S1.saveModel(infile, nBlock))
      return 3;

    // Save initial step
    if (!S1.saveStep(tp, nBlock))
      return false;

    // Solve for each time step up to final time
    for (int iStep = 1; advanceStep(); iStep++)
    {
      // Possible non-linear iteration loop
      int iter = 0;
      bool converged = false;
      while ((iter < maxIter) && (!converged)) {
	converged = S1.solveStep(tp);
	iter++;
      }

      // Return if not converged
      if (!converged) 
	return false;

      // Save solution
      if (!S1.saveStep(tp, nBlock))
        return false;
      if (exporter)
        exporter->dumpTimeLevel(&tp);
    }

    return true;
  }

  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyWord, std::istream& is)
  {
    return false;
  }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"timestepping") == 0)
      tp.parse(elem);

    return true;
  }

  const TimeStep& getTimePrm() const { return tp; }

protected:
  TimeStep tp; //<! Time stepping information
  T1& S1; //!< Solver
};

#endif
