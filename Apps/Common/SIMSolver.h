// $Id$
//==============================================================================
//!
//! \file SIMSolver.h
//!
//! \date Oct 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief SIM solver class template.
//!
//==============================================================================

#ifndef _SIM_SOLVER_H_
#define _SIM_SOLVER_H_

#include "SIMinput.h"
#include "DataExporter.h"
#include "TimeStep.h"
#include "IFEM.h"


//! \brief Structure for configuring a given simulator
//! \details Your SIM need to specialize this for its type
template<class T>
struct SolverConfigurator {
  //! \brief Configure a simulator
  //! \param sim The simulator to configure
  //! \param props The setup properties for the simulator
  //! \param infile The input file to parse
  int setup(T& sim, typename T::SetupProps& props, char* infile);
};


//! \brief Configuration template
template<class T>
int ConfigureSIM(T& t, typename T::SetupProps& p, char* infile)
{
  SolverConfigurator<T> setup;
  return setup.setup(t, p, infile);
}


/*!
  \brief Template class for simulator drivers.
  \details This template can be instanciated over any type implementing the
  ISolver interface. It provides a time stepping loop with data output.
*/

template<class T1> class SIMSolver : public SIMinput
{
public:
  //! \brief The constructor initializes the reference to the actual solver.
  SIMSolver(T1& s1) : SIMinput("Time integration driver"), S1(s1) {}
  //! \brief Empty destructor.
  virtual ~SIMSolver() {}

  //! \brief Advances the time step one step forward.
  bool advanceStep() { return tp.increment() && S1.advanceStep(tp); }
  //! \brief Advances the time step \a n steps forward.
  void fastForward(int n) { for (int i = 0; i < n; i++) this->advanceStep(); }

  //! \brief Solves the problem up to the final time.
  bool solveProblem(char* infile, DataExporter* exporter, const char* heading)
  {
    return this->solveProblem(infile,exporter,1,heading);
  }

  //! \brief Solves the problem up to the final time.
  virtual bool solveProblem(char* infile, DataExporter* exporter = NULL,
                            int maxIter = 1, const char* heading = NULL)
  {
    int geoBlk = 0;
    int nBlock = 0;

    // Save FE model to VTF file for visualization
    if (!S1.saveModel(infile,geoBlk,nBlock))
      return false;

    // Save the initial configuration to VTF file
    if (!S1.saveStep(tp,nBlock))
      return false;

    if (maxIter < 1)
      return true; // Data check, finish up without solving anything

    if (heading && myPid == 0)
    {
      // Write an application-specific heading, if provided
      std::string myHeading(heading);
      size_t n = myHeading.find_last_of('\n');
      if (n+1 < myHeading.size()) n = myHeading.size()-n;
      std::cout <<"\n\n" << myHeading <<"\n";
      for (size_t i = 0; i < n && i < myHeading.size(); i++) std::cout <<'=';
      std::cout << std::endl;
    }

    // Solve for each time step up to final time
    for (int iStep = 1; this->advanceStep(); iStep++)
    {
      // Possible staggered iteration loop
      bool converged = false;
      for (int iter = 0; iter < maxIter && !converged; iter++)
	converged = S1.solveStep(tp);

      // Return if not converged
      if (!converged)
	return false;

      // Save solution
      if (!S1.saveStep(tp,nBlock))
        return false;

      if (exporter)
        exporter->dumpTimeLevel(&tp);
      IFEM::pollControllerFifo();
    }

    return true;
  }

  //! \brief Parses a data section from an input stream.
  virtual bool parse(char*, std::istream&) { return false; }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"timestepping") == 0)
      return tp.parse(elem);

    return true;
  }

  //! \brief Returns a reference to the time stepping information.
  const TimeStep& getTimePrm() const { return tp; }

protected:
  TimeStep tp; //!< Time stepping information
  T1&      S1; //!< The actual solver
};

#endif
