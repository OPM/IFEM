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
#include "tinyxml.h"


//! \brief Structure for configuring a given simulator
//! \details Your SIM need to specialize this for its type
template<class T>
struct SolverConfigurator {
  //! \brief Configure a simulator
  //! \param sim The simulator to configure
  //! \param props The setup properties for the simulator
  //! \param infile The input file to parse
  int setup(T& sim, const typename T::SetupProps& props, char* infile);
};


//! \brief Configuration template
template<class T>
int ConfigureSIM(T& t, char* infile,
                 const typename T::SetupProps& p = typename T::SetupProps())
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

  //! \brief Returns a const reference to the time stepping information.
  const TimeStep& getTimePrm() const { return tp; }

  //! \brief Advances the time step one step forward.
  bool advanceStep() { return tp.increment() && S1.advanceStep(tp); }
  //! \brief Advances the time step \a n steps forward.
  void fastForward(int n) { for (int i = 0; i < n; i++) this->advanceStep(); }

  //! \brief Postprocesses the solution of current time step.
  void postSolve(const TimeStep& t, bool rst = false) { S1.postSolve(t,rst); }

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, DataExporter* exporter = nullptr,
                           const char* heading = nullptr, bool saveInit = true)
  {
    // Save FE model to VTF and HDF5 for visualization
    // Optionally save the initial configuration also
    int geoBlk = 0, nBlock = 0;
    if (!this->saveState(exporter,geoBlk,nBlock,true,infile,saveInit))
      return 2;

    this->printHeading(heading);

    // Solve for each time step up to final time
    for (int iStep = 1; this->advanceStep(); iStep++)
      if (!S1.solveStep(tp))
        return 3;
      else if (!this->saveState(exporter,geoBlk,nBlock))
        return 4;
      else
        IFEM::pollControllerFifo();

    return 0;
  }

protected:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyw, std::istream& is) { return tp.parse(keyw,is); }
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) { return tp.parse(elem); }

  //! \brief Writes an application-specific heading, if provided
  void printHeading(const char* heading)
  {
    if (heading)
    {
      std::string myHeading(heading);
      size_t n = myHeading.find_last_of('\n');
      if (n+1 < myHeading.size()) n = myHeading.size()-n;
      IFEM::cout <<"\n\n"<< myHeading <<"\n";
      for (size_t i = 0; i < n && i < myHeading.size(); i++) IFEM::cout <<'=';
      IFEM::cout << std::endl;
    }
  }

  //! \brief Saves geometry and results to VTF and HDF5 for current time step.
  bool saveState(DataExporter* exporter, int& geoBlk, int& nBlock,
                 bool newMesh = false, char* infile = nullptr,
                 bool saveRes = true)
  {
    if (newMesh && !S1.saveModel(infile,geoBlk,nBlock))
      return false;

    if (saveRes && !S1.saveStep(tp,nBlock))
      return false;

    if (saveRes && exporter)
      exporter->dumpTimeLevel(&tp,newMesh);

    return true;
  }

  TimeStep tp; //!< Time stepping information
  T1&      S1; //!< The actual solver
};

#endif
