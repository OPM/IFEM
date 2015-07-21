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
#include "AdaptiveSIM.h"


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
int ConfigureSIM(T& t, char* infile, const typename T::SetupProps& p=typename T::SetupProps())
{
  SolverConfigurator<T> setup;
  return setup.setup(t, p, infile);
}


//! \brief Structure for handling a time step for a given simulator
template<class T>
struct SolverHandler {
  //! \brief Constructor.
  SolverHandler(T& sim, char* infile) : m_sim(sim), m_infile(infile) {}

  bool pre(TimeStep& tp, int& geoBlk, int& nBlock)
  {
    // Save FE model to VTF file for visualization
    return m_sim.saveModel(m_infile,geoBlk,nBlock) &&
           m_sim.saveStep(tp,nBlock);
  }

  //! \brief Solve for a single time step.
  //! \param tp Time stepping parameters.
  //! \param nBlock Running VTF block counter.
  int solve(TimeStep& tp, int& nBlock, DataExporter* exporter)
  {
    if (!m_sim.solveStep(tp))
      return 1;

    // Save solution
    if (!m_sim.saveStep(tp,nBlock))
      return 1;

    if (exporter)
      exporter->dumpTimeLevel(&tp);

    return 0;
  }

  //! \brief Post step for a single time step.
  bool post(TimeStep& tp)
  {
    return false;
  }

  T& m_sim; //!< Reference to simulator to solve for.
  char* m_infile; //!< Input file to parse.
};


//! \brief Specialization for adaptive simulators.
template<>
struct SolverHandler<AdaptiveSIM> {
  //! \brief Constructor.
  SolverHandler<AdaptiveSIM>(AdaptiveSIM& sim, char* infile) :
    m_sim(sim), m_infile(infile), m_step(0)
  {
    m_norms = sim.getNoNorms();
  }

  bool pre(TimeStep& tp, int& geoBlk, int& nBlock)
  {
    return true;
  }

  //! \brief Solve for a adaptive level.
  //! \param tp Time stepping parameters.
  //! \param nBlock Running VTF block counter.
  int solve(TimeStep& tp, int& nBlock, DataExporter* exporter)
  {
    if (!m_sim.solveStep(m_infile, m_step))
      return 5;

    if (!m_sim.writeGlv(m_infile,m_step,m_norms))
      return 6;

    if (exporter)
      exporter->dumpTimeLevel(NULL, true);

    return 0;
  }

  //! \brief Post step for a single adaptive loop.
  bool post(TimeStep& tp)
  {
    ++m_step;

    if (m_step != 1 && !m_sim.adaptMesh(m_step)) {
      m_step = 0;
      tp.step--; // advanceStep will increment it
      return false;
    }

    return true;
  }

  AdaptiveSIM& m_sim; //!< Reference to adaptive simulator.
  char* m_infile;     //!< Input file to parse.
  size_t m_step;      //!< Current adaptive level.
  size_t m_norms;     //!< Number of norms in adaptive group.
};


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

  void postSolve(const TimeStep& tp, bool restart=false) { S1.postSolve(tp,restart); }
  //! \brief Advances the time step \a n steps forward.
  void fastForward(int n) { for (int i = 0; i < n; i++) this->advanceStep(); }

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, DataExporter* exporter = NULL,
                           const char* heading = NULL)
  {
    int geoBlk = 0;
    int nBlock = 0;
    SolverHandler<T1> handler(S1,infile);

    if (!handler.pre(tp,geoBlk,nBlock))
      return 1;

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
    for (int iStep = 1; handler.post(tp) || this->advanceStep(); iStep++)
    {
      // Solve
      int res = handler.solve(tp,nBlock,exporter);
      if (res)
        return res;

      IFEM::pollControllerFifo();
    }

    return 0;
  }

  //! \brief Parses a data section from an input stream.
  virtual bool parse(char*, std::istream&) { return true; }

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
