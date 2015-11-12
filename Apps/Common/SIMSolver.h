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
  SolverHandler(T& sim) : m_sim(sim) {}

  int pre(char* infile, int& nBlock, const TimeStep* tp = nullptr)
  {
    // Save FE model to VTF file for visualization
    int geoBlk = nBlock = 0;
    if (!m_sim.saveModel(infile,geoBlk,nBlock))
      return 1;

    // Optionally save the initial configuration
    return tp && !m_sim.saveStep(*tp,nBlock) ? 2 : 0;
  }

  //! \brief Solve for a single time step.
  //! \param tp Time stepping parameters.
  //! \param nBlock Running VTF block counter.
  int solve(TimeStep& tp, int& nBlock, DataExporter* exporter)
  {
    if (!m_sim.solveStep(tp))
      return 3;

    // Save solution
    if (!m_sim.saveStep(tp,nBlock))
      return 4;

    if (exporter)
      exporter->dumpTimeLevel(&tp);

    return 0;
  }

  //! \brief Post step for a single time step.
  bool post(TimeStep&) { return false; }

  T& m_sim; //!< Reference to simulator to solve for.
};


//! \brief Specialization for adaptive simulators.
template<>
struct SolverHandler<AdaptiveSIM> {
  //! \brief Constructor.
  SolverHandler<AdaptiveSIM>(AdaptiveSIM& sim) : m_sim(sim), m_infile(nullptr)
  {
    m_step = 0;
    m_norms = sim.getNoNorms();
  }

  int pre(char* infile, int&, const TimeStep* = nullptr)
  {
    m_infile = infile;
    return 0;
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
      exporter->dumpTimeLevel(nullptr, true);

    return 0;
  }

  //! \brief Post step for a single adaptive loop.
  bool post(TimeStep& tp)
  {
    if (++m_step == 1 || m_sim.adaptMesh(m_step))
      return true;

    m_step = 0;
    tp.step--; // advanceStep will increment it
    return false;
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

  //! \brief Returns a const reference to the time stepping information.
  const TimeStep& getTimePrm() const { return tp; }

  //! \brief Advances the time step one step forward.
  bool advanceStep() { return tp.increment() && S1.advanceStep(tp); }
  //! \brief Advances the time step \a n steps forward.
  void fastForward(int n) { for (int i = 0; i < n; i++) this->advanceStep(); }

  void postSolve(const TimeStep& tp, bool restart=false) { S1.postSolve(tp,restart); }

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, DataExporter* exporter = nullptr,
                           const char* heading = nullptr, bool saveInit = true)
  {
    int res, nBlock = 0;
    SolverHandler<T1> handler(S1);
    if (saveInit)
      res = handler.pre(infile,nBlock,&tp);
    else
      res = handler.pre(infile,nBlock);
    if (res) return res;

    if (heading)
    {
      // Write an application-specific heading, if provided
      std::string myHeading(heading);
      size_t n = myHeading.find_last_of('\n');
      if (n+1 < myHeading.size()) n = myHeading.size()-n;
      IFEM::cout <<"\n\n"<< myHeading <<"\n";
      for (size_t i = 0; i < n && i < myHeading.size(); i++) IFEM::cout <<'=';
      IFEM::cout << std::endl;
    }

    // Solve for each time step up to final time
    for (int iStep = 1; handler.post(tp) || this->advanceStep(); iStep++)
      if ((res = handler.solve(tp,nBlock,exporter)))
        return res;
      else
        IFEM::pollControllerFifo();

    return 0;
  }

protected:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyw, std::istream& is) { return tp.parse(keyw,is); }
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) { return tp.parse(elem); }

  TimeStep tp; //!< Time stepping information
  T1&      S1; //!< The actual solver
};

#endif
