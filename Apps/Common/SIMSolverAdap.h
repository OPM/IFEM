// $Id$
//==============================================================================
//!
//! \file SIMSolverAdap.h
//!
//! \date Feb 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Stationary, adaptive SIM solver class template.
//!
//==============================================================================

#ifndef _SIM_SOLVER_ADAP_H_
#define _SIM_SOLVER_ADAP_H_

#include "SIMSolver.h"
#include "AdaptiveSIM.h"
#include "TimeStep.h"


/*!
  \brief Adaptive simulator driver using the ISolver interface.
*/

template<class T1>
class AdaptiveISolver : public AdaptiveSIM
{
public:
  //! \brief The constructor initializes default adaptation parameters.
  //! \param sim The FE model
  //! \param sa If \e true, this is a stand-alone driver
  AdaptiveISolver(T1& sim, bool sa = true) : AdaptiveSIM(sim, sa), model(sim) {}

protected:
  //! \brief Assemble and solve the equation system.
  virtual bool assembleAndSolveSystem()
  {
    TimeStep dummy;
    model.init(dummy);
    bool result = model.solveStep(dummy);
    if (result)
      solution = model.getSolutions();

    rCond = 1.0;

    return result;
  }

  T1& model; //!< Reference to the actual sim
};


/*!
  \brief Template class for stationary adaptive simulator drivers.
  \details This template can be instantiated over any type implementing the
  ISolver interface. It provides an adaptive loop with data output.
*/

template<class T1, class AdapSim>
class SIMSolverAdapImpl : public SIMSolverStat<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMSolverAdapImpl(T1& s1) : SIMSolverStat<T1>(s1), aSim(s1,false)
  {
    this->S1.setSol(&aSim.getSolution());
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverAdapImpl() {}

  //! \brief Reads solver data from the specified input file.
  virtual bool read(const char* file) { return this->SIMadmin::read(file); }

  //! \brief Solves the problem on a sequence of adaptively refined meshes.
  virtual int solveProblem(char* infile, const char* = nullptr, bool = false)
  {
    if (!aSim.initAdaptor())
      return 1;

    for (int iStep = 1; aSim.adaptMesh(iStep); iStep++)
      if (!aSim.solveStep(infile,iStep))
        return 1;
      else if (!aSim.writeGlv(infile,iStep))
        return 2;
      else if (SIMSolverStat<T1>::exporter)
        SIMSolverStat<T1>::exporter->dumpTimeLevel(nullptr,true);

    return 0;
  }

protected:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyw, std::istream& is)
  {
    return aSim.parse(keyw,is);
  }
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    return aSim.parse(elem);
  }

  AdapSim aSim; //!< Adaptive simulation driver
};

template<class T1>
using SIMSolverAdap = SIMSolverAdapImpl<T1, AdaptiveSIM>; //!< Convenience alias template

#endif
