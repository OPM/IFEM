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
  //! \brief The constructor forwards to the parent class constructor.
  AdaptiveISolver(T1& sim, bool sa) : AdaptiveSIM(sim,sa), model(sim) {}

protected:
  //! \brief Assembles and solves the linearized FE equation system.
  bool assembleAndSolveSystem() override
  {
    TimeStep dummy;
    model.init(dummy);
    if (model.solveStep(dummy))
      solution = model.getSolutions();
    else
      return false;

    return true;
  }

  //! \brief Saves point results to output file for a given refinement step.
  //! \param[in] iStep Refinement step counter
  bool savePoints(int iStep) const override
  {
    return model.savePoints(0.0,iStep);
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
    this->S1.opt.saveNorms = true;
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverAdapImpl() {}

  //! \brief Reads solver data from the specified input file.
  virtual bool read(const char* file) { return this->SIMadmin::read(file); }

  //! \brief Solves the problem on a sequence of adaptively refined meshes.
  virtual int solveProblem(char* infile, const char* = nullptr)
  {
    if (!aSim.initAdaptor())
      return 1;

    if (SIMSolverStat<T1>::exporter)
      SIMSolverStat<T1>::exporter->setFieldValue(exporterName, &this->S1,
                                                 &aSim.getSolution(),
                                                 &aSim.getProjections(),
                                                 &aSim.getEnorm());

    for (int iStep = 1; aSim.adaptMesh(iStep); iStep++)
      if (!aSim.solveStep(infile,iStep))
        return 1;
      else if (!aSim.writeGlv(infile,iStep))
        return 2;
      else if (SIMSolverStat<T1>::exporter)
        SIMSolverStat<T1>::exporter->dumpTimeLevel(nullptr,true);

    return 0;
  }

  //! \brief Set name of data exporter registration to use.
  void setExporterName(const std::string& name) { exporterName = name; }

protected:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* kyw, std::istream& is) { return aSim.parse(kyw,is); }
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) { return aSim.parse(elem); }

  AdapSim aSim; //!< Adaptive simulation driver
  std::string exporterName = "u"; //!< Name for data exporter registration to use
};


//! Convenience alias template
template<class T1>
using SIMSolverAdap = SIMSolverAdapImpl<T1,AdaptiveSIM>;


/*!
  \brief Template class for adaptive simulator drivers
         which handle adaptation internally.
*/

template<class T1>
class SIMSolverAdapInternal : public SIMSolverStat<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMSolverAdapInternal(T1& s1) : SIMSolverStat<T1>(s1) {}

  //! \brief Empty destructor.
  virtual ~SIMSolverAdapInternal() {}

  //! \brief Reads solver data from the specified input file.
  bool read(const char* file) override { return this->SIMadmin::read(file); }

  //! \brief Solves the problem on a sequence of adaptively refined meshes.
  int solveProblem(char* infile, const char* = nullptr) override
  {
    if (!this->S1.initAdapPrm())
      return 1;

    int geoBlk = 0, nBlock = 0;
    for (int iStep = 1; this->S1.adaptMesh(iStep); iStep++) {
      IFEM::cout <<"\nAdaptive step "<< iStep << std::endl;
      if (!this->S1.solveStep(infile,iStep))
        return 1;
      else if (!this->S1.projectNorms(iStep))
        return 2;
      else if (!this->saveState(geoBlk,nBlock,iStep,infile))
        return 3;
    }

    return 0;
  }

  //! \brief Parse an element from an XML input file.
  bool parse(const TiXmlElement* elem) override
  {
    return this->S1.parse(elem);
  }

protected:
  //! \brief Saves geometry and results to VTF and HDF5 for current time step.
  bool saveState(int& geoBlk, int& nBlock, int iStep, char* infile)
  {
    if (!this->S1.saveModel(iStep == 1 ? infile : nullptr,geoBlk,nBlock))
      return false;

    TimeStep tp;
    tp.step = iStep;

    if (!this->S1.saveElmNorms(iStep,nBlock))
      return false;

    if (!this->S1.saveProjections(iStep,nBlock))
      return false;

    if (!this->S1.saveStep(tp,nBlock))
      return false;

    return this->exporter ? this->exporter->dumpTimeLevel(&tp,true) : true;
  }
};

#endif
