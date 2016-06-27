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


/*!
  \brief Template class for stationary adaptive simulator drivers.
  \details This template can be instanciated over any type implementing the
  ISolver interface. It provides an adaptive loop with data output.
*/

template<class T1> class SIMSolverAdap : public SIMSolver<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMSolverAdap(T1& s1) : SIMSolver<T1>(s1), aSim(s1,false) {}
  //! \brief Empty destructor.
  virtual ~SIMSolverAdap() {}

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, DataExporter* exporter = nullptr,
                           const char* heading = nullptr, bool = false)
  {
    if (exporter)
      exporter->setNormPrefixes(aSim.getNormPrefixes());

    this->S1.setSol(&aSim.getSolution());
    aSim.setupProjections();
    aSim.initAdaptor(0,2);

    this->printHeading(heading);

    for (int iStep = 1; aSim.adaptMesh(iStep); iStep++)
      if (!aSim.solveStep(infile,iStep))
        return 1;
      else if (!aSim.writeGlv(infile,iStep,aSim.getNoNorms()))
        return 2;
      else if (exporter)
        exporter->dumpTimeLevel(nullptr,true);

    return 0;
  }

protected:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyw, std::istream& is)
  {
    return this->SIMSolver<T1>::parse(keyw,is) && aSim.parse(keyw,is);
  }
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    return this->SIMSolver<T1>::parse(elem) && aSim.parse(elem);
  }

private:
  AdaptiveSIM aSim; //!< Adaptive simulation driver
};

#endif
