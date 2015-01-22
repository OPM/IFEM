// $Id$
//==============================================================================
//!
//! \file EigenModeSIM.h
//!
//! \date Jan 8 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for computing a time history from a set of eigen modes.
//!
//==============================================================================

#ifndef _EIGEN_MODE_SIM_H
#define _EIGEN_MODE_SIM_H

#include "MultiStepSIM.h"
#include "SIMbase.h"


/*!
  \brief Driver for computing a time history from a set of eigen mode shapes.
*/

class EigenModeSIM : public MultiStepSIM
{
public:
  //! \brief The constructor initializes the FE model reference.
  EigenModeSIM(SIMbase& sim);
  //! \brief Empty destructor.
  virtual ~EigenModeSIM() {}

  using MultiStepSIM::parse;
  //! \brief Parses a data section from an XML document.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out problem-specific data to the given stream.
  virtual void printProblem(std::ostream& os) const;

  //! \brief Initializes primary solution vectors.
  virtual bool initSol(size_t nSol);

  //! \brief Computes the primary solution at current time step.
  //! \param param Time stepping parameters
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual SIM::ConvStatus solveStep(TimeStep& param, SIM::SolutionMode,
                                    double zero_tolerance,
                                    std::streamsize outPre);

private:
  double              myStart;   //!< The start time of this simulator
  std::vector<Mode>   modes;     //!< Eigenmode shapes and frequencies
  std::vector<double> amplitude; //!< Eigenmode shape amplitudes
};

#endif
