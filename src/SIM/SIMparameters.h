// $Id$
//==============================================================================
//!
//! \file SIMparameters.h
//!
//! \date Feb 11 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for encapsulation of general simulation parameters.
//!
//==============================================================================

#ifndef SIM_PARAMETERS_H_
#define SIM_PARAMETERS_H_

#include "TimeDomain.h"
#include <vector>


/*!
  \brief Class for encapsulation of general simulation parameters.
*/

class SIMparameters
{
public:
  //! \brief The constructor initializes the counters to zero.
  SIMparameters() : step(0), iter(time.it) {}

  //! \brief Empty destructor.
  virtual ~SIMparameters() {}

  //! \brief Returns \e true if the simulation consists of several load steps.
  bool multiSteps() const;

  //! \brief Increments the time/load step.
  //! \return \e true, if we have reached the end of the simulation.
  bool increment();

  int    step;      //!< Load/time step counter
  int&   iter;      //!< Iteration counter
  double startTime; //!< Start (pseudo)time of simulation
  double stopTime;  //!< Stop (pseudo)time of simulation

  std::vector<double> tInc; //!< Time (or pseudo time) increment size(s)
  TimeDomain          time; //!< Time domain data
};

#endif
