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
#include <cstddef>


/*!
  \brief Class for encapsulation of general simulation parameters.
*/

class SIMparameters
{
public:
  //! \brief The constructor initializes the counters to zero.
  SIMparameters() : step(0), iter(time.it), lstep(0) { stepIt = mySteps.end(); }

  //! \brief Empty destructor.
  virtual ~SIMparameters() {}

  //! \brief Time stepping definition container.
  typedef std::vector< std::pair<std::vector<double>,double> > TimeSteps;

  //! \brief Initializes the time step definitions.
  void initTime(double start, double stop, const TimeSteps& steps);

  //! \brief Returns \e true if the simulation consists of several load steps.
  bool multiSteps() const;

  //! \brief Increments the time/load step.
  //! \return \e true, if we have reached the end of the simulation.
  bool increment();

  int        step; //!< Load/time step counter
  int&       iter; //!< Iteration counter
  TimeDomain time; //!< Time domain data

  double    starTime; //!< Start (pseudo)time of simulation
  double    stopTime; //!< Stop (pseudo)time of simulation
  TimeSteps mySteps;  //!< Time step definitions

private:
  size_t lstep; //!< Local step counter, i.e., within current \a *stepIt

  TimeSteps::iterator stepIt;  //!< Running iterator over the time steps
};

#endif
