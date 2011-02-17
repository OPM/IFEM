//==============================================================================
//!
//! \file SIMparameters.h
//!
//! \date Feb 11 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for encapsulation generic simulation parameters.
//!
//==============================================================================

#ifndef SIMPARAMETERS_H_
#define SIMPARAMETERS_H_

#include "MatVec.h"
#include "TimeDomain.h"

class SIMparameters {
  public:
    //! \brief The constructor initializes the counters to zero.
    SIMparameters();

    //! \brief Empty destructor
    virtual ~SIMparameters() {};

    //! \brief Returns \e true if the simulation consists of several load steps.
    bool multiSteps() const;

    //! \brief Increments the time/load step.
    //! \return \e true, if we have reached the end of the simulation.
    bool increment();

    int  step;        //!< Load/time step counter
    int& iter;        //!< Iteration counter
    double startTime; //!< Start (pseudo)time of simulation
    double stopTime;  //!< Stop (pseudo)time of simulation
    RealArray tInc;   //!< Time (or pseudo time) increment size(s)
    TimeDomain time;  //!< Time domain data
};


#endif // SIMPARAMETERS_H_
