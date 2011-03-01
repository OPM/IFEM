// $Id$
//==============================================================================
//!
//! \file SIMparameters.C
//!
//! \date Feb 11 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for encapsulation of general simulation parameters.
//!
//==============================================================================

#include "SIMparameters.h"


bool SIMparameters::multiSteps () const
{
  if (tInc.empty()) return false;

  const double epsT = 1.0e-6;
  return (startTime+(1.0+epsT)*tInc.front() < stopTime);
}


bool SIMparameters::increment ()
{
  const double epsT = 1.0e-6;

  if (++step <= (int)tInc.size() && step > 0)
    time.dt = tInc[step-1];

  if (time.t+time.dt*epsT >= stopTime)
    return false;

  time.t += time.dt;

  if (time.t > stopTime)
  {
    // Adjust the size of the last time step
    time.dt += stopTime - time.t;
    time.t = stopTime;
  }

  return true;
}
