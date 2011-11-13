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


void SIMparameters::initTime (double start, double stop, const TimeSteps& steps)
{
  starTime = start;
  stopTime = stop;
  mySteps  = steps;
  time.t   = start;
  time.dt  = steps.empty() ? stop-start : steps.front().first.front();
  stepIt   = mySteps.begin();
}


static const double epsT = 1.0e-6; //!< Tolerance parameter


bool SIMparameters::multiSteps () const
{
  if (mySteps.empty()) return false;

  return (starTime+(1.0+epsT)*mySteps.front().first.front() < stopTime);
}


bool SIMparameters::increment ()
{
  if (stepIt != mySteps.end())
    if (++lstep <= stepIt->first.size())
      time.dt = stepIt->first[lstep-1];

  if (time.t+time.dt*epsT >= stopTime)
    return false; // We've reached the end of the simulation

  ++step;
  time.t += time.dt;

  if (stepIt != mySteps.end())
  {
    if (stepIt->first.size() <= lstep)
      stepIt->first.push_back(time.dt);
    if (time.t+time.dt*epsT >= stepIt->second)
    {
      if (time.t != stepIt->second)
      {
	// Adjust the size of the last time step
	time.dt += stepIt->second - time.t;
	time.t = stepIt->second;
	stepIt->first.back() = time.dt;
      }
      lstep = 0;
      ++stepIt;
    }
  }

  if (time.t > stopTime)
  {
    // Adjust the size of the last time step
    time.dt += stopTime - time.t;
    time.t = stopTime;
  }

  return true;
}
