// $Id$
//==============================================================================
//!
//! \file TimeDomain.h
//!
//! \date May 27 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Time domain representation for time-dependent and nonlinear solvers.
//!
//==============================================================================

#ifndef _TIME_DOMAIN_H
#define _TIME_DOMAIN_H


/*!
  \brief Struct representing the time domain.
*/

struct TimeDomain
{
  double   t; //!< Current time (or pseudo time, load parameter)
  double  dt; //!< Current timestep (or load parameter) increment
  double dtn; //!< Previous timestep (or load parameter) increment
  double CFL; //!< Current CFL number (used by CFD simulators)
  double incNorm; //!< Norm of change between timesteps
  int     it; //!< Current iteration within current time/load step
  char first; //!< If \e true, this is the first load/time step

  //! \brief Default constructor.
  explicit TimeDomain(int i = 0, bool f = true) : it(i), first(f)
  { t = dt = dtn = CFL = incNorm = 0.0; }
  //! \brief Constructor for linear problems with fixed (non-zero) time.
  explicit TimeDomain(double t0) : t(t0), it(0), first(true)
  { dt = dtn = CFL = incNorm = 0.0; }
};

#endif
