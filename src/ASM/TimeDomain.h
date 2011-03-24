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
  \brief Struct for representing the time domain.
*/

struct TimeDomain
{
  double   t; //!< Current time (or pseudo time, load parameter)
  double  dt; //!< Current time (or load parameter) increment
  int     it; //!< Current iteration within current time/load step
  bool first; //!< If \e true, this is the first load/time step

  //! \brief Default constructor.
  TimeDomain(int i = 0, bool f = true) : it(i), first(f) { t = dt = 0.0; }
};

#endif
