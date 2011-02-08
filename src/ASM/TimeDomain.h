// $Id: TimeDomain.h,v 1.1 2010-06-01 09:40:25 kmo Exp $
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
  double t;  //!< Current time (or pseudo time)
  double dt; //!< Current time increment (or load parameter increment)
  int    it; //!< Current iteration within current time/load step

  //! \brief Default constructor.
  TimeDomain() { t = dt = 0.0; it = 0; }
};

#endif
