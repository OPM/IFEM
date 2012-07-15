// $Id$
//==============================================================================
//!
//! \file NonlinearDriver.h
//!
//! \date Jul 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear solution driver for finite deformation FEM analysis.
//!
//==============================================================================

#ifndef _NONLINEAR_DRIVER_H
#define _NONLINEAR_DRIVER_H

#include "NonLinSIM.h"


/*!
  \brief Nonlinear solution driver for finite deformation FEM analysis.
*/

class NonlinearDriver : public NonLinSIM
{
public:
  //! \brief Default constructor.
  //! \param sim Pointer to the spline FE model
  NonlinearDriver(SIMbase* sim = 0) : NonLinSIM(sim) {}
  //! \brief Empty destructor.
  virtual ~NonlinearDriver() {}

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(const TimeDomain& time, bool energyNorm,
			     double zero_tolerance = 1.0e-8,
			     std::streamsize outPrec = 0);
};

#endif
