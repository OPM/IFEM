// $Id$
//==============================================================================
//!
//! \file ForceIntegrator.h
//!
//! \date May 10 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for integration of boundary forces.
//!
//==============================================================================

#ifndef _FORCE_INTEGRATOR_H
#define _FORCE_INTEGRATOR_H

#include "MatVec.h"

class SIMbase;
class TimeDomain;
class Vec3;


namespace SIM
{
  //! \brief Integrate the force resultant on a specified boundary.
  //! \param[in] solution Primary solution vectors in DOF order
  //! \param[in] model The isogeometric finite element model
  //! \param[in] code Code indentifying the boundary subjected to integration
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X0 Pivot point for torque calculation
  Vector getBoundaryForce(const Vectors& solution, SIMbase* model, int code,
			  const TimeDomain& time, const Vec3* X0 = NULL);
}

#endif
