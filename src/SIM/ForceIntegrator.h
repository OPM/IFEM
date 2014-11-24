// $Id$
//==============================================================================
//!
//! \file ForceIntegrator.h
//!
//! \date May 10 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for integration of boundary and nodal forces.
//!
//==============================================================================

#ifndef _FORCE_INTEGRATOR_H
#define _FORCE_INTEGRATOR_H

#include "MatVec.h"

class SIMbase;
struct TimeDomain;
class Vec3;
class ForceBase;
class GlbForceVec;


namespace SIM
{
  //! \brief Integrates the force resultant on a specified boundary.
  //! \param[in] solution Primary solution vectors in DOF order
  //! \param[in] model The isogeometric finite element model
  //! \param[in] code Code indentifying the boundary subjected to integration
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X0 Pivot point for torque calculation
  //! \return The force (and torque) resultant
  Vector getBoundaryForce(const Vectors& solution, SIMbase* model, int code,
                          const TimeDomain& time, const Vector* X0);
  //! \brief Integrates the force resultant on a specified boundary.
  //! \param[in] solution Primary solution vectors in DOF order
  //! \param[in] model The isogeometric finite element model
  //! \param[in] code Code indentifying the boundary subjected to integration
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X0 Pivot point for torque calculation
  //! \return The force (and torque) resultant
  Vector getBoundaryForce(const Vectors& solution, SIMbase* model, int code,
                          const TimeDomain& time, const Vec3* X0 = NULL);

  //! \brief Integrates nodal forces on a specified boundary.
  //! \param[in] solution Primary solution vectors in DOF order
  //! \param[in] model The isogeometric finite element model
  //! \param[in] code Code indentifying the boundary subjected to integration
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param force Global nodal force container (compressed storage)
  bool getNodalForces(const Vectors& solution, SIMbase* model, int code,
                      const TimeDomain& time, GlbForceVec& force);

  //! \brief Detects the global nodes that reside on a specified boundary.
  //! \param[in] model The isogeometric finite element model
  //! \param[in] code Property code associated with the boundary
  //! \param[out] force Global nodal force container (compressed storage)
  bool initBoundaryNodeMap(SIMbase* model, int code, GlbForceVec& force);

  //! \brief Integrates a force integrand on a specified topology set.
  //! \param[in] solution Primary solution vectors in DOF order
  //! \param[in] model The isogeometric finite element model
  //! \param[in] code Code indentifying the boundary subjected to integration
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] forceInt The force integrand to integrate
  //! \param[in] force If non-NULL nodal forces are stored here
  //! \return True if integration succeeded
  bool integrate(const Vectors& solution, SIMbase* model, int code,
                 const TimeDomain& time, ForceBase* forceInt,
                 GlbForceVec* force=NULL);
}

#endif
