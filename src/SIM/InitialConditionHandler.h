// $Id$
//==============================================================================
//!
//! \file InitialConditionHandler.h
//!
//! \date Oct 31 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Functions for loading initial conditions
//!
//==============================================================================

#ifndef INITIAL_CONDITION_HANDLER_H
#define INITIAL_CONDITION_HANDLER_H

#include "SIMbase.h"
#include "SIMdependency.h"

namespace SIM
{
  //! \brief Set initial conditions
  //! \param sim The SIM to load initial conditions for
  //! \param fieldHolder The SIM to inject the initial conditions into. NULL means use sim
  bool setInitialConditions(SIMbase& sim, SIMdependency* fieldHolder=NULL);
}

#endif
