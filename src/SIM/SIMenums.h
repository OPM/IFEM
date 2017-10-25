// $Id$
//==============================================================================
//!
//! \file SIMenums.h
//!
//! \date Jul 8 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various enums for simulation scope.
//!
//==============================================================================

#ifndef _SIM_ENUMS_H
#define _SIM_ENUMS_H


namespace SIM //! Simulation scope
{
  //! \brief Enum defining various solution formulations that may occur.
  //! \details Used as a bit field - power of two values only.
  enum Formulation
  {
    NONE      = 0,
    LINEAR    = 1,
    NONLINEAR = 2
  };

  //! \brief Enum defining the various solution modes that may occur.
  enum SolutionMode
  {
    INIT = 0,
    STATIC = 1,
    ARCLEN = 2,
    DYNAMIC,
    VIBRATION,
    BUCKLING,
    STIFF_ONLY,
    MASS_ONLY,
    RHS_ONLY,
    INT_FORCES,
    RECOVERY,
    NORMS
  };

  //! \brief Enum defining the various convergence statuses that may occur.
  //! \note The order of the values are not random. It reflects the severity of
  //! the status somehow (the smaller the value, the further from convergence).
  enum ConvStatus
  {
    FAILURE,
    DIVERGED,
    SLOW,
    OK,
    CONVERGED
  };
}

#endif
