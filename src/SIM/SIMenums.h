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
  //! \details Used as a bit field - power of two values only
  enum Formulation
  {
    NONE      = 0,
    LINEAR    = 1,
    NONLINEAR = 2,
    LAPLACE   = 4,
    STRESS    = 8,
    RANS      = 16
  };

  //! \brief Enum defining the various solution modes that may occur.
  enum SolutionMode
  {
    INIT,
    STATIC,
    DYNAMIC,
    VIBRATION,
    BUCKLING,
    STIFF_ONLY,
    MASS_ONLY,
    RHS_ONLY,
    RECOVERY
  };
}

#endif
