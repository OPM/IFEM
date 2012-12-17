//==============================================================================
//!
//! \file CFDenums.h
//!
//! \date Nov 9 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various enums for CFD simulators.
//!
//==============================================================================

#ifndef _CFD_ENUMS_H
#define _CFD_ENUMS_H


namespace CFD //! CFD scope
{
  //! \brief Enum defining various solution formulations that may occur.
  //! \details Used as a bit field - power of two values only
  enum Formulation
  {
    NONE       = 0,
    LAPLACE    = 1,
    STRESS     = 2,
    RANS       = 4,
    BOUSSINESQ = 8,
    ALE        = 16,
    SUPG       = 32,
    VMS        = 64
  };
}

#endif
