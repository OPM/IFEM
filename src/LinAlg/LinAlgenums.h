// $Id$
//==============================================================================
//!
//! \file LinAlgenums.h
//!
//! \date Oct 23 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various enums for linear algebra scope.
//!
//==============================================================================

#ifndef _LINALG_ENUMS_H
#define _LINALG_ENUMS_H


namespace LinAlg //! Linear algebra scope
{
  //! \brief Enum defining linear system properties.
  enum LinearSystemType 
  {
    GENERAL_MATRIX = 0, //!< General matrix
    SYMMETRIC      = 1, //!< Symmetric matrix (value and structure)
    SPD            = 2  //!< Symmetric, positive definite matrix
  };
}

#endif
