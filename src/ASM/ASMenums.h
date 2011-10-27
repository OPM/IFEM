// $Id$
//==============================================================================
//!
//! \file ASMenums.h
//!
//! \date Oct 27 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various enums for assembly scope.
//!
//==============================================================================

#ifndef _ASM_ENUMS_H
#define _ASM_ENUMS_H


namespace ASM //! Assembly scope
{
  //! \brief Enum defining the available discretization methods.
  enum Discretization
  {
    Lagrange = 0,
    Spectral = 1,
    Spline   = 2,
    LRSpline = 3,
    SplineC1 = 4
  };
}

#endif
