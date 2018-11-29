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
    Triangle =-1,
    Lagrange = 0,
    Spectral = 1,
    // The spline entries need to be at the end and successively numbered
    Spline   = 2,
    SplineC1 = 3,
    LRSpline = 4,
    LRNurbs  = 5
  };

  //! \brief Operations to be applied after summing norm element contributions.
  enum FinalNormOp
  {
    NONE = 0,
    ABS  = 1,
    SQRT = 2
  };
}

#endif
