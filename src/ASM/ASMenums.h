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
    SuperElm =-2,
    Triangle =-1,
    Lagrange = 0,
    Spectral = 1,
    // The spline entries need to be at the end and successively numbered
    Spline   = 2,
    SplineC1 = 3,
    LRSpline = 4
  };

  //! \brief Operations to be applied after summing norm element contributions.
  enum FinalNormOp
  {
    NONE = 0,
    ABS  = 1,
    SQRT = 2
  };

  //! \brief Enum defining available basis function cache policies.
  enum CachePolicy {
    NO_CACHE,   //!< Cache is disabled - calculate on the fly
    PRE_CACHE,  //!< Cache basis function values up front, clear on assembly end
    ON_THE_FLY, //!< Cache basis functions on the fly
    FULL_CACHE  //!< Cache basis function values up front
  };

  extern CachePolicy cachePolicy; //!< Chosen basis function cache policy

  //! \brief Enumeration of different basis types.
  //! \details Entries should have non-positive values
  enum BasisType {
    GEOMETRY_BASIS     =  0, //!< Geometry basis
    PROJECTION_BASIS   = -1, //!< Projection basis
    PROJECTION_BASIS_2 = -2, //!< Second projection basis
    REFINEMENT_BASIS   = -3, //!< Refinement basis
    INTEGRATION_BASIS  = -4, //!< Integration basis
  };

  //! \brief Enum for cathegorization of result quantities.
  enum ResultClass {
    ANY       = 0,
    PRIMARY   = 1,
    SECONDARY = 2,
    PROJECTED = 3,
    ANALYTIC  = 4
  };
}

#endif
