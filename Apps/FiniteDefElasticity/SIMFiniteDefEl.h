// $Id: SIMFiniteDefEl.h,v 1.2 2010-12-29 18:51:53 kmo Exp $
//==============================================================================
//!
//! \file SIMFiniteDefEl.h
//!
//! \date Dec 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for NURBS-based finite deformation analysis.
//!
//==============================================================================

#ifndef _SIM_FINITE_DEF_EL_H
#define _SIM_FINITE_DEF_EL_H

#include "SIMLinEl2D.h"
#include "SIMLinEl3D.h"


namespace SIM
{
  //! \brief Enum defining various finite deformation formulations.
  enum NlFormulation
  {
    TOTAL_LAGRANGE   = 3,
    NEOHOOKE         = 4,
    NEOHOOKE_IV      = 5,
    UPDATED_LAGRANGE = 6,
    MIXED_QnPn1      = 7,
    MIXED_QnQn1      = 8
  };
};


/*!
  \brief Driver class for 2D isogeometric finite deformation analysis.
*/

class SIMFiniteDefEl2D : public SIMLinEl2D
{
public:
  //! \brief Default constructor.
  //! \param[in] form Problem formulation option
  //! \param[in] planeStress Plane stress/plane strain option
  //! \param[in] options Additional solution formulation options
  SIMFiniteDefEl2D(int form = SIM::TOTAL_LAGRANGE, bool planeStress = false,
		   const std::vector<int>& options = std::vector<int>());
  //! \brief Empty destructor.
  virtual ~SIMFiniteDefEl2D() {}
};


/*!
  \brief Driver class for 3D isogeometric finite deformation analysis.
*/

class SIMFiniteDefEl3D : public SIMLinEl3D
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  //! \param[in] form Problem formulation option
  //! \param[in] options Additional solution formulation options
  SIMFiniteDefEl3D(bool checkRHS = false, int form = SIM::TOTAL_LAGRANGE,
		   const std::vector<int>& options = std::vector<int>());
  //! \brief Empty destructor.
  virtual ~SIMFiniteDefEl3D() {}
};

#endif
