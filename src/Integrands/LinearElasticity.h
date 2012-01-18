// $Id$
//==============================================================================
//!
//! \file LinearElasticity.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for linear elasticity problems.
//!
//==============================================================================

#ifndef _LINEAR_ELASTICITY_H
#define _LINEAR_ELASTICITY_H

#include "Elasticity.h"


/*!
  \brief Class representing the integrand of the linear elasticity problem.
  \details Most methods of this class are inherited form the base class.
  Only the \a evalInt and \a evalBou methods, which are specific for linear
  elasticity problems (and not used in nonlinear problems) are implemented here.
*/

class LinearElasticity : public Elasticity
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, and axisymmetric 3D formulation is assumed
  LinearElasticity(unsigned short int n = 3, bool axS = false);
  //! \brief Empty destructor.
  virtual ~LinearElasticity() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;
};

#endif
