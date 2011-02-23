//==============================================================================
//!
//! \file ChorinVelPredBDF2.h
//!
//! \date Nov 23 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief BDF2 implementation for velocity prediction in Chorin's method
//!
//==============================================================================

#ifndef _CHORIN_VEL_PRED_BDF2_H
#define _CHORIN_VEL_PRED_BDF2_H

#include "ChorinVelPredBE.h"

/*!
  \brief Class representing the integrand of the velocity prediction step
  in Chorin's method.
  \details This class implements the velocity prediction step of Chorin's
  method using NURBS based FEM with equal order elements for velocity and
  pressure. Semi-implicit BDF2 time discretization.
*/

class ChorinVelPredBDF2 : public ChorinVelPredBE
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] form The solution formulation to use
  //! \param[in] itg The integrandtype to use
  ChorinVelPredBDF2(short int n, SIM::Formulation form, 
		    int itg, bool mixed = false);
  //! \brief The destructor frees dynamically allocated objects.
  virtual ~ChorinVelPredBDF2();

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used when \a getIntegrandType returns 1.
  bool evalInt(LocalIntegral*& elmInt,
	       const TimeDomain& time, double detJW,
	       const Vector& N, const Matrix& dNdX,
	       const Vec3& X) const;
  
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] d2NdX2 Basis function second derivatives
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] h Characteristic element length
  //!
  //! \details This interface is used when \a getIntegrandType returns 2.
  bool evalInt(LocalIntegral*& elmInt,
	       const TimeDomain& time, double detJW,
	       const Vector& N, const Matrix& dNdX,
	       const Matrix3D& d2NdX2, const Vec3& X,
	       double h = 0.0) const;

  //! \brief Evaluates the mixed field problem integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N1 Basis function values, field velocity
  //! \param[in] N2 Basis function values, field pressure
  //! \param[in] dN1dX Basis function gradients, velocity
  //! \param[in] dN2dX Basis function gradients, pressure
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral*& elmInt, 
	       const TimeDomain& time, double detJW,
	       const Vector& N1, const Vector& N2,
	       const Matrix& dN1dX, const Matrix& dN2dX,
	       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral*& elmInt, const TimeDomain& time,
	       double detJW, const Vector& N, const Matrix& dNdX,
	       const Vec3& X, const Vec3& normal) const;
};
#endif
