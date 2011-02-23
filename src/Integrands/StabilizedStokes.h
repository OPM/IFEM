// $Id: StabilizedStokes.h,v 1.4 2010-12-06 08:49:59 rho Exp $
//==============================================================================
//!
//! \file StabilizedStokes.h
//!
//! \date Feb 11 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for pressure stabilized Stokes problems.
//!
//==============================================================================

#ifndef _STABILIZED_STOKES_H
#define _STABILIZED_STOKES_H

#include "Stokes.h"

/*!
  \brief Class representing the integrand of the pressure-stabilized
  Stokes problem.
  \details This class supports pressure stabilized Stokes formulations
  using equal order elements for velocity and pressure.
*/

class StabilizedStokes : public Stokes
{
public:

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] form The solution formulation to use
  //! \param[in] itg The integrandtype to use
  StabilizedStokes(short int n, SIM::Formulation form = SIM::LAPLACE, 
		   int itg = 3);
  //! \brief The destructor frees dynamically allocated objects.
  virtual ~StabilizedStokes() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X) const;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] d2NdX2 Basis function second derivatives
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] h Characteristic element length
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Matrix3D& d2NdX2,
		       const Vec3& X, double h = 0.0) const;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] Navg Volume-averaged basis function values over the element
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vector& Navg, const Vec3& X) const;
};

#endif
