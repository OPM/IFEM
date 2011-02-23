//==============================================================================
//!
//! \file StabilizedNavierStokes.h
//!
//! \date Mar 17 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for pressure stabilized Navier-Stokes problems.
//!
//==============================================================================

#ifndef _STABILIZED_NAVIER_STOKES_H
#define _STABILIZED_NAVIER_STOKES_H

#include "StabilizedStokes.h"

/*!
  \brief Class representing the integrand of the pressure stabilized Navier-Stokes problem.
  \details This class supports pressure stabilized Navier-Stokes formulations using equal
  order elements for velocity and pressure.
*/


class StabilizedNavierStokes : public StabilizedStokes
{
 public:

  //! \brief The default constructor initializes all pointers to zero
  //! \param[in] n Number of spatial dimensions
  StabilizedNavierStokes(short int n, SIM::Formulation form, int itg = 3);
  //! \brief The destructor clears the static work arrays used internally.
  virtual ~StabilizedNavierStokes();

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
  //! \param[in] dNdX   Basis function gradients
  //! \param[in] d2NdX2 Basis function second derivatives
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] h characteristic element length
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Matrix3D& d2NdX2, 
		       const Vec3& X, double h = 0.0) const;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] Navg Average value of basis function over element
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vector& Navg, const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns which integrand to use
  int getIntegrandType() const { return 3; }
};

#endif
