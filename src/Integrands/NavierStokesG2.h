//==============================================================================
//!
//! \file NavierStokesG2.h
//!
//! \date May 28 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for Navier-Stokes equations with G2 stabilization.
//!
//==============================================================================

#ifndef _NAVIER_STOKES_G2_H
#define _NAVIER_STOKES_G2_H

#include "StabilizedStokes.h"

/*!
  \brief Class representing the integrand of the pressure stabilized Navier-Stokes problem.
  \details This class supports G2-stabilized Navier-Stokes formulations using equal
  order elements for velocity and pressure.
*/


class NavierStokesG2 : public StabilizedStokes
{
 public:

  //! \brief The default constructor initializes all pointers to zero
  //! \param[in] n Number of spatial dimensions
  NavierStokesG2(short int n, ProblemFormulation form, int itg = 3);
  //! \brief The destructor clears the static work arrays used internally.
  virtual ~NavierStokesG2();

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used when \a getIntegrandType returns 1.
  virtual bool evalInt(LocalIntegral*& elmInt,
                       const TimeDomain& time, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X) const
  { return false; }

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
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt,
                       const TimeDomain& time, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Matrix3D& d2NdX2, const Vec3& X,
                       double h = 0.0) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] Navg Volume-averaged basis function values over the element
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used when \a getIntegrandType returns 3.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt,
                       const TimeDomain& time, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vector& Navg, const Vec3& X) const
  { return true; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X) const
  { return false; }
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
                       const Vector& Navg, const Vec3& X) const
  { return false; }

  //! \brief Evaluates the secondary solution at current integration point.
  //! \param[out] s The solution field values
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] MNPC Matrix of nodal point correspondance
  virtual bool evalSol(Vector& s,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Returns which integrand to use
  int getIntegrandType() const { return 2; }
};

#endif
