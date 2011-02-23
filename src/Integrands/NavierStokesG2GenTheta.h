//==============================================================================
//!
//! \file NavierStokesG2GenTheta.h
//!
//! \date August 3 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for Navier-Stokes equations with G2 stabilization
//!        and generalized theta method for time integration
//!
//==============================================================================

#ifndef _NAVIER_STOKES_G2_GEN_THETA_H
#define _NAVIER_STOKES_G2_GEN_THETA_H

#include "NavierStokesG2.h"

/*!
  \brief Class representing the integrand of the pressure stabilized Navier-Stokes problem.
  \details This class supports G2-stabilized Navier-Stokes formulations using equal
  order elements for velocity and pressure. Time integration by the generalized theta method.
*/


class NavierStokesG2GenTheta : public NavierStokesG2
{
 public:

  //! \brief The default constructor initializes all pointers to zero
  //! \param[in] n Number of spatial dimensions
  NavierStokesG2GenTheta(short int n, SIM::Formulation form, int itg = 3);
  //! \brief The destructor clears the static work arrays used internally.
  virtual ~NavierStokesG2GenTheta();

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

 private:
  int   step;                              // Current step (1, 2 or 3)
  real  dt;                                // Current (sub-)timestep size
  //static const real theta  = 0.292893219;  // Time integration parameter
  //static const real thetam = 0.414213562;  // Time integration parameter
  //static const real thetat = 0.171572875;  // Time integration parameter
  //static const real alpha  = 0.585786438;  // Time integration parameter
  //static const real beta   = 0.414213562;  // Time integration parameter
  static const real theta  = 0.5;
  static const real thetam = 0.5;
  static const real thetat = 0.5;
  static const real alpha  = 1.0;
  static const real beta   = 1.0;
};

#endif
