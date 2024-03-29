//==============================================================================
//!
//! \file StabilizationUtils.h
//!
//! \date Oct 31 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various helpers for stabilized formulations.
//!
//==============================================================================

#ifndef STABILIZATIONUTILS_H_
#define STABILIZATIONUTILS_H_

#include <vector>
#include "Vec3.h"
#include "MatVec.h"

namespace StabilizationUtils {
  //! \brief Returns characteristic element size.
  //! \param XC The element corner coordinates
  //! \param nsd Number of spatial dimensions
  //! \details The size is taken as the shortest edge length
  double getElementSize(const std::vector<Vec3>& XC, int nsd);

  //! \brief Returns stabilization parameters for convection-diffusion equation.
  //! \param[in] dt The timestep size
  //! \param[in] mu Diffusion/viscosity parameter
  //! \param[in] U  Velocity vector
  //! \param[in] G  G matrix
  //! \param[in] Ct VMS parameter
  //! \param[in] Cl VMS parameter
  //! \details Stabilization parameter in integration point
  double getTauPt(double dt, double mu, const Vector& U, const Matrix& G,
                  const double Ct = 2.0, const double Cl = 36.0);

  //! \brief Computes stabilization parameters for Navier-Stokes equations.
  //! \param[in] dt The timestep size
  //! \param[in] mu Diffusion/viscosity parameter
  //! \param[in] U  Velocity vector
  //! \param[in] G  The G matrix
  //! \param[in] Ct VMS parameter
  //! \param[in] Cl VMS parameter
  //! \return Stabilization parameters in integration point
  std::pair<double,double>
  getTauNSPt(double dt, double mu, const Vector& U, const Matrix& G,
             const double Ct = 2.0, const double Cl = 36.0);

  //! \brief Computes stabilization parameters for Navier-Stokes equations.
  //! \param[in] dt The timestep size
  //! \param[in] mu Diffusion/viscosity parameter
  //! \param[in] U Velocity vector
  //! \param[in] G The G matrix
  //! \param[in] Ct VMS parameter
  //! \param[in] Cl VMS parameter
  //! \return Stabilization parameters in integration point
  std::pair<double,double>
  getTauNSALEPt(double dt, double mu, const Vector& U, const Matrix& G,
                const double Ct = 2.0, const double Cl = 36.0);

  //! \brief Computes variation of stability parameters with respect to velocity.
  //! \param[in] U  Velocity vector
  //! \param[in] G  The G matrix
  //! \param[in] tauM Stabilization parameter for momentum
  Vector getTauPtJac(const Vector& U, const Matrix& G,
                   const double tauM);

  //! \brief Computes variation of stability parameters with respect to velocity.
  //! \param[in] U  Velocity vector
  //! \param[in] G  The G matrix
  //! \param[in] tauM Stabilization parameter for momentum
  //! \param[in] tauC Stabilization parameter for continuity
  std::pair<Vector,Vector>
  getTauNSPtJac(const Vector& U, const Matrix& G,
                const double tauM, const double& tauC);
}

#endif
