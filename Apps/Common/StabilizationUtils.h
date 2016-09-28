#ifndef STABILIZATIONUTILS_H_
#define STABILIZATIONUTILS_H_

#include <vector>
#include "Vec3.h"
#include "MatVec.h"

namespace StabilizationUtils {
  //! \brief Returns characteristic element size
  //! \param XC The element corner coordinates
  //! \details The size is taken as the shortest edge length
  double getElementSize(const std::vector<Vec3>& XC, int nsd);

  //! \brief Returns stabilization parameters for convecction-diffusion equation
  //! \param[in] dt The timestep size
  //! \param[in] mu Diffusion/viscosity parameter
  //! \param[in] U  Velocity vector
  //! \param[in] G  G matrix
  //! \param[in] Ct VMS parameter
  //! \param[in] Cl VMS parameter
  //! \details Stabilization parameter in integration point
  double getTauPt(double dt, double mu, const Vector& U, const Matrix& G,
      const double Ct=2.0, const double Cl=36.0);
 		  
  //! \brief Computes stabilization parameters for Navier-Stokes equations
  //! \param[in] dt The timestep size
  //! \param[in] mu Diffusion/viscosity parameter
  //! \param[in] U  Velocity vector
  //! \param[in] G  The G matrix
  //! \param[out] tauM Stabilization parameter for momentum
  //! \param[out] tauC Stabilization parameter for continuity
  //! \param[in] Ct VMS parameter
  //! \param[in] Cl VMS parameter
  //! \details Stabilization parameters in integration point
  bool getTauNSPt(double dt, double mu, const Vector& U, const Matrix& G,
		  double& tauM, double& tauC, const double Ct=2.0, const double Cl=36.0);

  //! \brief Computes stabilization parameters for Navier-Stokes equations
  //! \param[in] dt The timestep size
  //! \param[in] mu Diffusion/viscosity parameter
  //! \param[in] U  Velocity vector
  //! \param[in] G  The G matrix
  //! \param[out] tauM Stabilization parameter for momentum
  //! \param[out] tauC Stabilization parameter for continuity
  //! \details Stabilization parameters in integration point
  bool getTauNSALEPt(double dt, double mu, const Vector& U, const Matrix& G,
		     double& tauM, double& tauC, const double Ct=2.0, const double Cl=36.0);

  //! \brief Computes variation of stability parameters with respect to velocity
  //! \param[in] U  Velocity vector
  //! \param[in] G  The G matrix
  //! \param[in] tauM Stabilization parameter for momentum
  //! \param[out] tauMjac Variation of tauM with respect to U
  bool getTauPtJac(const Vector& U, const Matrix& G, const double tauM, Vector& tauMjac);
  
  
  //! \brief Computes variation of stability parameters with respect to velocity
  //! \param[in] U  Velocity vector
  //! \param[in] G  The G matrix
  //! \param[in] tauM Stabilization parameter for momentum
  //! \param[in] tauC Stabilization parameter for continuity
  //! \param[out] tauMjac Variation of tauM with respect to U
  //! \param[out] tauCjac Variation of tauC with respect to U
  bool getTauNSPtJac(const Vector& U, const Matrix& G, const double tauM,
      const double& tauC, Vector& tauMjac, Vector& tauCjac);
}

#endif
