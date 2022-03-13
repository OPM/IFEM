// $Id$
//==============================================================================
//!
//! \file SplineField.h
//!
//! \date Mar 15 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utility class for spline-based finite element fields.
//!
//==============================================================================

#ifndef _SPLINE_FIELD_H
#define _SPLINE_FIELD_H

#include "MatVec.h"

#include <vector>

namespace Go {
  class SplineSurface;
  class SplineVolume;
}
class ItgPoint;


/*!
  \brief Utility class for spline-based finite element fields.

  \details This class implements the method to extract the geometry mapping
            and evaluate a basis.
*/

class SplineField
{
public:
  //! \brief Evaluate the Jacobian (and Hessian) of the geometry mapping.
  //! \param surf Surface to evaluate for
  //! \param nsd Number of spatial dimensions
  //! \param x Point to evaluate in
  //! \param ip Vector of basis function indices
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Basis function derivatives in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalMapping(const Go::SplineSurface& surf,
                          int nsd,
                          const ItgPoint& x,
                          std::vector<int>& ip,
                          Matrix& Xnod,
                          Matrix& Jac,
                          Matrix& dNdX,
                          Matrix3D* d2NdX2 = nullptr,
                          Matrix3D* Hess = nullptr);

  //! \brief Evaluate a basis.
  //! \param surf Surface to evaluate for
  //! \param x Point to evaluate in
  //! \param ip Vector of basis function indices
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Derivatives of basis functions in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalBasis (const Go::SplineSurface& surf,
                         const ItgPoint& x,
                         std::vector<int>& ip,
                         const Matrix& Xnod,
                         const Matrix& Jac,
                         Matrix& dNdX,
                         Matrix3D* d2NdX2 = nullptr,
                         Matrix3D* Hess = nullptr);

  //! \brief Evaluate the Jacobian (and Hessian) of the geometry mapping.
  //! \param vol Volume to evaluate for
  //! \param x Point to evaluate in
  //! \param ip Vector of basis function indices
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Basis function derivatives in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalMapping(const Go::SplineVolume& vol,
                          const ItgPoint& x,
                          std::vector<int>& ip,
                          Matrix& Xnod,
                          Matrix& Jac,
                          Matrix& dNdX,
                          Matrix3D* d2NdX2 = nullptr,
                          Matrix3D* Hess = nullptr);

  //! \brief Evaluate a basis.
  //! \param vol Volume to evaluate for
  //! \param x Point to evaluate in
  //! \param ip Vector of basis function indices
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Derivatives of basis functions in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalBasis (const Go::SplineVolume& vol,
                         const ItgPoint& x,
                         std::vector<int>& ip,
                         const Matrix& Xnod,
                         const Matrix& Jac,
                         Matrix& dNdX,
                         Matrix3D* d2NdX2 = nullptr,
                         Matrix3D* Hess = nullptr);
};

#endif
