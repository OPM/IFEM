// $Id$
//==============================================================================
//!
//! \file LRSplineField.h
//!
//! \date Mar 15 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utility class for LR spline-based finite element fields.
//!
//==============================================================================

#ifndef _LRSPLINE_FIELD_H
#define _LRSPLINE_FIELD_H

#include "MatVec.h"


namespace LR {
  class Element;
  class LRSplineSurface;
  class LRSplineVolume;
}
class ItgPoint;


/*!
  \brief Utility class for LR spline-based finite element fields.

  \details This class implements the method to extract the geometry mapping
            and evaluate a basis.
*/

class LRSplineField
{
public:
  //! \brief Evaluate the jacobian (and hessian) of the geometry mapping.
  //! \param surf Surface to evaluate
  //! \param x Point to evaluate in
  //! \param elm Element info for point
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Derivatives of basis functions in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalMapping(const LR::LRSplineSurface& surf,
                          const ItgPoint& x,
                          const LR::Element*& elm,
                          Matrix& Xnod,
                          Matrix& Jac,
                          Matrix& dNdX,
                          Matrix3D* d2NdX2 = nullptr,
                          Matrix3D* Hess = nullptr);

  //! \brief Evaluate a basis.
  //! \param surf Surface to evaluate
  //! \param x Point to evaluate in
  //! \param elm Element info for point
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Derivatives of basis functions in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalBasis(const LR::LRSplineSurface& surf,
                        const ItgPoint& x,
                        const LR::Element*& elm,
                        const Matrix& Xnod,
                        const Matrix& Jac,
                        Matrix& dNdX,
                        Matrix3D* d2NdX2 = nullptr,
                        Matrix3D* Hess = nullptr);

  //! \brief Evaluate the jacobian (and hessian) of the geometry mapping.
  //! \param vol Volume to evaluate
  //! \param x Point to evaluate in
  //! \param elm Element info for point
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Derivatives of basis functions in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalMapping(const LR::LRSplineVolume& vol,
                          const ItgPoint& x,
                          const LR::Element*& elm,
                          Matrix& Xnod,
                          Matrix& Jac,
                          Matrix& dNdX,
                          Matrix3D* d2NdX2 = nullptr,
                          Matrix3D* Hess = nullptr);

  //! \brief Evaluate a basis.
  //! \param vol Volume to evaluate
  //! \param x Point to evaluate in
  //! \param elm Element info for point
  //! \param Xnod Geometry coefficients in point
  //! \param Jac Jacobian of geometry mapping in point
  //! \param dNdX Derivatives of basis functions in point
  //! \param d2NdX2 Second derivatives of basis functions in point
  //! \param Hess Hessian of geometry mapping in point
  static bool evalBasis(const LR::LRSplineVolume& vol,
                        const ItgPoint& x,
                        const LR::Element*& elm,
                        const Matrix& Xnod,
                        const Matrix& Jac,
                        Matrix& dNdX,
                        Matrix3D* d2NdX2 = nullptr,
                        Matrix3D* Hess = nullptr);
};

#endif
