// $Id$
//==============================================================================
//!
//! \file LRSplineField2D.h
//!
//! \date Mar 9 2017
//!
//! \author Arne Moreten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element scalar fields in 2D.
//!
//==============================================================================

#ifndef _LRSPLINE_FIELD_2D_H
#define _LRSPLINE_FIELD_2D_H

#include "FieldBase.h"

class ASMu2D;

namespace LR {
  class Element;
  class LRSplineSurface;
}


/*!
  \brief Class for LR spline-based finite element scalar fields in 2D.

  \details This class implements the methods required to evaluate a 2D LR spline
  scalar field at a given point in parametrical or physical coordinates.
*/

class LRSplineField2D : public FieldBase
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] cmp Component to use from source field.
  //! Pass 0 to use vector as-is.
  //! \param[in] name Name of spline field
  LRSplineField2D(const ASMu2D* patch, const RealArray& v,
                  char basis = 1, char cmp = 1, const char* name = nullptr);
  //! \brief Construct directly from surface.
  //! \param[in] srf The spline surface to use
  //! \param[in] v Array of control point field values
  //! \param[in] name Name of spline field
  LRSplineField2D(const LR::LRSplineSurface* srf, const RealArray& v,
                  const char* name = nullptr);
  //! \brief Empty destructor.
  virtual ~LRSplineField2D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  virtual double valueNode(size_t node) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  virtual double valueFE(const ItgPoint& x) const;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  virtual double valueCoor(const Vec4& x) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const ItgPoint& x, Vector& grad) const;

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] H Hessian of solution in a given local coordinate
  virtual bool hessianFE(const ItgPoint& x, Matrix& H) const;

protected:
  const LR::LRSplineSurface* basis; //!< Spline basis description
  const LR::LRSplineSurface* surf;  //!< Spline geometry description
};

#endif
