// $Id$
//==============================================================================
//!
//! \file SplineField1D.h
//!
//! \date Nov 3 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for spline-based finite element scalar fields in 1D.
//!
//==============================================================================

#ifndef _SPLINE_FIELD_1D_H
#define _SPLINE_FIELD_1D_H

#include "FieldBase.h"

class ASMs1D;

namespace Go {
  class SplineCurve;
}


/*!
  \brief Class for spline-based finite element scalar fields in 2D.

  \details This class implements the methods required to evaluate a 2D
  spline scalar field at a given point in parametrical or physical coordinates.
*/

class SplineField1D : public FieldBase
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] name Name of spline field
  SplineField1D(const ASMs1D* patch, const RealArray& v, const char* name=nullptr);
  //! \brief Empty destructor.
  virtual ~SplineField1D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  virtual double valueNode(size_t node) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  virtual double valueFE(const FiniteElement& fe) const;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  virtual double valueCoor(const Vec4& x) const;

  //! \brief Computes the value at a grid of visualization points.
  //! \param[out] val Field values at the visualization points
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool valueGrid(RealArray& val, const int* npe) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const FiniteElement& fe, Vector& grad) const;

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] fe Finite element quantities
  //! \param[out] H Hessian of solution in a given local coordinate
  virtual bool hessianFE(const FiniteElement& fe, Matrix& H) const;

protected:
  const Go::SplineCurve* curv;  //!< Spline geometry description
};

#endif
