// $Id$
//==============================================================================
//!
//! \file SplineFields1D.h
//!
//! \date Nov 23 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 1D.
//!
//==============================================================================

#ifndef _SPLINE_FIELDS_1D_H
#define _SPLINE_FIELDS_1D_H

#include "Fields.h"

namespace Go {
  class SplineCurve;
}


/*!
  \brief Class for spline-based finite element vector fields in 1D.

  \details This class implements the methods required to evaluate a 1D
  spline vector field at a given point in parametrical or physical coordinates.
*/

class SplineFields1D : public Fields
{
public:
  //! \brief Construct directly from surface.
  //! \param[in] srf The spline surface to use
  //! \param[in] v Array of control point field values
  //! \param[in] ncmp Number of field components
  //! \param[in] name Name of spline field
  SplineFields1D(const Go::SplineCurve* srf, const RealArray& v, int ncmp,
                 const char* name = nullptr);

  //! \brief Empty destructor.
  virtual ~SplineFields1D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  virtual bool valueNode(size_t node, Vector& vals) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  //! \param[out] vals Values in local point in given element
  virtual bool valueFE(const FiniteElement& fe, Vector& vals) const;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[out] vals Values in given physical coordinate
  virtual bool valueCoor(const Vec4& x, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const FiniteElement& fe, Matrix& grad) const;

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] fe Finite element quantities
  //! \param[out] H Hessian of solution in a given local coordinate
  virtual bool hessianFE(const FiniteElement& fe, Matrix3D& H) const;

protected:
  const Go::SplineCurve* curv;  //!< Spline geometry description
};

#endif
