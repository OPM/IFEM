// $Id$
//==============================================================================
//!
//! \file SplineFields2D.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 2D.
//!
//==============================================================================

#ifndef _SPLINE_FIELDS_2D_H
#define _SPLINE_FIELDS_2D_H

#include "Fields.h"

class ASMs2D;

namespace Go {
  class SplineSurface;
}


/*!
  \brief Class for spline-based finite element vector fields in 2D.

  \details This class implements the methods required to evaluate a 2D
  spline vector field at a given point in parametrical or physical coordinates.
*/

class SplineFields2D : public Fields
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] nf Number of components in field
  //! \param[in] name Name of spline field
  SplineFields2D(const ASMs2D* patch, const RealArray& v,
                 char basis = 1, int nf = 0, const char* name = nullptr);

  //! \brief Construct directly from surface.
  //! \param[in] srf The spline surface to use
  //! \param[in] v Array of control point field values
  //! \param[in] ncmp Number of field components
  //! \param[in] name Name of spline field
  SplineFields2D(const Go::SplineSurface* srf, const RealArray& v, int ncmp,
                 const char* name = nullptr);

  //! \brief Empty destructor.
  virtual ~SplineFields2D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  bool valueNode(size_t node, Vector& vals) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  //! \param[out] vals Values in local point in given element
  bool valueFE(const FiniteElement& fe, Vector& vals) const;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[out] vals Values in given physical coordinate
  bool valueCoor(const Vec4& x, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const FiniteElement& fe, Matrix& grad) const;

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] fe Finite element quantities
  //! \param[out] H Hessian of solution in a given local coordinate
  virtual bool hessianFE(const FiniteElement& fe, Matrix3D& H) const;

protected:
  const Go::SplineSurface* basis; //!< Spline basis description
  const Go::SplineSurface* surf;  //!< Spline geometry description
};

#endif
