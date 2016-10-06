// $Id$
//==============================================================================
//!
//! \file SplineField3D.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element scalar fields in 3D.
//!
//==============================================================================

#ifndef _SPLINE_FIELD_3D_H
#define _SPLINE_FIELD_3D_H

#include "FieldBase.h"

class ASMs3D;

namespace Go {
  class SplineVolume;
}


/*!
  \brief Class for spline-based finite element scalar fields in 3D.

  \details This class implements the functions required to evaluate a 3D
  spline scalar field at a given point in parametrical or physical coordinates.
*/

class SplineField3D : public FieldBase
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] cmp Component to use from source field. Pass 0 to use vector as-is.
  //! \param[in] name Name of spline field
  SplineField3D(const ASMs3D* patch, const RealArray& v,
                char basis = 1, char cmp = 1, const char* name = nullptr);
  //! \brief Empty destructor.
  virtual ~SplineField3D() {}

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
  virtual double valueCoor(const Vec3& x) const;

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

  //! \brief Computes the gradient for a given global/physical coordinate.
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  virtual bool gradCoor(const Vec3& x, Vector& grad) const;

protected:
  const Go::SplineVolume* basis; //!< Spline basis description
  const Go::SplineVolume* vol;   //!< Spline geometry description
};

#endif
