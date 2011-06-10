// $Id$
//==============================================================================
//!
//! \file SplineFields3D.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 3D.
//!
//==============================================================================

#ifndef _SPLINE_FIELDS_3D_H
#define _SPLINE_FIELDS_3D_H

#include "Fields.h"

namespace Go {
  class SplineVolume;
}


/*!
  \brief Class for spline-based finite element vector fields in 3D.

  \details This class implements the functions required to evaluate a 3D
  spline vector field at a given point in parametrical or physical coordinates.
*/


class SplineFields3D : public Fields
{
public:
  //! \brief The constructor sets the field name.
  //! \param[in] geometry Spline volume geometry
  //! \param[in] name Name of spline field
  SplineFields3D(Go::SplineVolume *geometry, char* name = NULL);
  //! \brief Empty destructor.
  virtual ~SplineFields3D();

  // Methods to compute field values
  //================================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  bool valueNode(int node, Vector& vals) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  //! \param[out] vals Values in local point in given element
  bool valueFE(const FiniteElement& fe, Vector& vals) const;

  //! \brief Computed the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[in] vals Values in given physical coordinate
  bool valueCoor(const Vec3& x, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const FiniteElement& fe, Matrix& grad) const;

  //! \brief Computes the gradient for a given global/physical coordinate.
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  bool gradCoor(const Vec3& x, Matrix& grad) const;

protected:
  Go::SplineVolume* vol; //!< Spline volume geometry description
};

#endif
