//==============================================================================
//!
//! \file SplineField2D.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for spline based finite element field in 2D
//!
//==============================================================================

#ifndef _SPLINE_FIELD_2D_H
#define _SPLINE_FIELD_2D_H

namespace Go {
  class SplineSurface;
  class Point;
}

#include "SplineField.h"
#include "GoTools/geometry/SplineSurface.h"

/* 
   \brief Base class for spline-based finite element fields in 2D.

   \details This class implements the functions required to evaluate
   a 2D spline field at a given point in parametrical or physical 
   coordinates.
*/


class SplineField2D : public SplineField
{
 public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] geometry Spline volume geometry
  //! \param[in] name Name of spline field
  SplineField2D(Go::SplineSurface *geometry, char* name = NULL); 

  //! \brief Empty destructor
  virtual ~SplineField2D();

  // Methods to compute field values
  //================================================

  //! \brief Computes the value in a given node/control point
  //! \param[in] node Node number 
  double valueNode(int node) const;

  //! \brief Computes the value at a given local coordinate
  //! \param[in] fe Finite element definition
  double valueFE(const FiniteElement& fe) const;

  //! \brief Computed the value at a given global coordinate
  //! \param[in] x Global/physical coordinate for point
  double valueCoor(const Vec3& x) const;

  //! \brief Computes the gradient for a given local coordinate
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const FiniteElement& fe, Vector& grad) const;

  //! \brief Computes the gradient for a given global/physical coordinate
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  bool gradCoor(const Vec3& x, Vector& grad) const;

 protected:
  Go::SplineSurface *surf;    //!< Spline surface geometry description

};
   


#endif
