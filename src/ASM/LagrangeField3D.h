//==============================================================================
//!
//! \file LagrangeField3D.h
//!
//! \date Jun 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element scalar fields in 3D.
//!
//==============================================================================

#ifndef _LAGRANGE_FIELD_3D_H
#define _LAGRANGE_FIELD_3D_H

#include "Field.h"

/*!
  \brief Class for Lagrange-based finite element scalar fields in 3D.

  \details This class implements the functions required to evaluate a 3D
  Lagrange scalar field at a given point in parametrical or physical coordinates.
*/


class LagrangeField3D : public Field
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] geometry Spline volume geometry
  //! \param[in] name Name of spline field
  LagrangeField3D(Matrix X, int nx, int ny, int nz,
		  int px, int py, int pz, char* name = NULL);
  //! \brief Empty destructor.
  virtual ~LagrangeField3D() {}

  // Methods to compute field values
  //================================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  double valueNode(int node) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  double valueFE(const FiniteElement& fe) const;

  //! \brief Computed the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  double valueCoor(const Vec3& x) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const FiniteElement& fe, Vector& grad) const;

  //! \brief Computes the gradient for a given global/physical coordinate.
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  bool gradCoor(const Vec3& x, Vector& grad) const;

 protected:
  Matrix coord;   //!< Matrix of nodel coordinates
  int n1, n2, n3; //!< Number of nodes in each parameter direction
  int p1, p2, p3; //!< Element order in each parameter direction
};

#endif
