// $Id$
//==============================================================================
//!
//! \file LagrangeField2D.h
//!
//! \date June 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element scalar fields in 2D.
//!
//==============================================================================

#ifndef _LAGRANGE_FIELD_2D_H
#define _LAGRANGE_FIELD_2D_H

#include "Field.h"


/*!
  \brief Class for Lagrange-based finite element scalar fields in 2D.

  \details This class implements the methods required to evaluate a 2D Lagrange
  scalar field at a given point in parametrical or physical coordinates.
*/

class LagrangeField2D : public Field
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] X  Matrix of nodal coordinates
  //! \param[in] nx Number of nodes in first parameter direction
  //! \param[in] ny Number of nodes in second parameter direction
  //! \param[in] px Element order in first parameter direction
  //! \param[in] py Element order in second parameter direction
  //! \param[in] name Name of field
  LagrangeField2D(const Matrix& X,
		  int nx, int ny, int px, int py, char* name = NULL);
  //! \brief Empty destructor.
  virtual ~LagrangeField2D() {}

  // Methods to compute field values
  //================================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  double valueNode(int node) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  double valueFE(const FiniteElement& fe) const;

  //! \brief Computes the value at a given global coordinate.
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
  Matrix coord; //!< Matrix of nodal coordinates
  int n1; //!< Number of nodes in first parameter direction
  int n2; //!< Number of nodes in second parameter direction
  int p1; //!< Element order in first parameter direction
  int p2; //!< Element order in second parameter direction
};

#endif
