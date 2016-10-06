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

#include "FieldBase.h"

class ASMs2DLag;


/*!
  \brief Class for Lagrange-based finite element scalar fields in 2D.

  \details This class implements the methods required to evaluate a 2D Lagrange
  scalar field at a given point in parametrical or physical coordinates.
*/

class LagrangeField2D : public FieldBase
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] cmp Component to use
  //! \param[in] name Name of field
  LagrangeField2D(const ASMs2DLag* patch, const RealArray& v,
                  char basis = 1, char cmp = 1, const char* name = nullptr);
  //! \brief Empty destructor.
  virtual ~LagrangeField2D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  double valueNode(size_t node) const;

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
