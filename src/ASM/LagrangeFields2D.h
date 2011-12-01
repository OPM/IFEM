// $Id$
//==============================================================================
//!
//! \file LagrangeFields2D.h
//!
//! \date Jun 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element vector fields in 2D.
//!
//==============================================================================

#ifndef _LAGRANGE_FIELDS_2D_H
#define _LAGRANGE_FIELDS_2D_H

#include "Fields.h"

class ASMs2DLag;


/*!
  \brief Class for Lagrange-based finite element vector fields in 2D.

  \details This class implements the methods required to evaluate a 2D Lagrange
  vector field at a given point in parametrical or physical coordinates.
*/

class LagrangeFields2D : public Fields
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] name Name of field
  LagrangeFields2D(const ASMs2DLag* patch, const RealArray& v,
                   const char* name = NULL);
  //! \brief Empty destructor.
  virtual ~LagrangeFields2D() {}

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
  Matrix coord; //!< Matrix of nodal coordinates
  int n1; //!< Number of nodes in first parameter direction
  int n2; //!< Number of nodes in second parameter direction
  int p1; //!< Element order in first parameter direction
  int p2; //!< Element order in second parameter direction
};

#endif
