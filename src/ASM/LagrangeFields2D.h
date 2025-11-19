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
  //! \param[in] basis Basis to use from patch
  //! \param[in] name Name of field
  LagrangeFields2D(const ASMs2DLag* patch, const RealArray& v,
                   char basis = 1, const char* name = nullptr);
  //! \brief Empty destructor.
  virtual ~LagrangeFields2D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  bool valueNode(size_t node, Vector& vals) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] vals Values in local point in given element
  bool valueFE(const ItgPoint& x, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const ItgPoint& x, Matrix& grad) const;

protected:
  using IntMat = std::vector<std::vector<int>>; //!< Convenience type
  Matrix coord; //!< Matrix of nodal coordinates
  IntMat mnpc; //!< Matrix of element nodes
  int p1; //!< Element order in first parameter direction
  int p2; //!< Element order in second parameter direction
};

#endif
