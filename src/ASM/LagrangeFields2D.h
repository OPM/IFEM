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

/*!
  \brief Class for Lagrange-based finite element vector fields in 2D.

  \details This class implements the functions required to evaluate a 2D
  Lagrange vector field at a given point in parametrical or physical coordinates.
*/


class LagrangeFields2D : public Fields
{
public:
  //! \brief The constructor sets the field name.
  //! \param[in] X  Matrix of nodel coordinates
  //! \param[in] n1 Number of nodes in first parameter direction
  //! \param[in] n2 Number of nodes in second parameter direction
  //! \param[in] p1 Element order in first parameter direction
  //! \param[in] p2 Element order in second parameter direction
  //! \param[in] name Name of spline field
  LagrangeFields2D(Matrix X, int nx, int ny, int px, int py, char* name = NULL);
  //! \brief Empty destructor.
  virtual ~LagrangeFields2D() {}

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
  Matrix coord;  //!< Matrix of nodel coordinates
  int n1, n2;    //!< Number of nodes in each parameter direction
  int p1, p2;    //!< Element order in each parameter direction
};

#endif
