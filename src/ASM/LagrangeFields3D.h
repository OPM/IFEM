// $Id$
//==============================================================================
//!
//! \file LagrangeFields3D.h
//!
//! \date Jun 8 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for Lagrange-based finite element vector fields in 3D.
//!
//==============================================================================

#ifndef _LAGRANGE_FIELDS_3D_H
#define _LAGRANGE_FIELDS_3D_H

#include "Fields.h"

class ASMs3DLag;


/*!
  \brief Class for Lagrange-based finite element vector fields in 3D.

  \details This class implements the methods required to evaluate a 3D Lagrange
  vector field at a given point in parametrical or physical coordinates.
*/

class LagrangeFields3D : public Fields
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] name Name of field
  LagrangeFields3D(const ASMs3DLag* patch, const RealArray& v,
                   char basis = 1, const char* name = nullptr);
  //! \brief Empty destructor.
  virtual ~LagrangeFields3D() {}

  // Methods to compute field values
  //================================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  bool valueNode(size_t node, Vector& vals) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  //! \param[out] vals Values in local point in given element
  bool valueFE(const FiniteElement& fe, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const FiniteElement& fe, Matrix& grad) const;

protected:
  Matrix coord; //!< Matrix of nodel coordinates
  int n1; //!< Number of nodes in first parameter direction
  int n2; //!< Number of nodes in second parameter direction
  int n3; //!< Number of nodes in third parameter direction
  int p1; //!< Element order in first parameter direction
  int p2; //!< Element order in second parameter direction
  int p3; //!< Element order in third parameter direction
};

#endif
