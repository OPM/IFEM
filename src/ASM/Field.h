// $Id$
//==============================================================================
//!
//! \file Field.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for scalar fields.
//!
//==============================================================================

#ifndef _FIELD_H
#define _FIELD_H

#include "MatVec.h"
#include <string>

class ASMbase;
class ItgPoint;
class Vec4;


/*!
  \brief Interface class for scalar fields.
*/

class Field
{
protected:
  //! \brief The constructor sets the field name.
  //! \param[in] name Optional name of field
  explicit Field(const char* name = nullptr) { if (name) fname = name; }

public:
  //! \brief Empty destructor.
  virtual ~Field() {}

  //! \brief Creates a dynamically allocated field object.
  //! \param[in] pch The spline patch on which the field is to be defined on
  //! \param[in] v Array of nodal/control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] cmp Component to use for field
  //! \param[in] name Name of field
  static Field* create(const ASMbase* pch, const RealArray& v,
                       char basis = 1, char cmp = 1,
                       const char* name = nullptr);

  //! \brief Returns the name of field.
  const char* getFieldName() const { return fname.c_str(); }

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node 1-based node/control point index
  virtual double valueNode(size_t node) const = 0;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  virtual double valueFE(const ItgPoint& x) const = 0;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  virtual double valueCoor(const Vec4& x) const { return 0.0; }

  //! \brief Computes the value at a grid of visualization points.
  //! \param[out] val Field values at the visualization points
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool valueGrid(RealArray& val, const int* npe) const { return false; }

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const ItgPoint& x, Vector& grad) const = 0;

  //! \brief Computes the gradient for a given global/physical coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[out] grad Gradient of solution in a given global coordinate
  virtual bool gradCoor(const Vec4& x, Vector& grad) const { return false; }

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] H Hessian of solution in a given local coordinate
  virtual bool hessianFE(const ItgPoint& x, Matrix& H) const { return false; }

protected:
  std::string fname; //!< Name of the field
};

#endif
