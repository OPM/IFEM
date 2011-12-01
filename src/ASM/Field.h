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
class FiniteElement;
class Vec3;


/*!
  \brief Base class for scalar fields.

  \details This class incapsulates the data and methods needed to store and
  evaluate a scalar field. This is an abstract base class, and the fields
  associated with a specific field type are implemented as subclasses,
  for instance 1D, 2D and 3D spline formulations or Lagrange formulations.
*/

class Field
{
protected:
  //! \brief The constructor sets the number of space dimensions and field name.
  //! \param[in] n Number of space dimensions (1, 2 or 3)
  //! \param[in] name Name of field
  Field(unsigned char n, const char* name = 0) : nsd(n), nelm(0), nno(0)
  { if (name) fname = name; }

public:
  //! \brief Empty destructor.
  virtual ~Field() {}

  //! \brief Returns number of space dimensions.
  unsigned char getNoSpaceDim() const { return nsd; }

  //! \brief Returns number of elements.
  size_t getNoElm() const { return nelm; }

  //! \brief Returns number of nodal/control points.
  size_t getNoNodes() const { return nno; }

  //! \brief Returns the name of field.
  const char* getFieldName() const { return fname.c_str(); }

  //! \brief Creates a dynamically allocated field object.
  //! \param[in] pch The spline patch on which the field is to be defined on
  //! \param[in] v Array of nodal/control point field values
  //! \param[in] name Name of field
  static Field* create(const ASMbase* pch, const RealArray& v,
		       const char* name = NULL);

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  virtual double valueNode(size_t node) const = 0;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  virtual double valueFE(const FiniteElement& fe) const = 0;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  virtual double valueCoor(const Vec3& x) const = 0;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const FiniteElement& fe, Vector& grad) const = 0;

  //! \brief Computes the gradient for a given global/physical coordinate.
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  virtual bool gradCoor(const Vec3& x, Vector& grad) const = 0;

protected:
  unsigned char nsd; //!< Number of space dimensions
  size_t nelm;       //!< Number of elements/knot-spans
  size_t nno;        //!< Number of nodes/control points
  std::string fname; //!< Name of the field
  Vector values;     //!< Field values
};

#endif
