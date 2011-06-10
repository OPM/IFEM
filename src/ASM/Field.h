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

class FiniteElement;
class Vec3;


/*!
  \brief Base class for scalar fields.

  \details This class incapsulates the data and methods needed
  to store and evaluate a scalar field. This is an abstract base 
  class and the fields associated with a specific field type is 
  implemented as subclasses, for instance 1D, 2D and 3D spline 
  formulations or Lagrange formulations.
*/

class Field
{
protected:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] n Number of space dimensions (1, 2 or 3)
  //! \param[in] name Name of field
  Field(unsigned char n, char* name = NULL) : nsd(n), fieldname(name) {}

public:
  //! \brief Empty destructor
  virtual ~Field() {}

  //! \brief Returns number of space dimensions.
  unsigned char getNoSpaceDim() const { return nsd; }

  //! \brief Returns number of elements.
  int getNoElm() const { return nelm; }

  //! \brief Returns number of control points.
  int getNoNodes() const { return nno; }

  //! \brief Returns name of field.
  const char* getFieldName() const { return fieldname; }

  //! \brief Sets the name of the field.
  void setFieldName(const char* name) { fieldname = const_cast<char*>(name); }

  //! \brief Initializes the field values.
  void fill(const Vector& vec) { values = vec; }


  // Methods to compute field values
  //================================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  virtual double valueNode(int node) const = 0;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  virtual double valueFE(const FiniteElement& fe) const = 0;

  //! \brief Computed the value at a given global coordinate.
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
  int nelm;          //!< Number of elements/knot-spans
  int nno;           //!< Number of nodes/control points
  char* fieldname;   //!< Name of field
  Vector values;     //!< Field values
};

#endif
