// $Id$
//==============================================================================
//!
//! \file FieldBase.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for scalar fields.
//!
//==============================================================================

#ifndef _FIELD_BASE_H
#define _FIELD_BASE_H

#include "Field.h"


/*!
  \brief Base class for scalar fields.

  \details This class incapsulates the data and methods needed to store and
  evaluate a scalar field. This is an abstract base class, and the fields
  associated with a specific field type are implemented as subclasses,
  for instance 1D, 2D and 3D spline formulations or Lagrange formulations.
*/

class FieldBase : public Field
{
protected:
  //! \brief The constructor sets the field name.
  //! \param[in] name Optional name of field
  FieldBase(const char* name = nullptr) : Field(name) { nelm = nno = 0; }

public:
  //! \brief Empty destructor.
  virtual ~FieldBase() {}

  //! \brief Returns the number of elements.
  size_t getNoElm() const { return nelm; }
  //! \brief Returns the number of nodal/control points.
  size_t getNoNodes() const { return nno; }

protected:
  size_t nelm;   //!< Number of elements/knot-spans
  size_t nno;    //!< Number of nodes/control points
  Vector values; //!< Nodal field values
};

#endif
