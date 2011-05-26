// $Id$
//==============================================================================
//!
//! \file SplineFields.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for spline-based finite element vector fields.
//!
//==============================================================================

#ifndef _SPLINE_FIELDS_H
#define _SPLINE_FIELDS_H

#include "MatVec.h"

class FiniteElement;
class Vec3;


/*!
  \brief Base class for spline-based finite element vector fields.

  \details This class incapsulates the data and methods needed
  to store and evaluate a spline finite element vector field.
  This is an abstract base class and the fields associated with
  specific spline objects are implemented as subclasses,
  for instance 1D, 2D and 3D spline formulations.
*/

class SplineFields
{
protected:
  //! \brief The constructor sets the field name.
  //! \param[in] n Number of space dimensions (1, 2 or 3)
  //! \param[in] name Name of spline field
  SplineFields(unsigned char n, char* name = NULL) : nsd(n), fieldname(name) {}

public:
  //! \brief Empty destructor.
  virtual ~SplineFields() {}

  //! \brief Returns number of space dimensions.
  unsigned char getNoSpaceDim() const { return nsd; }

  //! \brief Returns number of fields.
  unsigned char getNoFields() const { return nf; }

  //! \brief Returns number of elements.
  int getNoElm() const { return nelm; }

  //! \brief Returns number of control points.
  int getNoNodes() const { return nno; }

  //! \brief Returns name of spline field.
  const char* getFieldName() const { return fieldname; }

  //! \brief Sets the name of the spline field.
  void setFieldName(char* name) { fieldname = name; }

  //! \brief Initializes the field values.
  void fill(const Vector& vec) { values = vec; nf = values.size()/nno; }


  // Methods to compute field values
  //================================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  virtual bool valueNode(int node, Vector& vals) const = 0;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] fe Finite element definition
  //! \param[out] vals Values in local point in given element
  virtual bool valueFE(const FiniteElement& fe, Vector& vals) const = 0;

  //! \brief Computed the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[in] vals Values in given physical coordinate
  virtual bool valueCoor(const Vec3& x, Vector& vals) const = 0;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const FiniteElement& fe, Matrix& grad) const = 0;

  //! \brief Computes the gradient for a given global/physical coordinate.
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  virtual bool gradCoor(const Vec3& x, Matrix& grad) const = 0;

protected:
  unsigned char nsd; //!< Number of space dimensions
  unsigned char nf;  //!< Number of fields
  int nelm;          //!< Number of elements/knot-spans
  int nno;           //!< Number of nodes/control points
  char* fieldname;   //!< Name of spline element field
  Vector values;     //!< Field values
};

#endif
