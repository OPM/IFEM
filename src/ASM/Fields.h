// $Id$
//==============================================================================
//!
//! \file Fields.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for vector fields.
//!
//==============================================================================

#ifndef _FIELDS_H
#define _FIELDS_H

#include "MatVec.h"
#include <string>

class FiniteElement;
class Vec3;


/*!
  \brief Base class for vector fields.

  \details This class incapsulates the data and methods needed to store and
  evaluate a vector field. This is an abstract base class, and the fields
  associated with a specific field type are implemented as subclasses,
  for instance 1D, 2D and 3D spline formulations or Lagrange formulations.
*/

class Fields
{
protected:
  //! \brief The constructor sets the number of space dimensions and field name.
  //! \param[in] n Number of space dimensions (1, 2 or 3)
  //! \param[in] name Name of field
  Fields(unsigned char n, const char* name = 0) : nsd(n), nf(0), nelm(0), nno(0)
  { if (name) fname = name; }

public:
  //! \brief Empty destructor.
  virtual ~Fields() {}

  //! \brief Returns number of space dimensions.
  unsigned char getNoSpaceDim() const { return nsd; }

  //! \brief Returns number of field components.
  unsigned char getNoFields() const { return nf; }

  //! \brief Returns number of elements.
  size_t getNoElm() const { return nelm; }

  //! \brief Returns number of control points.
  size_t getNoNodes() const { return nno; }

  //! \brief Returns name of field.
  const char* getFieldName() const { return fname.c_str(); }

  //! \brief Sets the name of the field.
  void setFieldName(const char* name) { fname = name; }

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

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[out] vals Values in given physical coordinate
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
  unsigned char nf;  //!< Number of field components
  size_t nelm;       //!< Number of elements/knot-spans
  size_t nno;        //!< Number of nodes/control points
  std::string fname; //!< Name of the field
  Vector values;     //!< Field values
};

#endif
