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

class ASMbase;
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
  //! \brief The constructor sets the field name.
  //! \param[in] name Name of field
  explicit Fields(const char* name = 0) : nf(0), nelm(0), nno(0)
  { if (name) fname = name; }

public:
  //! \brief Empty destructor.
  virtual ~Fields() {}

  //! \brief Returns the number of field components.
  unsigned char getNoFields() const { return nf; }

  //! \brief Returns the number of elements.
  size_t getNoElm() const { return nelm; }

  //! \brief Returns the number of nodal/control points.
  size_t getNoNodes() const { return nno; }

  //! \brief Returns the name of field.
  const char* getFieldName() const { return fname.c_str(); }

  //! \brief Creates a dynamically allocated field object.
  //! \param[in] pch The spline patch on which the field is to be defined on
  //! \param[in] v Array of nodal/control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] name Name of field
  static Fields* create(const ASMbase* pch, const RealArray& v,
                        char basis = 1, const char* name = nullptr);

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  virtual bool valueNode(size_t node, Vector& vals) const = 0;

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

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] fe Finite element
  //! \param[out] d2UdX2 Hessian of solution in a given local coordinate
  virtual bool hessianFE(const FiniteElement& fe, Matrix3D& d2UdX2) const
  { return false; }

  //! \brief Computes the gradient for a given global/physical coordinate.
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  virtual bool gradCoor(const Vec3& x, Matrix& grad) const = 0;

protected:
  unsigned char nf;  //!< Number of field components
  size_t nelm;       //!< Number of elements/knot-spans
  size_t nno;        //!< Number of nodes/control points
  std::string fname; //!< Name of the field
  Vector values;     //!< Field values
};

#endif
