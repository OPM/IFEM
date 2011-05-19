//==============================================================================
//!
//! \file SplineField.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for spline based finite element field
//!
//==============================================================================

#ifndef _SPLINE_FIELD_H
#define _SPLINE_FIELD_H

#include "Vec3.h"
#include "MatVec.h"
#include "FiniteElement.h"

/* 
   \brief Base class for spline-based finite element fields.

   \details This class incapsulates the data and methods needed
   to store and evaluate a spline finite element scalar field.
   This is an abstract base class and the fields associated with 
   specific spline objects are implemented as subclasses, for
   instance 1D, 2D and 3D spline formulations.
*/


class SplineField
{
 protected:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] nsd Number of space dimensions (1, 2 or 3)
  //! \param[in] name Name of spline field
 SplineField(unsigned char nsd_, char* name = NULL) 
   : fieldname(name) {}

 public:
  //! \brief Empty destructor
  virtual ~SplineField() {}

  // Returns number of space dimensions
  unsigned char getNoSpaceDim() const { return nsd; }

  // Returns number of space dimensions
  int getNoElm() const { return nelm; }

  // Returns number of control points
  int getNoNodes() const { return nno; }

  // Returns name of spline field
  const char* getFieldName() const { return fieldname; }

  // Sets the name of the spline field
  void setFieldName(char* name) { fieldname = name; }

  // Methods to initialize field
  virtual void fill(Vector& vec) { values = vec; }

  
  // Methods to compute field values
  //================================================

  //! \brief Computes the value in a given node/control point
  //! \param[in] node Node number 
  virtual double valueNode(int node) const = 0;

  //! \brief Computes the value at a given local coordinate
  //! \param[in] fe Finite element definition
  virtual double valueFE(const FiniteElement& fe) const = 0;
  
  //! \brief Computed the value at a given global coordinate
  //! \param[in] x Global/physical coordinate for point
  virtual double valueCoor(const Vec3 x) const = 0;

  //! \brief Computes the gradient for a given local coordinate
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const FiniteElement& fe, Vector& grad) const = 0;

  //! \brief Computes the gradient for a given global/physical coordinate
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  virtual bool gradCoor(const Vec3 x, Vector& grad) const = 0;

 protected:
  // Dimension of field
  unsigned char nsd;  //!< Number of space dimensions
  int nelm;           //!< Number of elements/knot-spans
  int nno;            //!< Number of nodes/control points

  // Fieldname
  char* fieldname;    //!< Name of spline element field

  // Field values
  Vector values;      //!< Field values
};
   


#endif
