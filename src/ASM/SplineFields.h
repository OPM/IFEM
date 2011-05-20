//==============================================================================
//!
//! \file SplineFields.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for spline based finite element vector field
//!
//==============================================================================

#ifndef _SPLINE_FIELDS_H
#define _SPLINE_FIELDS_H

#include "Vec3.h"
#include "MatVec.h"
#include "FiniteElement.h"

/* 
   \brief Base class for spline-based finite element vector fields.

   \details This class incapsulates the data and methods needed
   to store and evaluate a spline finite element vector field.
   This is an abstract base class and the fields associated with 
   specific spline objects are implemented as subclasses, for
   instance 1D, 2D and 3D spline formulations.
*/


class SplineFields
{
 protected:
  //! \brief The constructor sets the field name
  //! \param[in] nsd Number of space dimensions (1, 2 or 3)
  //! \param[in] name Name of spline field
 SplineFields(unsigned char nsd_, char* name = NULL) 
   : nsd(nsd_), fieldname(name) {}

 public:
  //! \brief Empty destructor
  virtual ~SplineFields() {}

  // Returns number of space dimensions
  unsigned char getNoSpaceDim() const { return nsd; }

  // Returns number of fields
  unsigned char getNoFields() const { return nf; }

  // Returns number of space dimensions
  int getNoElm() const { return nelm; }

  // Returns number of control points
  int getNoNodes() const { return nno; }

  // Returns name of spline field
  const char* getFieldName() const { return fieldname; }

  // Sets the name of the spline field
  void setFieldName(char* name) { fieldname = name; }

  // Methods to initialize field
  virtual void fill(Vector& vec) 
  { values = vec; nf = values.size()/nno; }

  
  // Methods to compute field values
  //================================================

  //! \brief Computes the value in a given node/control point
  //! \param[in] node Node number 
  //! \param[out] vals Node values 
  virtual bool valueNode(int node, Vector& vals) const = 0;

  //! \brief Computes the value at a given local coordinate
  //! \param[in] fe Finite element definition
  //! \param[out] vals Values in local point in given element
  virtual bool valueFE(const FiniteElement& fe, Vector& vals) const = 0;
  
  //! \brief Computed the value at a given global coordinate
  //! \param[in] x Global/physical coordinate for point
  //! \param[in] vals Values in given physical coordinate
  virtual bool valueCoor(const Vec3& x, Vector& vals) const = 0;

  //! \brief Computes the gradient for a given local coordinate
  //! \param[in] fe Finite element
  //! \param[out] grad Gradient of solution in a given local coordinate
  virtual bool gradFE(const FiniteElement& fe, Matrix& grad) const = 0;

  //! \brief Computes the gradient for a given global/physical coordinate
  //! \param[in] x Global coordinate
  //! \param[out] grad Gradient of solution in a given global coordinate
  virtual bool gradCoor(const Vec3& x, Matrix& grad) const = 0;

 protected:
  // Dimension of field
  unsigned char nsd;  //!< Number of space dimensions
  unsigned char nf;   //!< Number of fields
  int nelm;           //!< Number of elements/knot-spans
  int nno;            //!< Number of nodes/control points

  // Fieldname
  char* fieldname;    //!< Name of spline element field

  // Field values
  Vector values;      //!< Field values
};
   


#endif
