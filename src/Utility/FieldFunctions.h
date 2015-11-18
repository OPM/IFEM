// $Id$
//==============================================================================
//!
//! \file FieldFunctions.h
//!
//! \date Sep 20 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Field function implementations.
//!
//==============================================================================

#ifndef _FIELD_FUNCTIONS_H
#define _FIELD_FUNCTIONS_H

#include "Function.h"
#include <string>

class Field;
class ASMbase;


/*!
  \brief A scalar-valued spatial function, defined trough a scalar field.
*/

class FieldFunction : public RealFunc
{
  Field* field; //!< Pointer to the scalar field to be evaluated.
  ASMbase* pch; //!< Pointer to the patch on which the field is defined over

public:
  //! \brief Default constructor.
  FieldFunction(Field* f = nullptr) : field(f), pch(nullptr) {}
  //! \brief Constructor creating a field from a provided HDF5 file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  FieldFunction(const std::string& fileName,
                const std::string& basisName,
                const std::string& fieldName);
  //! \brief The destructor deletes the scalar field.
  virtual ~FieldFunction();

protected:
  //! \brief Evaluates the scalar field function.
  virtual Real evaluate(const Vec3& X) const;
};

#endif
