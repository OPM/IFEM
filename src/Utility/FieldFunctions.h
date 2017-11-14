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

#include "TensorFunction.h"
#include <string>

class Field;
class Fields;
class ASMbase;


/*!
  \brief A scalar-valued spatial function, defined through scalar fields.
*/

class FieldFunction : public RealFunc
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Time level to read for
  FieldFunction(const std::string& fileName,
                const std::string& basisName,
                const std::string& fieldName,
                size_t nPatches = 1, int level = 0);
  //! \brief The destructor deletes the scalar fields.
  virtual ~FieldFunction();

  //! \brief Sets the active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the scalar field function.
  virtual Real evaluate(const Vec3& X) const;

private:
  std::vector<Field*>   field; //!< The scalar field to be evaluated
  std::vector<ASMbase*> patch; //!< The patches on which the field is defined

  size_t pidx; //!< Current patch index
};


/*!
  \brief Base class for multi-valued spatial functions, defined through fields.
*/

class FieldsFuncBase
{
protected:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Time level to read for
  FieldsFuncBase(const std::string& fileName,
                 const std::string& basisName,
                 const std::string& fieldName,
                 size_t nPatches, int level);
  //! \brief The destructor deletes the vector fields.
  virtual ~FieldsFuncBase();

  std::vector<Fields*>  field; //!< The vector field to be evaluated
  std::vector<ASMbase*> patch; //!< The patches on which the field is defined

  size_t pidx; //!< Current patch index
};


/*!
  \brief A vector-valued spatial function, defined through a vector field.
*/

class VecFieldFunction : public VecFunc, private FieldsFuncBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Time level to read for
  VecFieldFunction(const std::string& fileName,
                   const std::string& basisName,
                   const std::string& fieldName,
                   size_t nPatches = 1, int level = 0);
  //! \brief Empty destructor.
  virtual ~VecFieldFunction() {}

  //! \brief Sets the active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the vectorial field function.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief A tensor-valued spatial function, defined through a vector field.
*/

class TensorFieldFunction : public TensorFunc, private FieldsFuncBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Time level to read for
  TensorFieldFunction(const std::string& fileName,
                      const std::string& basisName,
                      const std::string& fieldName,
                      size_t nPatches = 1, int level = 0);
  //! \brief Empty destructor.
  virtual ~TensorFieldFunction() {}

  //! \brief Sets the active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the tensorial field function.
  virtual Tensor evaluate(const Vec3& X) const;
};


/*!
  \brief A symmtensor-valued spatial function, defined through a vector field.
*/

class STensorFieldFunction : public STensorFunc, private FieldsFuncBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Time level to read for
  STensorFieldFunction(const std::string& fileName,
                       const std::string& basisName,
                       const std::string& fieldName,
                       size_t nPatches = 1, int level = 0);
  //! \brief Empty destructor.
  virtual ~STensorFieldFunction() {}

  //! \brief Sets the active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the tensorial field function.
  virtual SymmTensor evaluate(const Vec3& X) const;
};

#endif
