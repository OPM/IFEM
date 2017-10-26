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
#include "TensorFunction.h"
#include <string>

class Field;
class Fields;
class ASMbase;


/*!
  \brief A scalar-valued spatial function, defined trough a scalar field.
*/

class FieldFunction : public RealFunc
{
  std::vector<Field*> field; //!< Pointer to the scalar field to be evaluated
  std::vector<ASMbase*> patch; //!< Pointer to the patch on which the field is defined

public:
  //! \brief Default constructor.
  explicit FieldFunction(Field* f = nullptr) : field{f}, pidx(0) {}
  //! \brief Constructor creating a field from a provided HDF5 file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] npatches Number of patches to read for
  //! \param[in] level Level to read for
  FieldFunction(const std::string& fileName,
                const std::string& basisName,
                const std::string& fieldName,
                size_t npatches = 1, int level = 0);

  //! \brief The destructor deletes the scalar fields.
  virtual ~FieldFunction();

  //! \brief Set currently active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the scalar field function.
  virtual Real evaluate(const Vec3& X) const;

  size_t pidx; //!< Current patch index
};

/*!
  \brief Base class for multi-value spatial function, defined trough a vector field.
*/

class VecFieldFuncBase
{
protected:
  //! \brief Default constructor.
  VecFieldFuncBase(Fields* f = nullptr) : field{f} {}

  //! \brief Constructor creating a field from a provided HDF5 file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] patch Number of patches to read for
  //! \param[in] level Level to read for
  VecFieldFuncBase(const std::string& fileName,
                   const std::string& basisName,
                   const std::string& fieldName,
                   size_t npatches = 1, int level = 0);

  //! \brief The destructor deletes the vector field.
  virtual ~VecFieldFuncBase();

  std::vector<Fields*> field; //!< Pointer to the vector field to be evaluated
  std::vector<ASMbase*> patch; //!< Pointer to the patch on which the field is defined
  size_t pidx; //!< Current patch index
};


/*!
  \brief A vector-valued spatial function, defined trough a vector field.
*/

class VecFieldFunction : public VecFunc, public VecFieldFuncBase
{
public:
  //! \brief Default constructor.
  explicit VecFieldFunction(Fields* f = nullptr) : VecFieldFuncBase(f) {}

  //! \brief Constructor creating a field from a provided HDF5 file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Level to read for
  VecFieldFunction(const std::string& fileName,
                   const std::string& basisName,
                   const std::string& fieldName,
                   size_t nPatches = 1, int level = 0) :
    VecFieldFuncBase(fileName, basisName, fieldName, nPatches, level) {}

  //! \brief Empty destructor.
  virtual ~VecFieldFunction() {}

  //! \brief Set currently active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the vectorial field function.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief A tensor-valued spatial function, defined trough a vector field.
*/

class TensorFieldFunction : public TensorFunc, public VecFieldFuncBase
{
public:
  //! \brief Default constructor.
  explicit TensorFieldFunction(Fields* f = nullptr) : VecFieldFuncBase(f) {}

  //! \brief Constructor creating a field from a provided HDF5 file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Level to read for
  TensorFieldFunction(const std::string& fileName,
                      const std::string& basisName,
                      const std::string& fieldName,
                      size_t nPatches = 1, int level = 0) :
    VecFieldFuncBase(fileName, basisName, fieldName, nPatches, level) {}

  //! \brief Empty destructor.
  virtual ~TensorFieldFunction() {}

  //! \brief Set currently active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the tensorial field function.
  virtual Tensor evaluate(const Vec3& X) const;
};


/*!
  \brief A symmtensor-valued spatial function, defined trough a vector field.
*/

class STensorFieldFunction : public STensorFunc, public VecFieldFuncBase
{
public:
  //! \brief Default constructor.
  explicit STensorFieldFunction(Fields* f = nullptr) : VecFieldFuncBase(f) {}

  //! \brief Constructor creating a field from a provided HDF5 file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] nPatches Number of patches to read for
  //! \param[in] level Level to read for
  STensorFieldFunction(const std::string& fileName,
                       const std::string& basisName,
                       const std::string& fieldName,
                       size_t nPatches = 1, int level = 0) :
    VecFieldFuncBase(fileName, basisName, fieldName, nPatches, level) {}

  //! \brief Empty destructor.
  virtual ~STensorFieldFunction() {}

  //! \brief Set currently active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the tensorial field function.
  virtual SymmTensor evaluate(const Vec3& X) const;
};

#endif
