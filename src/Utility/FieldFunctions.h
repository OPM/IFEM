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
class HDF5Writer;
class ProcessAdm;


/*!
  \brief Base class for spatial functions, defined through patch-wise fields.
*/

class FieldFuncBase
{
protected:
  typedef std::vector<Real>      RealArray;  //!< Convenience type
  typedef std::vector<RealArray> RealArrays; //!< Convenience type

  //! \brief The constructor opens the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  explicit FieldFuncBase(const std::string& fileName);
  //! \brief The destructor deletes the patches and close the HDF5-file.
  virtual ~FieldFuncBase();
  //! \brief No copying of this class.
  FieldFuncBase(const FieldFuncBase&) = delete;

  //! \brief Loads field values for the specified time level.
  //! \param[in] fieldNames Name of the field components in the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] level Time level to read for
  //! \param[in] isScalar If \e true, assume this is a scalar field
  bool load(const std::vector<std::string>& fieldNames,
            const std::string& basisName, int level,
            bool isScalar = false);

  //! \brief Finds the level whose time is closest to the specified time.
  int findClosestLevel(double time) const;

  //! \brief Adds a patch-wise field with the given coefficient values.
  //! \param[in] pch The patch to define the field over
  //! \param[in] coefs Field values
  virtual void addPatchField(ASMbase* pch, const RealArrays& coefs) = 0;
  //! \brief Clears the field container.
  virtual void clearField() = 0;

private:
  HDF5Writer* hdf5; //!< The HDF5-file containing the field data
  ProcessAdm* pAdm; //!< Process administrator for the HDF5-file reader

  mutable int    lastLevel; //!< The last time level read from
  mutable double lastTime;  //!< The time of \a lastLevel

  std::vector<ASMbase*> patch; //!< The patches on which the field is defined

protected:
  size_t pidx; //!< Current patch index
};


/*!
  \brief A scalar-valued spatial function, defined through scalar fields.
*/

class FieldFunction : public RealFunc, private FieldFuncBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] level Time level to read for
  FieldFunction(const std::string& fileName,
                const std::string& basisName,
                const std::string& fieldName,
                int level = 0);
  //! \brief The destructor deletes the scalar fields.
  virtual ~FieldFunction() { this->clearField(); }

  //! \brief Sets the active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the scalar field function.
  virtual Real evaluate(const Vec3& X) const;

  //! \brief Adds a patch-wise field with the given coefficient values.
  //! \param[in] pch The patch to define the field over
  //! \param[in] coefs Field values
  virtual void addPatchField(ASMbase* pch, const RealArrays& coefs);
  //! \brief Clears the field container.
  virtual void clearField();

private:
  mutable int currentLevel; //!< Current time level to evaluate at

  std::string fName; //!< Name of field
  std::string bName; //!< Name of basis

  std::vector<Field*> field; //!< The scalar field to be evaluated
};


/*!
  \brief Base class for multi-valued spatial functions, defined through fields.
*/

class FieldsFuncBase : public FieldFuncBase
{
protected:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] level Time level to read for
  FieldsFuncBase(const std::string& fileName,
                 const std::string& basisName,
                 const std::string& fieldName,
                 int level);
  //! \brief The destructor deletes the vector fields.
  virtual ~FieldsFuncBase() { this->clearField(); }

  //! \brief Adds a patch-wise field with the given coefficient values.
  //! \param[in] pch The patch to define the field over
  //! \param[in] coefs Field values
  virtual void addPatchField(ASMbase* pch, const RealArrays& coefs);
  //! \brief Clears the field container.
  virtual void clearField();

  //! \brief Evaluates the field at the givent point \b X.
  RealArray getValues(const Vec3& X);

private:
  mutable int currentLevel; //!< Current time level to evaluate at

protected:
  std::vector<std::string> fName; //!< Name of field components
  std::string              bName; //!< Name of basis

  std::vector<Fields*> field; //!< The vector field to be evaluated
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
  //! \param[in] level Time level to read for
  VecFieldFunction(const std::string& fileName,
                   const std::string& basisName,
                   const std::string& fieldName,
                   int level = 0);
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
  //! \param[in] level Time level to read for
  TensorFieldFunction(const std::string& fileName,
                      const std::string& basisName,
                      const std::string& fieldName,
                      int level = 0);
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
  //! \param[in] level Time level to read for
  STensorFieldFunction(const std::string& fileName,
                       const std::string& basisName,
                       const std::string& fieldName,
                       int level = 0);
  //! \brief Empty destructor.
  virtual ~STensorFieldFunction() {}

  //! \brief Sets the active patch.
  virtual void initPatch(size_t pIdx) { if (pIdx < field.size()) pidx = pIdx; }

protected:
  //! \brief Evaluates the tensorial field function.
  virtual SymmTensor evaluate(const Vec3& X) const;
};

#endif
