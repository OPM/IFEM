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
class HDF5Reader;
class ProcessAdm;


/*!
  \brief Base class for spatial functions, defined through patch-wise fields.
*/

class FieldFuncBase
{
public:
  static constexpr int FIXED_LEVEL = 1 << 24; //!< Bit flag in level for using fixed

protected:
  //! \brief Default constructor.
  FieldFuncBase() : pidx(0), npch(0) {}
  //! \brief No copying of this class.
  FieldFuncBase(const FieldFuncBase&) = delete;
  //! \brief The destructor deletes the patches.
  virtual ~FieldFuncBase();

  //! \brief Sets the active patch.
  bool setPatch(size_t pIdx);

protected:
  std::vector<ASMbase*> patch; //!< The patches on which the field is defined

  size_t pidx; //!< Current patch index
  size_t npch; //!< Number of patches in the field
};


/*!
  \brief Base class for spatial functions, defined from a HDF5-file.
*/

class FieldFuncHDF5 : public FieldFuncBase
{
protected:
  //! \brief The constructor opens the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  explicit FieldFuncHDF5(const std::string& fileName);
  //! \brief No copying of this class.
  FieldFuncHDF5(const FieldFuncHDF5&) = delete;
  //! \brief The destructor closes the HDF5-file.
  virtual ~FieldFuncHDF5();

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
  //! \param[in] nf Number of field components
  //! \param[in] basis Basis to use
  virtual void addPatchField(ASMbase* pch,
                             const std::vector<Real>& coefs,
                             int nf, int basis) = 0;
  //! \brief Clears the field container.
  virtual void clearField() = 0;

  bool hasMultipleLevels; //!< True if we have multiple time levels
  mutable int lastLevel; //!< The last time level read from

private:
  HDF5Reader* hdf5; //!< The HDF5-file containing the field data
  ProcessAdm* pAdm; //!< Process administrator for the HDF5-file reader

  mutable double lastTime;  //!< The time of \a lastLevel
};


/*!
  \brief A scalar-valued spatial function, defined through scalar fields.
*/

class FieldFunction : public RealFunc, private FieldFuncHDF5
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
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

protected:
  //! \brief Evaluates the scalar field function.
  virtual Real evaluate(const Vec3& X) const;

  //! \brief Adds a patch-wise field with the given coefficient values.
  //! \param[in] pch The patch to define the field over
  //! \param[in] coefs Field values
  virtual void addPatchField(ASMbase* pch,
                             const std::vector<Real>& coefs,
                             int, int);
  //! \brief Clears the field container.
  virtual void clearField();

private:
  mutable int currentLevel; //!< Current time level to evaluate at

  std::string fName; //!< Name of field
  std::string bName; //!< Name of basis

  std::vector<Field*> field; //!< The scalar field to be evaluated
};


/*!
  \brief A scalar-valued spatial function, defined through scalar fields.
  \details This class reads its control point values from a stream.
*/

class FieldFuncStream : public RealFunc, private FieldFuncBase
{
public:
  //! \brief The constructor creates a field from the provided input stream.
  //! \param[in] patches List of patches the field is defined over
  //! \param[in] istr Input stream to read control point values from
  FieldFuncStream(const std::vector<ASMbase*>& patches, std::istream& istr);
  //! \brief The destructor deletes the scalar fields.
  virtual ~FieldFuncStream();

  //! \brief Sets the active patch.
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

protected:
  //! \brief Evaluates the scalar field function.
  virtual Real evaluate(const Vec3& X) const;

private:
  std::vector<Field*> field; //!< The scalar field to be evaluated
};


/*!
  \brief Base class for multi-valued spatial functions, defined through fields.
*/

class FieldsFuncBase : public FieldFuncHDF5
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
  //! \param[in] nf Number of field components
  //! \param[in] basis Basis to use
  virtual void addPatchField(ASMbase* pch,
                             const std::vector<Real>& coefs,
                             int nf, int basis);
  //! \brief Clears the field container.
  virtual void clearField();

  //! \brief Evaluates the field at the givent point \b X.
  std::vector<Real> getValues(const Vec3& X);

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
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

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
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

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
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

protected:
  //! \brief Evaluates the tensorial field function.
  virtual SymmTensor evaluate(const Vec3& X) const;
};

#endif
