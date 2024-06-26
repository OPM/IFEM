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
  //! \brief Default constructor.
  FieldFuncHDF5() : hasMultipleLevels(false), lastLevel(-1),
                    hdf5(nullptr), pAdm(nullptr)
  {}
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
  \brief Base class for field functions derived from a scalar field.
*/

class FieldFuncScalarBase : protected FieldFuncHDF5
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] level Time level to read for
  FieldFuncScalarBase(const std::string& fileName,
                      const std::string& basisName,
                      const std::string& fieldName,
                      int level = 0);
  //! \brief The destructor deletes the scalar fields.
  virtual ~FieldFuncScalarBase() { this->clearField(); }

protected:
  //! \brief Default constructor
  FieldFuncScalarBase() : currentLevel(-1) {}

  //! \brief Adds a patch-wise field with the given coefficient values.
  //! \param[in] pch The patch to define the field over
  //! \param[in] coefs Field values
  virtual void addPatchField(ASMbase* pch,
                             const std::vector<Real>& coefs,
                             int, int);
  //! \brief Clears the field container.
  virtual void clearField();

  mutable int currentLevel; //!< Current time level to evaluate at

  std::string fName; //!< Name of field
  std::string bName; //!< Name of basis

  std::vector<Field*> field; //!< The scalar field to be evaluated
};


/*!
  \brief A scalar-valued spatial function, defined through scalar fields.
*/

class FieldFunction : public RealFunc, public FieldFuncScalarBase
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
  //! \brief Wrap existing field.
  //! \param fields Field to evaluate
  //!
  //! \details Takes ownership of given fields.
  FieldFunction(const std::vector<Field*>& fields);
  //! \brief Empty destructor.
  virtual ~FieldFunction() {}

  //! \brief Sets the active patch.
  bool initPatch(size_t pIdx) override { return this->setPatch(pIdx); }

  //! \brief Evaluates first derivatives of the function.
  Vec3 gradient(const Vec3& X) const override;

  //! \brief Evaluates second derivatives of the function.
  SymmTensor hessian(const Vec3& X) const override;

protected:
  //! \brief Evaluates the scalar field function.
  Real evaluate(const Vec3& X) const override;
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
  //! \brief Construct from vector of fields.
  //! \param fields Fields to use
  FieldsFuncBase(const std::vector<Fields*>& fields);

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

  //! \brief Evaluates the field gradient at the givent point \b X.
  std::vector<Real> getGradient(const Vec3& X) const;

  //! \brief Evaluates the field gradient at the givent point \b X.
  std::vector<Real> getHessian(const Vec3& X) const;

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

  //! \brief Wrap existing fields.
  //! \param fields Fields to evaluate
  //!
  //! \details Takes ownership of given fields.
  VecFieldFunction(const std::vector<Fields*>& fields);

  //! \brief Empty destructor.
  virtual ~VecFieldFunction() {}

  //! \brief Sets the active patch.
  bool initPatch(size_t pIdx) override { return this->setPatch(pIdx); }

protected:
  //! \brief Evaluates the vectorial field function.
  Vec3 evaluate(const Vec3& X) const override;

  //! \brief Returns the gradient of the function as a 1D array.
  std::vector<Real> evalGradient(const Vec3& X) const override
  { return this->getGradient(X); }
  //! \brief Returns the hessian of the function as a 1D array.
  std::vector<Real> evalHessian(const Vec3& X) const override
  { return this->getHessian(X); }
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
  \brief A vector-valued spatial function, defined through gradient of a scalar field.
*/

class ScalarGradFieldFunction : public VecFunc, private FieldFuncScalarBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] level Time level to read for
  ScalarGradFieldFunction(const std::string& fileName,
                          const std::string& basisName,
                          const std::string& fieldName,
                          int level = 0);

  //! \brief Wrap existing fields.
  //! \param[in] fields Fields to evaluate
  //! \param[in] dim Dimension of fields
  //!
  //! \details Takes ownership of given fields.
  ScalarGradFieldFunction(const std::vector<Field*>& fields, int dim);

  //! \brief Empty destructor.
  virtual ~ScalarGradFieldFunction() {}

  //! \brief Sets the active patch.
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

protected:
  //! \brief Evaluates the vectorial field function.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief A vector-valued spatial function, defined through laplacian of a scalar field.
*/

class ScalarLaplacianFieldFunction : public VecFunc, private FieldFuncScalarBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] level Time level to read for
  ScalarLaplacianFieldFunction(const std::string& fileName,
                               const std::string& basisName,
                               const std::string& fieldName,
                               int level = 0);
  //! \brief Empty destructor.
  virtual ~ScalarLaplacianFieldFunction() {}

  //! \brief Sets the active patch.
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

protected:
  //! \brief Evaluates the vectorial field function.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief A tensor-valued spatial function, defined through gradient of a vector field.
*/

class VecGradFieldFunction : public TensorFunc, private FieldsFuncBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] level Time level to read for
  VecGradFieldFunction(const std::string& fileName,
                       const std::string& basisName,
                       const std::string& fieldName,
                       int level = 0);
  //! \brief Empty destructor.
  virtual ~VecGradFieldFunction() {}

  //! \brief Sets the active patch.
  virtual bool initPatch(size_t pIdx) { return this->setPatch(pIdx); }

protected:
  //! \brief Evaluates the tensorial field function.
  virtual Tensor evaluate(const Vec3& X) const;
};


/*!
  \brief A tensor-valued spatial function, defined through laplacian of a vector field.
  \details Laplacian refers to the hessian without cross terms.
*/
class VecLaplacianFieldFunction : public TensorFunc, private FieldsFuncBase
{
public:
  //! \brief The constructor creates a field from the provided HDF5-file.
  //! \param[in] fileName Name of the HDF5-file
  //! \param[in] basisName Name of the basis which the field values refer to
  //! \param[in] fieldName Name of the field in the HDF5-file
  //! \param[in] level Time level to read for
  VecLaplacianFieldFunction(const std::string& fileName,
                            const std::string& basisName,
                            const std::string& fieldName,
                            int level = 0);
  //! \brief Empty destructor.
  virtual ~VecLaplacianFieldFunction() {}

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
