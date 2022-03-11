// $Id$
//==============================================================================
//!
//! \file PythonFunctions.h
//!
//! \date Sep 7 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Python function implementations.
//!
//==============================================================================

#ifndef _PYTHON_FUNCTIONS_H
#define _PYTHON_FUNCTIONS_H

#ifdef HAS_PYTHON

#include <pybind11/embed.h>

#include "Function.h"
#include "TensorFunction.h"

/*!
 * \brief RAII style class for managing the interpreter in use.
*/

class InterpreterRAII {
public:
  //! \brief The constructor allocates the interpreter.
  InterpreterRAII()
  {
    pybind11::initialize_interpreter();
  }

  //! \brief The destructor finalizes the interpreter.
  ~InterpreterRAII()
  {
    pybind11::finalize_interpreter();
  }
};


extern std::shared_ptr<InterpreterRAII> pyInterp;


/*!
  \brief Base class for python module functions.
*/

class __attribute__((visibility("hidden"))) PythonBaseFunc {
protected:
  //! \brief The constructor sets module to use and its parameters.
  //! \param module Name of python module to use
  //! \param params Parameters for python module in JSON format
  PythonBaseFunc(const char* module, const char* params);

  //! \brief Empty destructor.
  virtual ~PythonBaseFunc() = default;

  std::shared_ptr<InterpreterRAII> myInterp; //!< Ref-counted pointer to global interpreter

  pybind11::module myModule; //!< Imported python module
  pybind11::object myInstance; //!< Instance of module
  pybind11::object setParams; //!< Method to set parameters in instance
  pybind11::object eval; //!< Method to evaluate instance
};


/*!
  \brief A scalar-valued function, python module.
*/

class __attribute__((visibility("hidden"))) PythonFunc :
  public ScalarFunc, public PythonBaseFunc
{
public:
  //! \brief The constructor sets the module to use and its parameters.
  //! \param module Name of python module to use
  //! \param params Parameters for python module in JSON format
  PythonFunc(const char* function, const char* params) :
    PythonBaseFunc(function,params)
  {}

  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~PythonFunc() = default;

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

protected:
  //! \brief Evaluates the function expression.
  Real evaluate(const Real& X) const override;
};


/*!
  \brief A scalar-valued spatial function, python module.
*/

class __attribute__((visibility("hidden"))) PythonFunction :
  public RealFunc, public PythonBaseFunc
{
public:
  //! \brief The constructor sets the module to use and its parameters.
  //! \param module Name of python module to use
  //! \param params Parameters for python module in JSON format
  PythonFunction(const char* function, const char* params) :
    PythonBaseFunc(function,params)
  {}

  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~PythonFunction() = default;

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

protected:
  //! \brief Evaluates the function expression.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A vector-valued spatial function, python module.
*/

class __attribute__((visibility("hidden"))) PythonVecFunc :
  public VecFunc, public PythonBaseFunc
{
public:
  //! \brief The constructor sets the module to use and its parameters.
  //! \param module Name of python module to use
  //! \param params Parameters for python module in JSON format
  PythonVecFunc(const char* function, const char* params) :
    PythonBaseFunc(function,params)
  {}

  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~PythonVecFunc() = default;

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

protected:
  //! \brief Evaluates the function expression.
  Vec3 evaluate(const Vec3& X) const override;
};


/*!
  \brief A tensor-valued spatial function, python module.
*/

class __attribute__((visibility("hidden"))) PythonTensorFunc :
  public TensorFunc, public PythonBaseFunc
{
public:
  //! \brief The constructor sets the module to use and its parameters.
  //! \param n Number of spatial dimensions
  //! \param module Name of python module to use
  //! \param params Parameters for python module in JSON format
  PythonTensorFunc(const char* function, const char* params, size_t n = 0) :
    TensorFunc(n), PythonBaseFunc(function,params)
  {}

  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~PythonTensorFunc() = default;

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

protected:
  //! \brief Evaluates the function expression.
  Tensor evaluate(const Vec3& X) const override;
};


/*!
  \brief A symmetric tensor-valued spatial function, python module.
*/

class __attribute__((visibility("hidden"))) PythonSTensorFunc :
  public STensorFunc, public PythonBaseFunc
{
public:
  //! \brief The constructor sets the module to use and its parameters.
  //! \param n Number of spatial dimensions
  //! \param module Name of python module to use
  //! \param params Parameters for python module in JSON format
  PythonSTensorFunc(size_t n, const char* function, const char* params, bool with33 = false) :
    STensorFunc(n,with33), PythonBaseFunc(function,params), Result(n,with33)
  {}

  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~PythonSTensorFunc() = default;

  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

protected:
  //! \brief Evaluates the function expression.
  SymmTensor evaluate(const Vec3& X) const override;

  SymmTensor Result; //!< Cloned for results
};



#endif

#endif
