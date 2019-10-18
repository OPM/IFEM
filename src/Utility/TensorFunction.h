// $Id$
//==============================================================================
//!
//! \file TensorFunction.h
//!
//! \date May 10 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Spatial tensor-valued functions.
//!
//==============================================================================

#ifndef _TENSOR_FUNCTION_H
#define _TENSOR_FUNCTION_H

#include "Function.h"
#include "Tensor.h"


/*!
  \brief Tensor-valued unary function of a spatial point.
*/

class TensorFunc : public utl::SpatialFunction<Tensor>, public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  explicit TensorFunc(size_t n = 0) : utl::SpatialFunction<Tensor>(Tensor(n))
  {
    ncmp = zero.size();
  }

public:
  //! \brief Empty destructor.
  virtual ~TensorFunc() {}

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const { return 3; }

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3& X) const
  {
    return this->evaluate(X);
  }

  //! \brief Returns a representative scalar equivalent of the function value.
  virtual Real getScalarValue(const Vec3& X) const
  {
    return this->evaluate(X).trace();
  }
};


/*!
  \brief Symmetric tensor-valued unary function of a spatial point.
*/

class STensorFunc : public utl::SpatialFunction<SymmTensor>,
                    public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  STensorFunc(size_t n = 0, bool with33 = false)
    : utl::SpatialFunction<SymmTensor>(SymmTensor(n,with33))
  {
    ncmp = zero.size();
  }

public:
  //! \brief Empty destructor.
  virtual ~STensorFunc() {}

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const { return 3; }

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3& X) const
  {
    return this->evaluate(X);
  }

  //! \brief Returns a representative scalar equivalent of the function value.
  virtual Real getScalarValue(const Vec3& X) const
  {
    return this->evaluate(X).trace();
  }
};

#endif
