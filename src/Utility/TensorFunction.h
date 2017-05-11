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

class TensorFunc : public utl::SpatialFunction<Tensor>
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  TensorFunc(size_t n = 0) : utl::SpatialFunction<Tensor>(Tensor(n)) {}
public:
  //! \brief Empty destructor.
  virtual ~TensorFunc() {}
};


/*!
  \brief Symmetric tensor-valued unary function of a spatial point.
*/

class STensorFunc : public utl::SpatialFunction<SymmTensor>
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  STensorFunc(size_t n = 0) : utl::SpatialFunction<SymmTensor>(SymmTensor(n)) {}
public:
  //! \brief Empty destructor.
  virtual ~STensorFunc() {}
};

#endif
