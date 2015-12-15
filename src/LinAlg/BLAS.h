// $Id$
//==============================================================================
//!
//! \file BLAS.h
//!
//! \date Jan 10 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief BLAS support for various platforms.
//!
//==============================================================================

#ifndef UTL_BLAS_H
#define UTL_BLAS_H

#if defined(USE_MKL)
#define HAS_BLAS 2
#include <mkl_cblas.h>
#elif defined(USE_ACCELERATE)
#define HAS_BLAS 3
#include <Accelerate/Accelerate.h>
#elif defined(USE_CBLAS)
#define HAS_BLAS 1
extern "C"
{
#include <cblas.h>
}
#endif

#endif
