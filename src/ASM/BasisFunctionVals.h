// $Id$
//==============================================================================
//!
//! \file BasisFunctionVals.h
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function values container.
//!
//==============================================================================

#ifndef _BASIS_FUNCTION_VALS_H
#define _BASIS_FUNCTION_VALS_H

#include "MatVec.h"


/*!
  \brief Struct holding basis function values and derivatives.
*/

struct BasisFunctionVals
{
  Vector     N;    //!< Basis function values
  Matrix    dNdu;  //!< Basis function derivatives
  Matrix3D d2Ndu2; //!< Second order derivatives of basis functions
  Matrix4D d3Ndu3; //!< Third order derivatives of basis functions
};


//! \brief Convenience type alias
using BasisValuesPtrs = std::vector<const BasisFunctionVals*>;


//! \brief Utility class holding a vector of basis function values.
//! \details This is a helper class that allows us to avoid pointer
//!          semantics in code, while allowing to pass as
//!          pointers into mapping functions. This is necessary
//!          to avoid additional copying in code using the
//!          basis function cache.
//!          While it overrides most common vector operations,
//!          be careful with advanced insertions.
class BasisValues : public std::vector<BasisFunctionVals>
{
public:
  //! \brief Constructor resizing to a given size.
  BasisValues(size_t size) : std::vector<BasisFunctionVals>(size)
  {
    pointers.reserve(this->size());
    for (const BasisFunctionVals& val : *this)
      pointers.push_back(&val);
  }

   //! \brief Destructor.
   //! \details Need to manually implement to call parent dtor.
  ~BasisValues()
  {
    this->clear();
  }

  //! \brief Clears the container.
  void clear() noexcept
  {
    this->std::vector<BasisFunctionVals>::clear();
    pointers.clear();
  }

  //! \brief Push back for reference.
  void push_back(const BasisFunctionVals& v)
  {
    this->std::vector<BasisFunctionVals>::push_back(v);
    pointers.push_back(&this->back());
  }

  //! \brief Push back for r-reference.
  void push_back(BasisFunctionVals&& v)
  {
    this->std::vector<BasisFunctionVals>::push_back(v);
    pointers.push_back(&this->back());
  }

  //! \brief Cast to a vector with pointers to our values.
  operator const BasisValuesPtrs& () const { return pointers; }

private:
  BasisValuesPtrs pointers; //!< Vector of pointers to our values
};

#endif
