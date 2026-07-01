// $Id$
//==============================================================================
//!
//! \file FunctionSum.h
//!
//! \date Apr 16 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unary spatial function as a sum of other spatial functions.
//!
//==============================================================================

#ifndef _FUNCTION_SUM_H
#define _FUNCTION_SUM_H

#include "Function.h"


/*!
  \brief Unary spatial function as a sum of other spatial functions.
*/

class FunctionSum : public virtual FunctionBase
{
protected:
  //! \brief Constructor to allow empty initialization in subclasses.
  explicit FunctionSum(bool ownFn = false) : ownFunc(ownFn) {}

public:
  //! \brief The constructor specifies the first function to sum.
  //! \param[in] f Pointer to a function to sum
  //! \param[in] w Weighting factor. If negative, take max value instead.
  explicit FunctionSum(FunctionBase* f, double w = 1.0) : ownFunc(false)
  {
    this->add(f,w);
  }

  //! \brief The destructor deletes the function components if flagged so.
  virtual ~FunctionSum();

  //! \brief Adds a function to the list of functions to sum.
  //! \param[in] f Pointer to a function to sum
  //! \param[in] w Weighting factor. If negative, take max value instead.
  bool add(FunctionBase* f, double w = 1.0);

  //! \brief Returns the function type flag.
  unsigned char getType() const override;

  //! \brief Checks if a specified point is within the function domain.
  bool inDomain(const Vec3& X) const override;
  //! \brief Returns \e true if current patch is affected by this function.
  bool initPatch(size_t idx) override;

  //! \brief Returns the function value as an array.
  std::vector<double> getValue(const Vec3& X) const override;
  //! \brief Returns a representative scalar equivalent of the function value.
  double getScalarValue(const Vec3& X) const override;

private:
  using WeightedFunc = std::pair<FunctionBase*,double>; //!< Convenience type

  std::vector<WeightedFunc> comps; //!< List of weighted functions to sum

  bool ownFunc; //!< If \e true, the destructor deletes the function components
};


/*!
  \brief A real-valued spatial function as a sum of other spatial functions.
*/

class RealFuncSum : public FunctionSum, public RealFunc
{
protected:
  //! \brief The constructor is protected to allow subclass instances only.
  RealFuncSum() : FunctionSum(true) {}

  //! \brief Evaluates the function.
  double evaluate(const Vec3& X) const override
  {
    return this->FunctionSum::getScalarValue(X);
  }

public:
  //! \copydoc FunctionSum::getScalarValue()
  double getScalarValue(const Vec3& X) const override
  {
    return this->FunctionSum::getScalarValue(X);
  }

  //! \copydoc FunctionSum::getValue()
  std::vector<double> getValue(const Vec3& X) const override
  {
    return this->FunctionSum::getValue(X);
  }

  //! \copydoc FunctionSum::getType()
  unsigned char getType() const override
  {
    return this->FunctionSum::getType();
  }
};


/*!
  \brief A sum of spatial dirac functions.
*/

class DiracSum : public RealFuncSum
{
public:
  //! \brief The constructor reads the function definition from a string.
  //! \param[in] input The string to parse for function parameters
  //! \param[in] tol Radius of non-zero-valued function domains
  //! \param[in] nsd Number of spatial dimensions
  DiracSum(const char* input, double tol, int nsd);
};

#endif
