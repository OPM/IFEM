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

class FunctionSum : public FunctionBase
{
  typedef std::pair<FunctionBase*,double> WeightedFunc; //!< Convenience type

protected:
  //! \brief Constructor to allow empty initialization in subclasses.
  FunctionSum() {}

public:
  //! \brief The constructor specifies the first function to sum.
  //! \param[in] f Pointer to a function to sum
  //! \param[in] w Weighting factor. If negative, take max value instead.
  FunctionSum(FunctionBase* f, double w = 1.0) { this->add(f,w); }

  //! \brief Empty destructor.
  virtual ~FunctionSum() {}

  //! \brief Adds a function to the list of functions to sum.
  //! \param[in] f Pointer to a function to sum
  //! \param[in] w Weighting factor. If negative, take max value instead.
  bool add(FunctionBase* f, double w = 1.0);

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const;

  //! \brief Checks if a specified point is within the function domain.
  virtual bool inDomain(const Vec3& X) const;
  //! \brief Returns \e true if current patch is affected by this function.
  virtual bool initPatch(size_t idx);

  //! \brief Returns the function value as an array.
  virtual std::vector<double> getValue(const Vec3& X) const;
  //! \brief Returns a representative scalar equivalent of the function value.
  virtual double getScalarValue(const Vec3& X) const;

private:
  std::vector<WeightedFunc> comps; //!< List of weighted functions to sum
};

#endif
