// $Id$
//==============================================================================
//!
//! \file Chebyshev.h
//!
//! \date Jul 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Evaluation of Chebyshev polynomials.
//!
//==============================================================================

#ifndef _CHEBYSHEV_H
#define _CHEBYSHEV_H

#include "Function.h"
#include "TensorFunction.h"

#include <array>
#include <functional>
#include <memory>


/*!
  \brief Evaluation of Chebyshev polynomials.
*/

namespace Chebyshev
{
  //! \brief Evaluates a 1D Chebyshev polynomial of first kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  Real evalPol1(int polnum, Real xi);
  //! \brief Evaluates a 1D Chebyshev polynomial of second kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  Real evalPol2(int polnum, Real xi);

  //! \brief Evaluates the first derivative of a 1D Chebyshev polynomial of first kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  Real evalDer1(int polnum, Real xi);
  //! \brief Evaluates the first derivative of a 1D Chebyshev polynomial of second kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  Real evalDer2(int polnum, Real xi);

  //! \brief Evaluates the second derivative of a 1D Chebyshev polynomial of first kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  Real eval2Der1(int polnum, Real xi);
}


/*!
  \brief A scalar-valued spatial function, chebyshev polynomials.
*/

class ChebyshevFunc : public RealFunc
{
  std::vector<Real> coefs; //!< Function coefficients
  std::array<int,3> n;     //!< Number of coefficients
  std::array<std::array<double,2>,3> domain; //!< Physical domain

public:
  //! \brief The constructor initializes the function parameters from a file.
  //! \param[in] input Name of file or string to read \ref coefs from
  //! \param[in] file True if input is a file name
  ChebyshevFunc(const std::string& input, bool file);

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return coefs.empty(); }

  //! \brief Returns a const reference to the coefficients.
  const std::vector<Real>& getCoefs() const { return coefs; }
  //! \brief Returns the number of polynomials in each parameter direction.
  const std::array<int,3>& getSize() const { return n; }

  //! \brief Returns a first-derivative of the function.
  Real deriv(const Vec3& X, int c) const override;

  //! \brief Returns a second-derivative of the function.
  Real dderiv(const Vec3& X, int c, int c2) const override;

protected:
  //! \brief Evaluates the function at point \a X.
  Real evaluate(const Vec3& X) const override;

private:
  //! \brief Struct defining function and weight in one direction.
  struct Func
  {
    std::function<Real(int, Real)> f; //!< Function to evaluate
    double w; //!< Weight
  };

  //! \brief Performs the actual tensor-product evaluation.
  Real evaluateTP (const Vec3& X, const std::array<Func,3>& funcs) const;

  //! \brief Reads input from a stream.
  //! \param in Stream to read from
  void read(std::istream& in);
};


/*!
  \brief A vector-valued spatial function, chebyshev polynomials.
*/

class ChebyshevVecFunc : public VecFunc
{
  std::array<std::unique_ptr<ChebyshevFunc>,3> f; //!< Functions
  bool secondDer; //!< True to take second derivatives

public:
  //! \brief The constructor initializes the function parameters from files.
  //! \param[in] input Name of files or strings to read coefs from
  //! \param[in] file True if input is file names
  //! \param[in] second True to take second derivatives
  ChebyshevVecFunc(const std::vector<std::string>& input,
                   bool file, bool second = false);

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return !f[0] || f[0]->isZero(); }

protected:
  //! \brief Evaluates the function at point \a X.
  Vec3 evaluate(const Vec3& X) const override;
  //! \brief Returns the gradient of the function as a 1D array.
  std::vector<Real> evalGradient(const Vec3&) const override;
  //! \brief Returns the hessian of the function as a 1D array.
  std::vector<Real> evalHessian(const Vec3& X) const override;
};


/*!
  \brief A tensor-valued spatial function, chebyshev polynomials.

  \details If 2 or 3 functions: Take the derivative of the interpolants.
           If 4 or 9 functions: Each component has their own interpolant.
*/

class ChebyshevTensorFunc : public TensorFunc
{
  std::array<std::unique_ptr<ChebyshevVecFunc>,3> f; //!< Array of vector components

public:
  //! \brief The constructor initializes the function parameters from files.
  //! \param[in] input Name of files or strings to read coefs from
  //! \param[in] file True if input is file names
  //! \param[in] second True to take second derivatives
  ChebyshevTensorFunc(const std::vector<std::string>& input,
                      bool file, bool second);

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return !(f[0] || f[1] || f[2]); }

protected:
  //! \brief Evaluates the function at point \a X.
  virtual Tensor evaluate(const Vec3& X) const;
};

#endif
