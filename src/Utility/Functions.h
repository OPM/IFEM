// $Id: Functions.h,v 1.7 2011-02-08 12:55:33 rho Exp $
//==============================================================================
//!
//! \file Functions.h
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Specific function implementations.
//!
//==============================================================================

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "Function.h"
#include <math.h>


/*!
  \brief A linear scalar function.
*/

class LinearFunc : public ScalarFunc
{
  real scale; //!< Scaling factor

public:
  //! \brief Constructor initializing the function parameter.
  LinearFunc(real s = real(1)) { scale = s; }

protected:
  //! \brief Evaluates the scalar function.
  virtual real evaluate(const real& x) const { return scale*x; }
};


/*!
  \brief A sinusoidal scalar function.
*/

class SineFunc : public ScalarFunc
{
  real scale; //!< Amplitude of the sine function
  real freq;  //!< Angular frequency of the sine function
  real phase; //!< Phase shift of the sine function

public:
  //! \brief Constructor initializing the function parameters.
  SineFunc(real s = real(1), real f = real(1), real p = real(0))
  { scale = s; freq = f; phase = p; }

protected:
  //! \brief Evaluates the scalar function.
  virtual real evaluate(const real& x) const { return scale*sin(freq*x+phase); }
};


/*!
  \brief A scalar function, constant in space and time.
*/

class ConstFunc : public RealFunc
{
  real fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  ConstFunc(real v) { fval = v; }

protected:
  //! \brief Evaluates the constant function.
  virtual real evaluate(const Vec3&) const { return fval; }
};


/*!
  \brief A scalar function, constant in space but varying in time.
*/

class ConstTimeFunc : public RealFunc
{
  const ScalarFunc* tfunc; //!< The time dependent function value

public:
  //! \brief Constructor initializing the function value.
  ConstTimeFunc(const ScalarFunc* f) : tfunc(f) {}
  //! \brief The destructor frees the time function.
  virtual ~ConstTimeFunc() { delete tfunc; }

protected:
  //! \brief Evaluates the time-varying function.
  virtual real evaluate(const Vec3& x) const;
};


/*!
  \brief A scalar function, varying in space and time.
  \details The function value is defined as a product between one
  space-dependent component and one time-dependent component.
*/

class SpaceTimeFunc : public RealFunc
{
  const RealFunc*   sfunc; //!< The space-dependent term
  const ScalarFunc* tfunc; //!< The time-dependent term

public:
  //! \brief Constructor initializing the function terms.
  SpaceTimeFunc(const RealFunc* s, const ScalarFunc* t) : sfunc(s), tfunc(t) {}
  //! \brief The destructor frees the space and time functions.
  virtual ~SpaceTimeFunc() { delete sfunc; delete tfunc; }

protected:
  //! \brief Evaluates the space-time function.
  virtual real evaluate(const Vec3& x) const;
};


/*!
  \brief A scalar function, linear in \a t up to \a Tinit.
*/

class LinearTinitFunc : public RealFunc
{
  real value; //!< Max function value
  real Tinit; //!< Function is linear from 0 to t = Tinit

public:
  //! \brief Constructor initializing the function parameters.
  LinearTinitFunc(real value_, real Tinit_) { value = value_; Tinit = Tinit_; }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& x) const;
};


/*!
  \brief A scalar function, linear in \a x.
*/

class LinearXFunc : public RealFunc
{
  real a; //!< The function derivative
  real b; //!< The function value at \a x = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearXFunc(real A, real B = real(0)) { a = A; b = B; }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, linear in \a y.
*/

class LinearYFunc : public RealFunc
{
  real a; //!< The function derivative
  real b; //!< The function value at \a y = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearYFunc(real A, real B = real(0)) { a = A; b = B; }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, linear in \a z.
*/

class LinearZFunc : public RealFunc
{
  real a; //!< The function derivative
  real b; //!< The function value at \a z = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearZFunc(real A, real B = real(0)) { a = A; b = B; }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, quadratic in \a x.
*/

class QuadraticXFunc : public RealFunc
{
  real max;   // Max value of function
  real a, b;  // Roots where function is \a 0

 public:
  //! \brief Constructor initializing the function parameters.
  QuadraticXFunc(real MAX, real A, real B) { max = MAX; a = A; b = B; }

 protected:
  //! \brief Evaluates the quadratic function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, quadratic in \a y.
*/

class QuadraticYFunc : public RealFunc
{
  real max;   // Max value of function
  real a, b;  // Roots where function is \a 0

 public:
  //! \brief Constructor initializing the function parameters.
  QuadraticYFunc(real MAX, real A, real B) { max = MAX; a = A; b = B; }

 protected:
  //! \brief Evaluates the quadratic function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, quadratic in \a x.
*/

class QuadraticZFunc : public RealFunc
{
  real max;   // Max value of function
  real a, b;  // Roots where function is \a 0

 public:
  //! \brief Constructor initializing the function parameters.
  QuadraticZFunc(real MAX, real A, real B) { max = MAX; a = A; b = B; }

 protected:
  //! \brief Evaluates the quadratic function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, defining a linear rotation about the global Z-axis.
  \details The time component of the function argument multiplied with the
  function parameter \a A, is interpreted as the angle of rotation (in radians)
  about the Z-axis passing through the point \a x0, \a y0.
  The function then returns the translation in either \a x or \a y direction
  (depending on the \a retX argument to the constructor) of the global point
  { \a X.x, \a X.y } corresponding to this rotation.
  \note If the function is passed a Vec3 object as argument (and not a Vec4),
  it will always return zero.
*/

class LinearRotZFunc : public RealFunc
{
  bool rX; //!< Flag telling whether to return the X (true) or Y component
  real A;  //!< Magnitude of the rotation
  real x0; //!< Global x-coordinate of rotation centre
  real y0; //!< Global y-coordinate of rotation centre

public:
  //! \brief Constructor initializing the function parameters.
  LinearRotZFunc(bool retX, real a, real x_0 = real(0), real y_0 = real(0))
  { rX = retX; A = a; x0 = x_0; y0 = y_0; }

protected:
  //! \brief Evaluates the rotation function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, step in \a x.
*/

class StepXFunc : public RealFunc
{
  real fv; //!< The non-zero function value
  real x0; //!< Function is zero for \a x < \a x0
  real x1; //!< Function is zero for \a x > \a x1

public:
  //! \brief Constructor initializing the function parameters.
  StepXFunc(real v, real X0 = real(0), real X1 = real(1))
  { fv = v; x0 = X0; x1 = X1; }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar function, step in \a x and \a y.
*/

class StepXYFunc : public RealFunc
{
  real fv; //!< The non-zero function value
  real x0; //!< Function is zero for \a x < \a x0
  real y0; //!< Function is zero for \a y < \a y0
  real x1; //!< Function is zero for \a x > \a x1
  real y1; //!< Function is zero for \a y > \a y1

public:
  //! \brief Constructor initializing the function parameters.
  StepXYFunc(real v,
	     real X1 = real(1), real Y1 = real(1),
	     real X0 = real(-1), real Y0 = real(-1))
  { fv = v; x0 = X0; y0 = Y0; x1 = X1; y1 = Y1; }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


namespace utl
{
  //! \brief Creates a real function by parsing data from a character string.
  const RealFunc* parseRealFunc(char* cline, real A = real(0));
}

#endif
