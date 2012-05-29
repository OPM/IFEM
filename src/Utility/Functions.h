// $Id$
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
#include "ExprFunctions.h"
#include "Vec3.h"


/*!
  \brief A scalar-valued linear function.
*/

class LinearFunc : public ScalarFunc
{
  real scale; //!< Scaling factor

public:
  //! \brief Constructor initializing the function parameter.
  LinearFunc(real s = real(1)) : scale(s) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return scale == real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual real evaluate(const real& x) const { return scale*x; }
};


/*!
  \brief A scalar-valued ramp function, linear up to \a xmax.
*/

class RampFunc : public ScalarFunc
{
  real fval; //!< Max function value
  real xmax; //!< Function is linear from \a x = 0 to \a x = \a xmax

public:
  //! \brief Constructor initializing the function parameters.
  RampFunc(real f = real(1), real x = real(1)) : fval(f), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fval == real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual real evaluate(const real& x) const;
};


/*!
  \brief A scalar-valued dirac function.
*/

class DiracFunc : public ScalarFunc
{
  real amp;  //!< The amplitude of the dirac function
  real xmax; //!< Associated \a x value

public:
  //! \brief Constructor initializing the function parameters.
  DiracFunc(real a = real(1), real x = real(0)) : amp(a), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return amp == real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual real evaluate(const real& x) const;
};


/*!
  \brief A scalar-valued step function.
*/

class StepFunc : public ScalarFunc
{
  real amp;  //!< The amplitude of the step function
  real xmax; //!< Associated \a x value

public:
  //! \brief Constructor initializing the function parameters.
  StepFunc(real a, real x = real(0)) : amp(a), xmax(x) {}
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return amp == real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual real evaluate(const real& x) const;
};


/*!
  \brief A scalar-valued sinusoidal function.
*/

class SineFunc : public ScalarFunc
{
  real scale; //!< Amplitude of the sine function
  real freq;  //!< Angular frequency of the sine function
  real phase; //!< Phase shift of the sine function

public:
  //! \brief Constructor initializing the function parameters.
  SineFunc(real s = real(1), real f = real(1), real p = real(0))
    : scale(s), freq(f), phase(p) {}
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return scale == real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual real evaluate(const real& x) const;
};


/*!
  \brief A scalar-valued spatial function, constant in space and time.
*/

class ConstFunc : public RealFunc
{
  real fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  ConstFunc(real v) : fval(v) {}
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fval == real(0); }

protected:
  //! \brief Evaluates the constant function.
  virtual real evaluate(const Vec3&) const { return fval; }
};


/*!
  \brief A scalar-valued spatial function, constant in space, varying in time.
*/

class ConstTimeFunc : public RealFunc
{
  const ScalarFunc* tfunc; //!< The time dependent function value

public:
  //! \brief Constructor initializing the function value.
  ConstTimeFunc(const ScalarFunc* f) : tfunc(f) {}
  //! \brief The destructor frees the time function.
  virtual ~ConstTimeFunc() { delete tfunc; }
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return tfunc->isZero(); }

protected:
  //! \brief Evaluates the time-varying function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, varying in space and time.
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
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return sfunc->isZero() || tfunc->isZero(); }

protected:
  //! \brief Evaluates the space-time function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear in \a x.
*/

class LinearXFunc : public RealFunc
{
  real a; //!< The function derivative
  real b; //!< The function value at \a x = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearXFunc(real A, real B = real(0)) : a(A), b(B) {}
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return a == real(0) && b == real(0); }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear in \a y.
*/

class LinearYFunc : public RealFunc
{
  real a; //!< The function derivative
  real b; //!< The function value at \a y = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearYFunc(real A, real B = real(0)) : a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return a == real(0) && b == real(0); }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear in \a z.
*/

class LinearZFunc : public RealFunc
{
  real a; //!< The function derivative
  real b; //!< The function value at \a z = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearZFunc(real A, real B = real(0)) : a(A), b(B) {}
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return a == real(0) && b == real(0); }

protected:
  //! \brief Evaluates the linear function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a x.
*/

class QuadraticXFunc : public RealFunc
{
  real max; //!< Max value of function
  real a;   //!< First root where function is zero
  real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticXFunc(real MAX, real A, real B) : max(MAX), a(A), b(B) {}
  
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return max == real(0); }

protected:
  //! \brief Evaluates the quadratic function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a y.
*/

class QuadraticYFunc : public RealFunc
{
  real max; //!< Max value of function
  real a;   //!< First root where function is zero
  real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticYFunc(real MAX, real A, real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return max == real(0); }

protected:
  //! \brief Evaluates the quadratic function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a z.
*/

class QuadraticZFunc : public RealFunc
{
  real max; //!< Max value of function
  real a;   //!< First root where function is zero
  real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticZFunc(real MAX, real A, real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return max == real(0); }

protected:
  //! \brief Evaluates the quadratic function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, defining a rotation about the Z-axis.
  \details The time component of the function argument multiplied with the
  function parameter \a A, is interpreted as the angle of rotation (in radians)
  about the global Z-axis passing through the point \a x0, \a y0.
  The function then returns the translation in either X- or Y-direction
  (depending on the \a retX argument to the constructor) of the global point
  { \a X.x, \a X.y } corresponding to this rotation.
  \note If the function is passed a Vec3 object as argument (and not a Vec4),
  it will always return zero.
*/

class LinearRotZFunc : public RealFunc
{
  bool rX; //!< Flag telling whether to return the X- (true) or Y-component
  real A;  //!< Magnitude of the rotation
  real x0; //!< Global X-coordinate of rotation centre
  real y0; //!< Global Y-coordinate of rotation centre

public:
  //! \brief Constructor initializing the function parameters.
  LinearRotZFunc(bool retX, real a, real x_0 = real(0), real y_0 = real(0))
    : rX(retX), A(a), x0(x_0), y0(y_0) {}

protected:
  //! \brief Evaluates the rotation function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, step in \a x.
*/

class StepXFunc : public RealFunc
{
  real fv; //!< The non-zero function value
  real x0; //!< Function is zero for \a x < \a x0
  real x1; //!< Function is zero for \a x > \a x1

public:
  //! \brief Constructor initializing the function parameters.
  StepXFunc(real v, real X0 = real(0), real X1 = real(1))
    : fv(v), x0(X0), x1(X1) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fv == real(0); }

protected:
  //! \brief Evaluates the step function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, step in \a x and \a y.
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
    : fv(v), x0(X0), y0(Y0), x1(X1), y1(Y1) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fv == real(0); }

protected:
  //! \brief Evaluates the step function.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear interpolation.
*/

class Interpolate1D : public RealFunc
{
  std::vector<real> grid;   //!< The (1D) grid the data is associated with
  std::vector<real> values; //!< The (scalar) data values
  int               dir;    //!< In which direction to perform the interpolation

public:
  //! \brief The constructor initializes the function parameters from a file.
  //! \param[in] file Name of file to read grid and function data from
  //! \param[in] dir_ Coordinate direction of the spatial variation
  Interpolate1D(const char* file, int dir_);

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return values.empty(); }

protected:
  //! \brief Evaluates the function by interpolating the 1D grid.
  virtual real evaluate(const Vec3& X) const;
};


/*!
  \brief A vector-valued spatial function, constant in space and time.
*/

class ConstVecFunc : public VecFunc
{
  Vec3 fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  ConstVecFunc(const Vec3& v) : fval(v) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fval.isZero(0.0); }

protected:
  //! \brief Evaluates the constant function.
  virtual Vec3 evaluate(const Vec3&) const { return fval; }
};


namespace utl
{
  //! \brief Creates a scalar-valued function by parsing a character string.
  const RealFunc* parseRealFunc(char* cline, real A = real(1));

  //! \brief Creates a scalar-valued function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type Function definition type flag
  RealFunc* parseRealFunc(const std::string& func, const std::string& type);

  //! \brief Creates a vector-valued function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type Function defintion type flag
  VecFunc* parseVecFunc(const std::string& func, const std::string& type);

  //! \brief Creates a vector-valued function defining a surface traction.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type Function defintion type flag
  //! \param[in] dir Coordinate direction of the traction (0=normal direction)
  TractionFunc* parseTracFunc(const std::string& func,
			      const std::string& type, int dir);
}

#endif
