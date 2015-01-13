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
#include "FieldFunctions.h"
#include "Vec3.h"


/*!
  \brief A scalar-valued constant function.
*/

class ConstantFunc : public ScalarFunc
{
  Real fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  ConstantFunc(Real v) : fval(v) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fval == Real(0); }

protected:
  //! \brief Evaluates the constant function.
  virtual Real evaluate(const Real&) const { return fval; }
};


/*!
  \brief A scalar-valued linear function.
*/

class LinearFunc : public ScalarFunc
{
  Real scale; //!< Scaling factor

public:
  //! \brief Constructor initializing the function parameter.
  LinearFunc(Real s = Real(1)) : scale(s) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return scale == Real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual Real evaluate(const Real& x) const { return scale*x; }
};


/*!
  \brief A scalar-valued ramp function, linear up to \a xmax.
*/

class RampFunc : public ScalarFunc
{
  Real fval; //!< Max function value
  Real xmax; //!< Function is linear from \a x = 0 to \a x = \a xmax

public:
  //! \brief Constructor initializing the function parameters.
  RampFunc(Real f = Real(1), Real x = Real(1)) : fval(f), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fval == Real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual Real evaluate(const Real& x) const;
};


/*!
  \brief A scalar-valued dirac function.
*/

class DiracFunc : public ScalarFunc
{
  Real amp;  //!< The amplitude of the dirac function
  Real xmax; //!< Associated \a x value

public:
  //! \brief Constructor initializing the function parameters.
  DiracFunc(Real a = Real(1), Real x = Real(0)) : amp(a), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return amp == Real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual Real evaluate(const Real& x) const;
};


/*!
  \brief A scalar-valued step function.
*/

class StepFunc : public ScalarFunc
{
  Real amp;  //!< The amplitude of the step function
  Real xmax; //!< Associated \a x value

public:
  //! \brief Constructor initializing the function parameters.
  StepFunc(Real a, Real x = Real(0)) : amp(a), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return amp == Real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual Real evaluate(const Real& x) const;
};


/*!
  \brief A scalar-valued sinusoidal function.
*/

class SineFunc : public ScalarFunc
{
  Real scale; //!< Amplitude of the sine function
  Real freq;  //!< Angular frequency of the sine function
  Real phase; //!< Phase shift of the sine function

public:
  //! \brief Constructor initializing the function parameters.
  SineFunc(Real s = Real(1), Real f = Real(1), Real p = Real(0))
    : scale(s), freq(f), phase(p) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return scale == Real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  virtual Real evaluate(const Real& x) const;
};


/*!
  \brief A scalar-valued spatial function, constant in space and time.
*/

class ConstFunc : public RealFunc
{
  Real fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  ConstFunc(Real v) : fval(v) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fval == Real(0); }

protected:
  //! \brief Evaluates the constant function.
  virtual Real evaluate(const Vec3&) const { return fval; }
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
  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return tfunc->isZero(); }

protected:
  //! \brief Evaluates the time-varying function.
  virtual Real evaluate(const Vec3& X) const;
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
  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return this->isZero(); }

protected:
  //! \brief Evaluates the space-time function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear in \a x.
*/

class LinearXFunc : public RealFunc
{
  Real a; //!< The function derivative
  Real b; //!< The function value at \a x = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearXFunc(Real A, Real B = Real(0)) : a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return a == Real(0) && b == Real(0); }

protected:
  //! \brief Evaluates the linear function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear in \a y.
*/

class LinearYFunc : public RealFunc
{
  Real a; //!< The function derivative
  Real b; //!< The function value at \a y = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearYFunc(Real A, Real B = Real(0)) : a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return a == Real(0) && b == Real(0); }

protected:
  //! \brief Evaluates the linear function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear in \a z.
*/

class LinearZFunc : public RealFunc
{
  Real a; //!< The function derivative
  Real b; //!< The function value at \a z = 0

public:
  //! \brief Constructor initializing the function parameters.
  LinearZFunc(Real A, Real B = Real(0)) : a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return a == Real(0) && b == Real(0); }

protected:
  //! \brief Evaluates the linear function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a x.
*/

class QuadraticXFunc : public RealFunc
{
  Real max; //!< Max value of function
  Real a;   //!< First root where function is zero
  Real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticXFunc(Real MAX, Real A, Real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return max == Real(0); }

protected:
  //! \brief Evaluates the quadratic function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a y.
*/

class QuadraticYFunc : public RealFunc
{
  Real max; //!< Max value of function
  Real a;   //!< First root where function is zero
  Real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticYFunc(Real MAX, Real A, Real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return max == Real(0); }

protected:
  //! \brief Evaluates the quadratic function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a z.
*/

class QuadraticZFunc : public RealFunc
{
  Real max; //!< Max value of function
  Real a;   //!< First root where function is zero
  Real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticZFunc(Real MAX, Real A, Real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return max == Real(0); }

protected:
  //! \brief Evaluates the quadratic function.
  virtual Real evaluate(const Vec3& X) const;
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
  Real A;  //!< Magnitude of the rotation
  Real x0; //!< Global X-coordinate of rotation centre
  Real y0; //!< Global Y-coordinate of rotation centre

public:
  //! \brief Constructor initializing the function parameters.
  LinearRotZFunc(bool retX, Real a, Real x_0 = Real(0), Real y_0 = Real(0))
    : rX(retX), A(a), x0(x_0), y0(y_0) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return A == Real(0); }
  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return A == Real(0); }

protected:
  //! \brief Evaluates the rotation function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, step in \a x.
*/

class StepXFunc : public RealFunc
{
  Real fv; //!< The non-zero function value
  Real x0; //!< Function is zero for \a x < \a x0
  Real x1; //!< Function is zero for \a x > \a x1

public:
  //! \brief Constructor initializing the function parameters.
  StepXFunc(Real v, Real X0 = Real(0), Real X1 = Real(1))
    : fv(v), x0(X0), x1(X1) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fv == Real(0); }

protected:
  //! \brief Evaluates the step function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, step in \a x and \a y.
*/

class StepXYFunc : public RealFunc
{
  Real fv; //!< The non-zero function value
  Real x0; //!< Function is zero for \a x < \a x0
  Real y0; //!< Function is zero for \a y < \a y0
  Real x1; //!< Function is zero for \a x > \a x1
  Real y1; //!< Function is zero for \a y > \a y1

public:
  //! \brief Constructor initializing the function parameters.
  StepXYFunc(Real v,
	     Real X1 = Real(1), Real Y1 = Real(1),
	     Real X0 = Real(-1), Real Y0 = Real(-1))
    : fv(v), x0(X0), y0(Y0), x1(X1), y1(Y1) {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return fv == Real(0); }

protected:
  //! \brief Evaluates the step function.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A scalar-valued spatial function, linear interpolation.
*/

class Interpolate1D : public RealFunc
{
  std::vector<Real> grid;   //!< The (1D) grid the data is associated with
  std::vector<Real> values; //!< The (scalar) data values
  int               dir;    //!< In which direction to perform the interpolation
  Real              time;   //!< Ramp-up time

public:
  //! \brief The constructor initializes the function parameters from a file.
  //! \param[in] file Name of file to read grid and function data from
  //! \param[in] dir_ Coordinate direction of the spatial variation
  //! \param[in] col Which column of the file to read function values from
  //! \param[in] ramp Ramp-up time
  Interpolate1D(const char* file, int dir_, int col = 2, Real ramp = Real(0));

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return values.empty(); }
  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return time <= Real(0) || grid.size() < 2; }

protected:
  //! \brief Evaluates the function by interpolating the 1D grid.
  virtual Real evaluate(const Vec3& X) const;
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
  //! \brief Creates a time function by parsing a character string.
  const ScalarFunc* parseTimeFunc(const char* type, char* cline = NULL,
				  Real C = Real(1));

  //! \brief Creates a scalar-valued function by parsing a character string.
  const RealFunc* parseRealFunc(char* cline, Real A = Real(1));

  //! \brief Creates a time function by parsing a character string.
  ScalarFunc* parseTimeFunc(const char* func, const std::string& type);

  //! \brief Creates a scalar-valued function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type Function definition type flag
  RealFunc* parseRealFunc(const std::string& func, const std::string& type);

  //! \brief Creates a vector-valued function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type Function defintion type flag
  //! \param[in] variables Variable definition for expression functions
  VecFunc* parseVecFunc(const std::string& func, const std::string& type,
                        const std::string& variables="");

  //! \brief Creates a vector-valued function defining a surface traction.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type Function defintion type flag
  //! \param[in] dir Coordinate direction of the traction (0=normal direction)
  TractionFunc* parseTracFunc(const std::string& func,
			      const std::string& type, int dir);
}

#endif
