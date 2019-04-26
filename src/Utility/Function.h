// $Id$
//==============================================================================
//!
//! \file Function.h
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General functions with arbitrary argument and value type.
//!
//==============================================================================

#ifndef UTL_FUNCTION_H
#define UTL_FUNCTION_H

#include "Vec3.h"
#include <functional>
#include <cstddef>

class TensorFunc;
class STensorFunc;


namespace utl
{
  /*!
    \brief Base class for unary functions of arbitrary result and argument type.
  */

  template<class Arg, class Result>
  class Function : public std::unary_function<const Arg&,Result>
  {
  protected:
    //! \brief The constructor is protected to allow sub-class instances only.
    Function() {}

  public:
    //! \brief Empty destructor.
    virtual ~Function() {}

    //! \brief Returns whether the function is identically zero or not.
    virtual bool isZero() const { return false; }
    //! \brief Returns whether the function is time-independent or not.
    virtual bool isConstant() const { return true; }

  protected:
    //! \brief Evaluates the function for the argument \a x.
    virtual Result evaluate(const Arg& x) const = 0;

  public:
    //! \brief Operator returning the function value for the given argument.
    Result operator()(const Arg& x) const { return this->evaluate(x); }

    typedef Arg    Input;  //!< Input type
    typedef Result Output; //!< Output type
  };


  /*!
    \brief Base class for binary function of arbitrary result and argument type.
    \details The two arguments have to be of the same type.
  */

  template<class Arg, class Result>
  class Function2 : public std::binary_function<const Arg&,const Arg&,Result>
  {
  protected:
    //! \brief The constructor is protected to allow sub-class instances only.
    Function2() {}

  public:
    //! \brief Empty destructor.
    virtual ~Function2() {}

    //! \brief Returns whether the function is identically zero or not.
    virtual bool isZero() const { return false; }

  protected:
    //! \brief Evaluates the function for the arguments \a x and \a y.
    virtual Result evaluate(const Arg& x, const Arg& y) const = 0;

  public:
    //! \brief Operator returning the function value for the given arguments.
    Result operator()(const Arg& x, const Arg& y) const
    { return this->evaluate(x,y); }

    typedef Arg    Input;  //!< Input type
    typedef Result Output; //!< Output type
  };


  /*!
    \brief Base class for unary spatial function of arbitrary result type.
    \details Includes interfaces for evaluation of first and second derivates.
  */

  template<class Result>
  class SpatialFunction : public Function<Vec3,Result>
  {
  protected:
    //! \brief The constructor is protected to allow sub-class instances only.
    explicit SpatialFunction(const Result& val) : zero(val) {}

  public:
    //! \brief Empty destructor.
    virtual ~SpatialFunction() {}

    //! \brief Returns a first-derivative of the function.
    virtual Result deriv(const Vec3&, int) const { return zero; }
    //! \brief Returns a second-derivative of the function.
    virtual Result dderiv(const Vec3&, int, int) const { return zero; }

  protected:
    Result zero; //!< Return value for default implementations of derivatives
  };
}


/*!
  \brief Scalar-valued unary function of a scalar value.
*/

class ScalarFunc : public utl::Function<Real,Real>
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  ScalarFunc() {}

public:
  //! \brief Empty destructor.
  virtual ~ScalarFunc() {}

  //! \brief Returns the first-derivative of the function.
  virtual Real deriv(Real) const { return Real(0); }
};


/*!
  \brief Base class for unary spatial functions of arbitrary result type.
  \details Includes an interface for returning the function value as an array.
*/

class FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  FunctionBase() : ncmp(1) {}

public:
  //! \brief Empty destructor.
  virtual ~FunctionBase() {}

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3&) const = 0;

  //! \brief Returns the number of components of the return value.
  size_t dim() const { return ncmp; }

  //! \brief Sets the active patch.
  virtual bool initPatch(size_t) { return true; }

  //! \brief Checks if a specified point is within the function domain.
  virtual bool inDomain(const Vec3&) const { return true; }

protected:
  size_t ncmp; //!< Number of components in the return value
};


/*!
  \brief Scalar-valued unary function of a spatial point.
*/

class RealFunc : public utl::SpatialFunction<Real>, public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  RealFunc() : utl::SpatialFunction<Real>(Real(0)) {}

public:
  //! \brief Empty destructor.
  virtual ~RealFunc() {}

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3& X) const
  {
    return std::vector<Real>(1,this->evaluate(X));
  }
};


/*!
  \brief Vector-valued unary function of a spatial point.
*/

class VecFunc : public utl::SpatialFunction<Vec3>, public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  explicit VecFunc(size_t n = 3) : utl::SpatialFunction<Vec3>(Vec3())
  {
    ncmp = n;
  }

public:
  //! \brief Empty destructor.
  virtual ~VecFunc() {}

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3& X) const
  {
    return this->evaluate(X).vec(ncmp);
  }
};


/*!
  \brief Vector-valued binary function of a spatial point and normal vector.
*/

class TractionFunc : public utl::Function2<Vec3,Vec3>
{
public:
  //! \brief Returns whether the traction is always normal to the face or not.
  virtual bool isNormalPressure() const { return false; }
};


/*!
  \brief Traction field based on a given pressure function.
*/

class PressureField : public TractionFunc
{
  const RealFunc* pressure; //!< Scalar field to derive the traction field from
  char            pdir;     //!< The global pressure direction (0...3)

public:
  //! \brief Constructor initializing a constant pressure field.
  //! \param[in] p The constant pressure value
  //! \param[in] dir The global direction the pressure is acting in
  PressureField(Real p, int dir = 0);
  //! \brief Constructor initializing the scalar pressure field function.
  //! \param[in] p The scalar field defining the spatial pressure distribution
  //! \param[in] dir The global direction the pressure is acting in
  PressureField(const RealFunc* p, int dir = 0) : pressure(p), pdir(dir) {}
  //! \brief The destructor frees the scalar field function.
  virtual ~PressureField() { delete pressure; }

  //! \brief Returns whether the traction is always normal to the face or not.
  virtual bool isNormalPressure() const { return pdir < 1 || pdir > 3; }
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return pressure ? pressure->isZero() : true; }

protected:
  //! \brief Evaluates the traction at point \a x and surface normal \a n.
  virtual Vec3 evaluate(const Vec3& x, const Vec3& n) const;
};


/*!
  \brief Traction field based on a given stress tensor function.
*/

class TractionField : public TractionFunc
{
  const STensorFunc* sigma; //!< Symmetric tensor field to derive tractions from
  const TensorFunc* sigmaN; //!< Tensor field to derive tractions from

public:
  //! \brief Constructor initializing the symmetric tensor function pointer.
  explicit TractionField(const STensorFunc& field);
  //! \brief Constructor initializing the tensor function pointer.
  explicit TractionField(const TensorFunc& field);
  //! \brief Empty destructor.
  virtual ~TractionField() {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const;

protected:
  //! \brief Evaluates the traction at point \a x and surface normal \a n.
  virtual Vec3 evaluate(const Vec3& x, const Vec3& n) const;
};

#endif
