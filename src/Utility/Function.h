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

#include <functional>
#include <cstddef>


namespace utl
{
  /*!
    \brief Base class for unary function of arbitrary result and argument type.
  */

  template<class Arg, class Result>
  class Function : public std::unary_function<const Arg&,Result>
  {
  protected:
    //! \brief Empty constructor.
    Function() {}
  public:
    //! \brief Empty destructor.
    virtual ~Function() {}

    //! \brief Returns whether the function is identically zero or not.
    virtual bool isZero() const { return false; }

  protected:
    //! \brief Evaluates the function for the argument \a x.
    virtual Result evaluate(const Arg& x) const = 0;

  public:
    //! \brief Operator returning the function value for the given argument.
    Result operator()(const Arg& x) const { return this->evaluate(x); }
  };


  /*!
    \brief Base class for binary function of arbitrary result and argument type.
    \details The two arguments have to be of the same type.
  */

  template<class Arg, class Result>
  class Function2 : public std::binary_function<const Arg&,const Arg&,Result>
  {
  protected:
    //! \brief Empty constructor.
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
  };
}


class Vec3;
class Tensor;
class SymmTensor;


//! \brief Scalar-valued unary function of a scalar value.
typedef utl::Function<real,real> ScalarFunc;

//! \brief Scalar-valued unary function of a spatial point.
typedef utl::Function<Vec3,real> RealFunc;

//! \brief Vector-valued unary function of a spatial point.
typedef utl::Function<Vec3,Vec3> VecFunc;

//! \brief Tensor-valued unary function of a spatial point.
typedef utl::Function<Vec3,Tensor> TensorFunc;

//! \brief Symmetric tensor-valued unary function of a spatial point.
typedef utl::Function<Vec3,SymmTensor> STensorFunc;


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
  PressureField(real p, int dir = 0);
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
  TractionField(const STensorFunc& field) : sigma(&field), sigmaN(NULL) {}
  //! \brief Constructor initializing the tensor function pointer.
  TractionField(const TensorFunc& field) : sigma(NULL), sigmaN(&field) {}
  //! \brief Empty destructor.
  virtual ~TractionField() {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const;

protected:
  //! \brief Evaluates the traction at point \a x and surface normal \a n.
  virtual Vec3 evaluate(const Vec3& x, const Vec3& n) const;
};

#endif
