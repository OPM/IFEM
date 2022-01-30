// $Id$
//==============================================================================
//!
//! \file TractionField.h
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Function interfaces for representation of explicit traction fields.
//!
//==============================================================================

#ifndef _TRACTION_FIELD_H
#define _TRACTION_FIELD_H

#include "Function.h"
#include "Tensor.h"

class TensorFunc;
class STensorFunc;


/*!
  \brief Traction field based on a given stress tensor function.

  \details This class defines an explicit traction field function based on
  a specified stress tensor field. The traction field is then evaluated simply
  as the inner-product \b &sigma; &sdot; \b n of the stress tensor value
  \b &sigma; and the outward-directed surface/edge normal vector \b n, which is
  provided as the second argument to the TractionField::evaluate() method.
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

  //! \brief Returns the time-derivative of the function.
  //! \param[in] x Global coordinates of evaluation point
  //! \param[in] n Outward-directed unit normal vector at evaluation point
  virtual Vec3 deriv(const Vec3& x, const Vec3& n) const;

protected:
  //! \brief Evaluates the traction field function at the specified point.
  //! \param[in] x Global coordinates of evaluation point
  //! \param[in] n Outward-directed unit normal vector at evaluation point
  virtual Vec3 evaluate(const Vec3& x, const Vec3& n) const;
};


/*!
  \brief Traction field based on a given scalar pressure function.

  \details This class defines an explicit traction field function based on
  a specified scalar pressure function and a traction direction.
  The traction direction can either be one of the global coordinate axes,
  the surface normal vector (water pressure), or it can also be an
  arbitrary spatial and/or temporal function (the class member #pdfn),
  which must evaluate to a vector of unit length.
*/

class PressureField : public TractionFunc
{
  const RealFunc* pressure; //!< Scalar field to derive the traction field from
  char            pdir;     //!< The global pressure direction (0...3)
  const VecFunc*  pdfn;     //!< Varying global pressure direction

public:
  //! \brief Constructor initializing a constant pressure field.
  //! \param[in] p The constant pressure value
  //! \param[in] dir The global direction the pressure is acting in
  explicit PressureField(Real p, int dir = 0);
  //! \brief Constructor initializing the scalar pressure field function.
  //! \param[in] p The scalar field defining the spatial pressure distribution
  //! \param[in] dir The global direction the pressure is acting in
  explicit PressureField(const RealFunc* p, int dir = 0)
    : pressure(p), pdir(dir), pdfn(nullptr) {}
  //! \brief Constructor initializing the scalar pressure field function.
  //! \param[in] p The scalar field defining the spatial pressure distribution
  //! \param[in] dir The global direction the pressure is acting in
  PressureField(const RealFunc* p, const VecFunc* dir)
    : pressure(p), pdir(-1), pdfn(dir) {}
  //! \brief The destructor frees the scalar and vector field functions.
  virtual ~PressureField() { delete pressure; delete pdfn; }

  //! \brief Returns whether the traction is always normal to the face or not.
  virtual bool isNormalPressure() const { return !pdfn && (pdir<1 || pdir>3); }
  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return pressure ? pressure->isZero() : true; }

  //! \brief Returns the time-derivative of the function.
  //! \param[in] x Global coordinates of evaluation point
  //! \param[in] n Outward-directed unit normal vector at evaluation point
  virtual Vec3 deriv(const Vec3& x, const Vec3& n) const;

protected:
  //! \brief Evaluates the traction field function at the specified point.
  //! \param[in] x Global coordinates of evaluation point
  //! \param[in] n Outward-directed unit normal vector at evaluation point
  virtual Vec3 evaluate(const Vec3& x, const Vec3& n) const;
};


/*!
  \brief Traction field based on a given force resultant and orientation.

  \details This class defines an explicit traction field function based on
  a time-dependent #force magnitude function, a time-dependent vector-valued
  function (#fdir) giving either the force direction as a unit vector or the
  three rotation angles (roll, pitch and yaw) about the local coordinate axes.
  Finally, a spatial #shape function is used to describe the distribution of
  the force resultant over the loaded face.

  The global #force is acting along the Y-axis of the local coordinate
  system, which is defined by the specified origin #X0 and the
  local-to-global transformation #Tlg. This local coordinate system is then
  assigned a time-dependent rotation as given by #fdir. The traction value at
  a given global point \b X and time \a t is therefore computed as
  \f[
    {\bf t}({\bf X},t) = {\bf T}_{lg}\left\{f(t)S({\bf x}(t)){\bf d}(t)\right\}
  \f] where
  \f$f(t)\f$ denotes the value for the #force magnitude function,
  \f$S({\bf x}(t))\f$ is the shape function value at the local point \b x
  and time \a t, and
  \f${\bf d}(t)\f$ is the time-dependent force direction vector of unit length.
  The local, time-dependent coordinates \b x are
  computed from the global coordinate \b X through
  \f[
    {\bf x}(t) = {\bf T}_{rot}(t)\left\{({\bf X}-{\bf X}_0){\bf T}_{lg}\right\}
  \f]
  where \f${\bf T}_{rot}(t)\f$ is the time-dependent rotation matrix
  derived from either the vector-valued function #fdir, or the scalar-valued
  function #angle and the #rotAxis parameter. The Local force direction vector
  \f${\bf d}(t)\f$ is then taken as the second column of this
  transformation matrix (i.e., the Y-axis).
*/

class ForceDirField : public TractionFunc
{
  typedef utl::Function<Real,Vec3> Vec3Func; //!< Convenience type

  Vec3   X0;  //!< Local origin
  Tensor Tlg; //!< Local-to-global transformation

  const ScalarFunc* force; //!< Force resultant magnitude
  const ScalarFunc* angle; //!< Force angle about #rotAxis
  const Vec3Func*   fdir;  //!< Force direction/angles
  const RealFunc*   shape; //!< Shape function for force distribution

  bool dirVec;  //!< If \e true, #fdir is the force direction vector
  char rotAxis; //!< Which local axis the #angle is referring to

public:
  //! \brief The constructor initializes the function pointers.
  //! \param[in] f The scalar function defining the force magnitude
  //! \param[in] d The direction/angles the force is acting in
  //! \param[in] s Shape function defining the force distribution
  //! \param[in] Xaxis Local X-axis direction
  //! \param[in] Zaxis Local Z-axis direction
  //! \param[in] Xorig Origin of the local coordinate system
  //! \param[in] v If \e true, #fdir is the force direction vector
  ForceDirField(const ScalarFunc* f, const Vec3Func* d, const RealFunc* s,
                const Vec3& Xaxis = Vec3(Real(1),Real(0),Real(0)),
                const Vec3& Zaxis = Vec3(Real(0),Real(0),Real(1)),
                const Vec3& Xorig = Vec3(), bool v = false)
    : X0(Xorig), Tlg(Xaxis,Zaxis,false,true), force(f), angle(nullptr), fdir(d),
      shape(s), dirVec(v), rotAxis(0) {}

  //! \brief Alternative constructor for single-axis rotation.
  //! \param[in] f The scalar function defining the force magnitude
  //! \param[in] a The scalar function defining the force angle
  //! \param[in] x Which local axis the angle is referring to
  //! \param[in] s Shape function defining the force distribution
  //! \param[in] Xaxis Local X-axis direction
  //! \param[in] Zaxis Local Z-axis direction
  //! \param[in] Xorig Origin of the local coordinate system
  ForceDirField(const ScalarFunc* f, const ScalarFunc* a, char x,
                const RealFunc* s,
                const Vec3& Xaxis = Vec3(Real(1),Real(0),Real(0)),
                const Vec3& Zaxis = Vec3(Real(0),Real(0),Real(1)),
                const Vec3& Xorig = Vec3())
    : X0(Xorig), Tlg(Xaxis,Zaxis,false,true), force(f), angle(a), fdir(nullptr),
      shape(s), dirVec(false), rotAxis(x) {}

  //! \brief The destructor frees the force functions.
  virtual ~ForceDirField();

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return force ? force->isZero() : true; }

protected:
  //! \brief Evaluates the traction field function at the specified point.
  //! \param[in] x Global coordinates of evaluation point
  virtual Vec3 evaluate(const Vec3& x, const Vec3&) const;
};

#endif
