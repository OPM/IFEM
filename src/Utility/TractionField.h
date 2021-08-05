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

class TensorFunc;
class STensorFunc;


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
  //! \brief Evaluates the traction at point \a x and with surface normal \a n.
  virtual Vec3 evaluate(const Vec3& x, const Vec3& n) const;
};


/*!
  \brief Traction field based on a given pressure function.
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

protected:
  //! \brief Evaluates the traction at point \a x and with surface normal \a n.
  virtual Vec3 evaluate(const Vec3& x, const Vec3& n) const;
};

#endif
