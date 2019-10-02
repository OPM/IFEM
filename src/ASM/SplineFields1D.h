// $Id$
//==============================================================================
//!
//! \file SplineFields1D.h
//!
//! \date Nov 23 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 1D.
//!
//==============================================================================

#ifndef _SPLINE_FIELDS_1D_H
#define _SPLINE_FIELDS_1D_H

#include "Fields.h"

namespace Go {
  class SplineCurve;
}


/*!
  \brief Class for spline-based finite element vector fields in 1D.

  \details This class implements the methods required to evaluate a 1D
  spline vector field at a given point in parametrical or physical coordinates.
*/

class SplineFields1D : public Fields
{
public:
  //! \brief Construct directly from curve.
  //! \param[in] crv The spline curve to use
  //! \param[in] v Array of control point field values
  //! \param[in] ncmp Number of field components
  //! \param[in] name Name of spline field
  SplineFields1D(const Go::SplineCurve* crv, const RealArray& v, int ncmp,
                 const char* name = nullptr);

  //! \brief Empty destructor.
  virtual ~SplineFields1D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] vals Values at local point in given element
  virtual bool valueFE(const ItgPoint& x, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] grad Gradient at local point in given element
  virtual bool gradFE(const ItgPoint& x, Matrix& grad) const;

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] H Hessian at local point in given element
  virtual bool hessianFE(const ItgPoint& x, Matrix3D& H) const;

protected:
  const Go::SplineCurve* curv; //!< Spline geometry description

private:
  unsigned char nsd; //!< Number of space dimensions
};

#endif
