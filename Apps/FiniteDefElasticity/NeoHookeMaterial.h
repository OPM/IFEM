// $Id$
//==============================================================================
//!
//! \file NeoHookeMaterial.h
//!
//! \date Mar 08 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Neo-Hookean hyperelastic material model.
//!
//==============================================================================

#ifndef _NEO_HOOKE_MATERIAL_H
#define _NEO_HOOKE_MATERIAL_H

#include "LinIsotropic.h"


/*!
  \brief Class representing a Neo-Hookean hyperelastic material model.
  \details This class is a wrapper for the FORTRAN material routines of FENRIS,
  implemented in by the subroutine CONS3D and subroutines called from CONS3D.
  Which material routine to use is governed by the \a mVER parameter.
  Only isotropic material properties are supported by this class.
  \note In 2D, only plane strain is supported by this class.
*/

class NeoHookeMaterial : public LinIsotropic
{
public:
  //! \brief Default constructor.
  NeoHookeMaterial();
  //! \brief Constructor initializing the material parameters.
  NeoHookeMaterial(double E, double v = 0.0, double density = 0.0, int ver = 1);
  //! \brief Empty destructor.
  virtual ~NeoHookeMaterial() {}

  //! \brief Prints out material parameters to the given output stream.
  virtual void print(std::ostream&) const;

  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] F Deformation gradient at current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option;
  //!  0 : Calculate the consitutive matrix only,
  //!  1 : Cauchy stresses and the tangent constitutive matrix,
  //!  2 : 2nd Piola-Kirchhoff stresses and the tangent constitutive matrix,
  //!  3 : Calculate strain energy density only.
  //! \param[in] Fpf Deformation gradient for push-forward transformation
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
			const Vec3& X, const Tensor& F, const SymmTensor& eps,
			char iop = 1, const TimeDomain* = 0,
			const Tensor* Fpf = 0) const;

private:
  int    mTYP;     //!< Material type
  int    mVER;     //!< Material version
  double pmat[13]; //!< Material properties
};

#endif
