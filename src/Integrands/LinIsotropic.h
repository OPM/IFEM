// $Id$
//==============================================================================
//!
//! \file LinIsotropic.h
//!
//! \date Mar 01 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Isotropic linear elastic material model.
//!
//==============================================================================

#ifndef _LIN_ISOTROPIC_H
#define _LIN_ISOTROPIC_H

#include "MaterialBase.h"


/*!
  \brief Class representing an isotropic linear elastic material model.
*/

class LinIsotropic : public Material
{
public:
  //! \brief Default constructor.
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  LinIsotropic(bool ps = false, bool ax = false);
  //! \brief Constructor initializing the material parameters.
  //! \param[in] E Young's modulus
  //! \param[in] v Poisson's ratio
  //! \param[in] density Mass density
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  LinIsotropic(double E, double v = 0.0, double density = 0.0,
	       bool ps = false, bool ax = false)
    : Emod(E), nu(v), rho(density), planeStress(ps), axiSymmetry(ax) {}
  //! \brief Empty destructor.
  virtual ~LinIsotropic() {}

  //! \brief Returns \e false if plane stress in 2D.
  virtual bool isPlaneStrain() const { return !planeStress; }

  //! \brief Prints out material parameters to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return rho; }

  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option;
  //!  -1 : Calculate the inverse constitutive matrix only,
  //!   0 : Calculate the constitutive matrix only,
  //!   1 : Calculate Cauchy stresses and the constitutive matrix,
  //!   3 : Calculate the strain energy density only.
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
			const Vec3&, const Tensor&, const SymmTensor& eps,
			char iop = 1, const TimeDomain* = 0,
			const Tensor* = 0) const;

protected:
  // Material properties (constant)
  double Emod;        //!< Young's modulus
  double nu;          //!< Poisson's ratio
  double rho;         //!< Mass density
  bool   planeStress; //!< Plane stress/strain option for 2D problems
  bool   axiSymmetry; //!< Axi-symmetric option
};

#endif
