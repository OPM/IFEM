// $Id$
//==============================================================================
//!
//! \file LinearMaterial.h
//!
//! \date Mar 08 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General linear elastic material with push-forward transformations.
//!
//==============================================================================

#ifndef _LINEAR_MATERIAL_H
#define _LINEAR_MATERIAL_H

#include "MaterialBase.h"


/*!
  \brief Class representing a isotropic linear elastic material model.
*/

class LinearMaterial : public Material
{
public:
  //! \brief The constructor initializes the material properties pointer.
  LinearMaterial(const Material* mat) : material(mat) {}
  //! \brief The destructor deletes the material properties object.
  virtual ~LinearMaterial() { delete material; }

  //! \brief Prints out material parameters to the given output stream.
  virtual void print(std::ostream& os) const { material->print(os); }

  //! \brief Evaluates the mass density at an integration point.
  virtual double getMassDensity(const Vec3& X) const;

  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] F Deformation gradient at current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option;
  //!  -1 : Calculate the inverse constitutive matrix only,
  //!   0 : Calculate the constitutive matrix only,
  //!   1 : Calculate Cauchy stresses and the constitutive matrix.
  //!   2 : 2nd Piola-Kirchhoff stresses and the tangent constitutive matrix.
  //! \param[in] prm Nonlinear solution algorithm parameters
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
			const Vec3& X, const Tensor& F, const SymmTensor& eps,
			char iop = 1, const TimeDomain* prm = 0) const;

private:
  const Material* material; //!< Pointer to actual material properties object
};

#endif
