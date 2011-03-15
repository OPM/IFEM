// $Id$
//==============================================================================
//!
//! \file PlasticMaterial.h
//!
//! \date Mar 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Plasticity material models.
//!
//==============================================================================

#ifndef _PLASTIC_MATERIAL_H
#define _PLASTIC_MATERIAL_H

#include "MaterialBase.h"
#include "Tensor.h"


/*!
  \brief Class containing parameters of a plasticity material model.
*/

class PlasticPrm : public Material
{
public:
  //! \brief Constructor initializing the material parameters.
  PlasticPrm(const RealArray& p, double dens = 0.0) : pMAT(p), rho(dens) {}
  //! \brief Empty destructor.
  virtual ~PlasticPrm() {}

  //! \brief Prints out material parameters to the given output stream.
  virtual void print(std::ostream&) const;

  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return rho; }

  //! \brief Dummy method (should not be invoked).
  virtual bool evaluate(Matrix&, SymmTensor&, double&,
                        const Vec3&, const Tensor&, const SymmTensor&,
                        char, const TimeDomain*) const { return false; }

private:
  RealArray pMAT; //!< Material property parameters
  double    rho;  //!< Mass density

  friend class PlasticMaterial;
};


/*!
  \brief Class representing a plasticity material model.
*/

class PlasticMaterial : public Material
{
public:
  //! \brief Constructor initializing the material parameters.
  PlasticMaterial(const PlasticPrm* prm, unsigned short int n = 0);
  //! \brief Empty destructor.
  virtual ~PlasticMaterial() {}

  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[in] F Deformation gradient at current point
  //! \param[in] prm Nonlinear solution algorithm parameters
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double&,
                        const Vec3&, const Tensor& F, const SymmTensor&,
                        char, const TimeDomain* prm) const;

  //! \brief Returns a reference to the previous deformation gradient.
  Tensor& defGrad() { return Fp; }

  //! \brief Updates the internal history variables after convergence.
  void updateHistoryVars();

private:
  const RealArray& pMAT; //!< Material property parameters

  double HVc[10]; //!< History variables, current configuration
  double HVp[10]; //!< History variables, previous configuration
  Tensor Fp;      //!< Deformation gradient, previous configuration
};

#endif
