// $Id$
//==============================================================================
//!
//! \file PlasticMaterial.h
//!
//! \date Mar 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Elasto-plastic material models.
//!
//==============================================================================

#ifndef _PLASTIC_MATERIAL_H
#define _PLASTIC_MATERIAL_H

#include "MaterialBase.h"
#include "Tensor.h"


/*!
  \brief Class containing parameters of an elasto-plastic material model.
*/

class PlasticPrm : public Material
{
public:
  //! \brief Constructor initializing the material parameters.
  PlasticPrm(const RealArray&);
  //! \brief Empty destructor.
  virtual ~PlasticPrm() {}

  //! \brief Prints out material parameters to the given output stream.
  virtual void print(std::ostream&) const;

  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return pMAT[3]; }

  //! \brief Dummy method (should not be invoked).
  virtual bool evaluate(Matrix&, SymmTensor&, double&,
                        const Vec3&, const Tensor&, const SymmTensor&,
                        char, const TimeDomain*) const { return false; }

private:
  RealArray pMAT; //!< Material property parameters

  friend class PlasticMaterial;
};


/*!
  \brief Class representing an elasto-plastic material point.
*/

class PlasticMaterial : public Material
{
public:
  //! \brief Constructor initializing the material parameters.
  //! \param[in] prm Pointer to actual material parameters object.
  //! \param[in] n Number of space dimensions
  PlasticMaterial(const PlasticPrm* prm, unsigned short int n);
  //! \brief Empty destructor.
  virtual ~PlasticMaterial() {}

  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return pMAT[3]; }

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
