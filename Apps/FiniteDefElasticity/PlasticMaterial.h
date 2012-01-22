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

class PlasticMaterial;


/*!
  \brief Class representing an history-dependent elasto-plastic material model.
  \details This class is a wrapper for the plasticity routines of FENRIS,
  implemented in by the FORTRAN subroutine PLAS3D.

  The plasticity models have history variables at the integration points that
  need to be "remembered" from one iteration/increment to the next. This is
  maintained inside this class as a vector of integration point objects, and
  a global counter keeping track of which integration point we are calculating.
  This counter is initialized in the beginning of each iteration, and then
  incremented by each invokation of the \a evaluate method. Therefore, it is
  paramount that this method is invoked only once per integration point per
  iteration, and in the same order in every iteration. Otherwise it will
  output incorrect results.

  A separate "integration point" vector is dedicated for results points.
  This is needed because the results points (for visualization, etc.) are not
  the same as the integration points used in the tangent evaluation.
*/

class PlasticMaterial : public Material
{
  //! \brief Class representing an elasto-plastic material point.
  class PlasticPoint
  {
  public:
    //! \brief Constructor initializing the material parameters.
    //! \param[in] prm Pointer to actual material model object.
    //! \param[in] n Number of space dimensions
    PlasticPoint(const PlasticMaterial* prm, unsigned short int n);

    //! \brief Updates the internal history variables after convergence.
    bool updateHistoryVars();

    //! \brief Evaluates the constitutive relation at this point.
    //! \param[out] C Constitutive matrix at current point
    //! \param[out] sigma Stress tensor at current point
    //! \param[in] F Deformation gradient at current point
    //! \param[in] prm Nonlinear solution algorithm parameters
    bool evaluate(Matrix& C, SymmTensor& sigma,
                  const Tensor& F, const TimeDomain& prm) const;

    //! \brief Updates the path integral of the strain energy density.
    //! \param[in] S Stress tensor at current configuration
    //! \param[in] E Strain tensor at current configuration
    //! \return Updated strain energy density
    double energyIntegral(const SymmTensor& S, const SymmTensor& E);

    //! \brief Returns a history variable.
    double getVariable(size_t i) const { return HVc[i]; }

  private:
    const RealArray& pMAT; //!< Material property parameters

    double HVc[10]; //!< History variables, current configuration
    double HVp[10]; //!< History variables, previous configuration
    bool   updated; //!< Flag indicating whether history variables are updated

    // Data for path integral of strain energy
    SymmTensor Ep; //!< Strain tensor, previous configuration
    SymmTensor Sp; //!< Stress tensor, previous configuration
    double     Up; //!< Strain energy density

  public:
    Tensor Fp; //!< Deformation gradient, previous configuration
  };

public:
  //! \brief Constructor initializing the material parameters.
  PlasticMaterial(const RealArray&);
  //! \brief The destructor frees the dynamically allocated data.
  virtual ~PlasticMaterial();

  //! \brief Prints out material parameters to the given output stream.
  virtual void print(std::ostream&) const;

  //! \brief Initializes the material with the number of integration points.
  virtual void initIntegration(size_t nGP);
  //! \brief Initializes the material model for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration over the entire spatial domain, and is used to reset the
  //! global integration point counter, and to update the history variables.
  virtual void initIntegration(const TimeDomain& prm);
  //! \brief Initializes the material model for a new result point loop.
  virtual void initResultPoints();

  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return pMAT[3]; }

  //! \brief Evaluates the constitutive relation at current integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] iP1 Global index for current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] F Deformation gradient at current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option;
  //!  0 : Calculate the consitutive matrix only,
  //!  1 : Cauchy stresses and the tangent constitutive matrix,
  //!  2 : 2nd Piola-Kirchhoff stresses and the tangent constitutive matrix,
  //!  3 : Calculate strain energy density only.
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] Fpf Deformation gradient for push-forward transformation
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U, size_t iP1,
                        const Vec3& X, const Tensor& F, const SymmTensor& eps,
                        char iop, const TimeDomain* prm,
                        const Tensor* Fpf = 0) const;

  //! \brief Returns number of internal result variables of the material model.
  virtual int getNoIntVariables() const;
  //! \brief Returns an internal variable associated with the material model.
  //! \param[in] index Index of the internal variable
  //! \param[out] label Name of the internal variable (for result presentation)
  virtual double getInternalVariable(int index, char* label = 0) const;

private:
  friend class PlasticPoint;

  RealArray     pMAT; //!< Material property parameters
  mutable size_t iP2; //!< Global result point counter

  bool iAmIntegrating; //!< Flag indicating integration or result evaluation

  mutable std::vector<PlasticPoint*> itgPoints; //!< Integration point data
  mutable std::vector<PlasticPoint*> resPoints; //!< Result point data
};

#endif
