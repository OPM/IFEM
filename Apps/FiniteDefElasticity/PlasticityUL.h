// $Id$
//==============================================================================
//!
//! \file PlasticityUL.h
//!
//! \date Sep 21 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for finite deformation plasticity problems.
//!
//==============================================================================

#ifndef _PLASTICITY_UL_H
#define _PLASTICITY_UL_H

#include "NonlinearElasticityUL.h"

class PlasticMaterial;
struct TimeDomain;


/*!
  \brief Class representing the integrand of the nonlinear plasticity problem.
  \details This class implements an Updated Lagrangian formulation. It inherits
  most of the NonlinearElasticityUL methods, but reimplements \a evalInt using
  an internal buffer to store history variables in every integration point
  between the iterations.
*/

class PlasticityUL : public NonlinearElasticityUL
{
public:
  //! \brief The default constructor invokes the parent class constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] lop Load option (0=on initial length, 1=on updated length)
  PlasticityUL(unsigned short int n = 3, char lop = 0);
  //! \brief The destructor frees the dynamically allocated integration buffers.
  virtual ~PlasticityUL();

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Initializes the integrand for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration over the entire spatial domain, and is used to reset the
  //! internal integration point counter, and to update the history variables.
  virtual void initIntegration(const TimeDomain& prm);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt,
		       const TimeDomain& prm, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const;

  //! \brief Returns a null pointer for solution norm evaluation.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const { return 0; }

private:
  mutable size_t                        iP;   //!< Integration point counter
  mutable std::vector<PlasticMaterial*> pBuf; //!< Integration point buffers
};

#endif
