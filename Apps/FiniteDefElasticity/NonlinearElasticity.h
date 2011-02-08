// $Id: NonlinearElasticity.h,v 1.2 2011-02-08 09:06:02 kmo Exp $
//==============================================================================
//!
//! \file NonlinearElasticity.h
//!
//! \date May 25 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_H
#define _NONLINEAR_ELASTICITY_H

#include "NonlinearElasticityTL.h"
#include "Tensor.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements a Total Lagrangian formulation, tensorial form.
  It reimplements most of the Elasticity methods, except for the methods
  \a kinematics and \a evalBou, which are inherited from NonlinearElasticityTL.

  \note This class is obsolete, as it gives exactly the same results as the
  NonlinearElasticityTL class, but is less efficient. It is retained mostly for
  educational purposes and historical reasons, and such that the performance of
  the two ways of implementing the total Lagrangian formulation can be compared
  on a varity of problems.
*/

class NonlinearElasticity : public NonlinearElasticityTL
{
public:
  //! \brief The default constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] ps If \e true, assume plane stress in 2D
  NonlinearElasticity(unsigned short int n = 3, bool ps = true);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticity() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const Matrix& dNdX, const Vec3& X) const;

protected:
  //! \brief Forms tangential tensorial quantities needed by the evalInt method.
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] Ctan Tangential constitutive tensor at current point
  //! \param[out] S 2nd Piola-Kirchhoff stress tensor at current point
  virtual bool formTangent(Matrix& Ctan, SymmTensor& S, const Vec3& X) const;

  //! \brief Forms the 2nd Piola-Kirchhoff stress tensor.
  //! \param[in] dNdX Basis function gradients at current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] S 2nd Piola-Kirchhoff stress tensor at current point
  virtual bool formStressTensor(const Matrix& dNdX, const Vec3& X,
				SymmTensor& S) const;

protected:
  bool        fullCmat; //!< If \e true, assume a full (but symmetric) C-matrix
  mutable SymmTensor E; //!< Green-Lagrange strain tensor
};

#endif
