// $Id$
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
  NonlinearElasticity(unsigned short int n = 3);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticity() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
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
  //! \param[in] eV Element solution vector
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const Vector& eV,
		       const Matrix& dNdX, const Vec3& X) const;

protected:
  //! \brief Forms tangential tensorial quantities needed by the evalInt method.
  //! \param[out] Ctan Tangential constitutive tensor at current point
  //! \param[out] S 2nd Piola-Kirchhoff stress tensor at current point
  //! \param[in] iGP Global integration point counter
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] F Deformation gradient at current integration point
  //! \param[in] E Green-Lagrange strain tensor at current integration point
  virtual bool formTangent(Matrix& Ctan, SymmTensor& S, size_t iGP,
			   const Vec3& X, const Tensor& F,
			   const SymmTensor& E) const;

  //! \brief Forms the 2nd Piola-Kirchhoff stress tensor.
  //! \param[in] eV Element solution vector
  //! \param[in] dNdX Basis function gradients at current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[out] S 2nd Piola-Kirchhoff stress tensor at current point
  virtual bool formStressTensor(const Vector& eV,
				const Matrix& dNdX, const Vec3& X,
				SymmTensor& S) const;

protected:
  bool fullCmat; //!< If \e true, assume a full (but symmetric) C-matrix
};

#endif
