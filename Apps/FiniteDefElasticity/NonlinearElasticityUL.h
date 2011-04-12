// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityUL.h
//!
//! \date Sep 21 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_UL_H
#define _NONLINEAR_ELASTICITY_UL_H

#include "Elasticity.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements an Updated Lagrangian formulation. It inherits
  most of the Elasticity methods, but reimplements the \a kinematics method
  for calculating the deformation gradient and the associated Green-Lagrange
  strain tensor. The \a evalInt and \a evalBou methods are also reimplemented
  to account for the updated geometry.
*/

class NonlinearElasticityUL : public Elasticity
{
public:
  //! \brief The default constructor invokes the parent class constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] lop Load option (0=on initial length, 1=on updated length)
  NonlinearElasticityUL(unsigned short int n = 3, char lop = 0);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityUL() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes the integrand for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  virtual void initIntegration(const TimeDomain& prm);
  //! \brief Initializes the integrand for a new result point loop.
  //! \param[in] lambda Load parameter
  virtual void initResultPoints(double lambda);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This method is reimplemented in this class to account for
  //! possibly with-rotated traction fields (non-conservative loads).
  //! For uni-directional (conservative) loads, it is similar to the
  //! \a LinearElasticity::evalBou method.
  virtual bool evalBou(LocalIntegral*& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const;

  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[out] F Deformation gradient at current point
  //! \param[out] E Green-Lagrange strain tensor at current point
  virtual bool kinematics(const Matrix& dNdX, Tensor& F, SymmTensor& E) const;

protected:
  //! \brief Calculates the deformation gradient at current point.
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[out] F Deformation gradient at current point
  bool formDefGradient(const Matrix& dNdX, Tensor& F) const;

protected:
  mutable Matrix dNdx; //!< Basis function gradients in current configuration
  mutable Matrix CB;   //!< Result of the matrix-matrix product C*B

private:
  char loadOp; //!< Load option
  double plam; //!< Load parameter of the previous result evaluation

  friend class ElasticityNormUL;
};


/*!
  \brief Class representing the integrand of the elasticity energy norm.
  \details This class reimplements the \a evalInt method to use the strain
  energy density value returned by the nonlinear constitutive model.
  It also reimplements the \a evalBou method with a path integral of
  the external energy due to boundary tractions.
*/

class ElasticityNormUL : public ElasticityNorm
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  ElasticityNormUL(NonlinearElasticityUL& p) : ElasticityNorm(p), iP(0) {}
  //! \brief Empty destructor.
  virtual ~ElasticityNormUL() {}

  //! \brief Initializes the integrand for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  virtual void initIntegration(const TimeDomain& prm);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

private:
  // Data for path-integral of the external energy due to boundary tractions
  mutable size_t            iP; //!< Global integration point counter
  mutable RealArray         Ux; //!< External energy densities
  mutable std::vector<Vec3> up; //!< Previous displacements
  mutable std::vector<Vec3> tp; //!< Previous tractions
};

#endif
