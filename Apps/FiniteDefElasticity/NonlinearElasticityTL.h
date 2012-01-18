// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityTL.h
//!
//! \date May 25 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_TL_H
#define _NONLINEAR_ELASTICITY_TL_H

#include "Elasticity.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements a Total Lagrangian formulation in matrix form.
  It inherits most of the Elasticity methods, but reimplements the \a kinematics
  method, for calculating the nonlinear variant of the strain-displacement
  matrix, \b B, and the associated Green-Lagrange strain tensor, \b E.
  The \a evalBou method is also reimplemented to account for with-rotated loads.
*/

class NonlinearElasticityTL : public Elasticity
{
public:
  //! \brief The default constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, and axisymmetric 3D formulation is assumed
  NonlinearElasticityTL(unsigned short int n = 3, bool axS = false);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityTL() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of DOFs on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This method is reimplemented in this class to account for
  //! possibly with-rotated traction fields in the Total-Lagrangian setting.
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

protected:
  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[in] F Deformation gradient at current point
  //! \param[out] Bmat The strain-displacement matrix
  //! \param[out] E Green-Lagrange strain tensor at current point
  //!
  //! \details The deformation gradient \b F and the nonlinear
  //! strain-displacement matrix \b B are established.
  //! The B-matrix is formed only when the variable \a formB is true.
  virtual bool kinematics(const Vector& eV,
			  const Vector& N, const Matrix& dNdX, double r,
			  Matrix& Bmat, Tensor& F, SymmTensor& E) const;

protected:
  bool formB; //!< Flag determining whether we need to form the B-matrix
};

#endif
