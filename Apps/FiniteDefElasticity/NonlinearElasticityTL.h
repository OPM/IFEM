// $Id: NonlinearElasticityTL.h,v 1.2 2011-02-08 09:06:02 kmo Exp $
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
  //! \param[in] ps If \e true, assume plane stress in 2D
  NonlinearElasticityTL(unsigned short int n = 3, bool ps = true);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityTL() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This method is reimplemented in this class to account for
  //! possibly with-rotated traction fields in the Total-Lagrangian setting.
  virtual bool evalBou(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X, const Vec3& normal) const;

protected:
  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[out] E Green-Lagrange strain tensor at current point
  //!
  //! \details The deformation gradient \b F and the nonlinear
  //! strain-displacement matrix \b B are established. The latter matrix
  //! is stored in the mutable class member \a Bmat of the parent class.
  //! The B-matrix is formed only when the variable \a formB is true.
  virtual bool kinematics(const Matrix& dNdX, SymmTensor& E) const;

protected:
  bool       formB; //!< Flag determining whether we need to form the B-matrix
  mutable Matrix F; //!< Deformation gradient
};

#endif
