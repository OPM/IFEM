// $Id$
//==============================================================================
//!
//! \file Integrand.h
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for classes representing FEM integrands.
//!
//==============================================================================

#ifndef _INTEGRAND_H
#define _INTEGRAND_H

#include "MatVec.h"

struct TimeDomain;
class LocalIntegral;
class FiniteElement;
class MxFiniteElement;
class Vec3;


/*!
  \brief Abstract base class representing a system level integrated quantity.
  \details This class defines the interface between the finite element (FE)
  assembly drivers of the ASM-hierarchy and the problem-dependent classes
  containing all physical properties for the problem to be solved.

  The interface consists of methods for evaluating the integrand at interior
  integration points (\a evalInt), and at boundary integration points
  (\a evalBou). The latter are used for Neumann boundary conditions, typically.
  The integrand evaluation methods have access to the FE basis function values
  and derivatives through the FiniteElement argument. There are also a set
  of methods dedicated for mixed field interpolation problems, which take
  a MxFiniteElement object as argument instead.
*/

class Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  Integrand() {}

public:
  //! \brief Empty destructor.
  virtual ~Integrand() {}


  // Global initialization interface
  // ===============================

  //! \brief Initializes the integrand for a new integration loop.
  //! \details This method is invoked once before starting the numerical
  //! integration over the entire spatial domain.
  virtual void initIntegration(const TimeDomain&) = 0;


  // Element-level initialization interface
  // ======================================

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  //! The default implementation forwards to an overloaded method not taking
  //! \a X0 and \a nPt as arguments.
  //! Reimplement this method for problems requiring the element center and/or
  //! the number of integration points during/before the integrand evaluations.
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3& X0, size_t nPt)
  {
    return this->initElement(MNPC);
  }
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  //! Reimplement this method for problems \e not requiring the
  //! the element center nor the number of integration points before the
  //! integration loop is started.
  virtual bool initElement(const std::vector<int>& MNPC) = 0;
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  virtual bool initElement(const std::vector<int>& MNPC1,
			   const std::vector<int>& MNPC2, size_t n1) = 0;

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC) = 0;
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  virtual bool initElementBou(const std::vector<int>& MNPC1,
			      const std::vector<int>& MNPC2, size_t n1) = 0;


  // Integrand evaluation interface
  // ==============================

  //! \brief Defines which FE quantities are needed by the integrand.
  //! \return 1: The basis functions and their gradients are needed
  //! \return 2: As 1, but in addition the second derivatives of the basis
  //! functions and the characteristic element size is needed
  //! \return 3: As 1, but in addition the volume-averaged basis functions
  //! are needed
  //! \return 4: As 1, but in addition the element center coordinates are needed
  virtual int getIntegrandType() const { return 1; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const TimeDomain& time, const Vec3& X) const
  {
    return this->evalInt(elmInt,fe,X);
  }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Mixed finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used for mixed formulations only.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalIntMx(LocalIntegral*& elmInt, const MxFiniteElement& fe,
			 const TimeDomain& time, const Vec3& X) const
  {
    return this->evalIntMx(elmInt,fe,X);
  }

  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  //! It can also be used to implement multiple integration point loops within
  //! the same element, provided the necessary integration point values are
  //! stored internally in the object during the first integration loop.
  virtual bool finalizeElement(LocalIntegral*&, const TimeDomain&)
  {
    return true;
  }

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalBou(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const TimeDomain& time,
		       const Vec3& X, const Vec3& normal) const
  {
    return this->evalBou(elmInt,fe,X,normal);
  }

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Mixed finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This interface is used for mixed formulations.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalBouMx(LocalIntegral*& elmInt, const MxFiniteElement& fe,
			 const TimeDomain& time,
			 const Vec3& X, const Vec3& normal) const
  {
    return this->evalBouMx(elmInt,fe,X,normal);
  }

protected:
  //! \brief Evaluates the integrand at interior points for stationary problems.
  virtual bool evalInt(LocalIntegral*&, const FiniteElement& fe,
		       const Vec3&) const { return false; }
  //! \brief Evaluates the integrand at interior points for stationary problems.
  virtual bool evalIntMx(LocalIntegral*&, const MxFiniteElement& fe,
			 const Vec3&) const { return false; }

  //! \brief Evaluates the integrand at boundary points for stationary problems.
  virtual bool evalBou(LocalIntegral*&, const FiniteElement&,
		       const Vec3&, const Vec3&) const { return false; }
  //! \brief Evaluates the integrand at boundary points for stationary problems.
  virtual bool evalBouMx(LocalIntegral*&, const MxFiniteElement&,
			 const Vec3&, const Vec3&) const { return false; }

public:

  // Solution field evaluation interface
  // ===================================

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] N1 Basis function values at current point, field 1
  //! \param[in] N2 Basis function values at current point, field 2
  //! \param[in] dN1dX Basis function gradients at current point, field 1
  //! \param[in] dN2dX Basis function gradients at current point, field 2
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  virtual bool evalSol(Vector& s,
		       const Vector& N1, const Vector& N2,
		       const Matrix& dN1dX, const Matrix& dN2dX, const Vec3& X,
		       const std::vector<int>& MNPC1,
		       const std::vector<int>& MNPC2) const;
};

#endif
