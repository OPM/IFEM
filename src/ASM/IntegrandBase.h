// $Id: IntegrandBase.h,v 1.20 2011-02-08 12:10:25 rho Exp $
//==============================================================================
//!
//! \file IntegrandBase.h
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for classes representing FEM integrands.
//!
//==============================================================================

#ifndef _INTEGRAND_BASE_H
#define _INTEGRAND_BASE_H

#include "SIMenums.h"
#include "Function.h"
#include "MatVec.h"

struct TimeDomain;
class LocalIntegral;
class NormBase;
class AnaSol;
class Vec3;
class VTF;


/*!
  \brief Abstract base class representing a system level integrated quantity.
  \details This class defines the interface between the finite element
  assembly drivers of the ASM-hierarchy and the problem-dependent classes
  containing all physical properties for the problem to be solved.

  The interface has several overloaded versions of the \a evalInt method, for
  evaluation of the coefficient matrix contributions at an element-interior
  integration point. Only one of these methods needs to be implemented, and
  which one to use is governed by the \a getIntegrandType method.
  There is also one method intended for mixed field interpolation problems.
*/

class Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  Integrand() {}

public:
  //! \brief Empty destructor.
  virtual ~Integrand() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream&) const {}

  //! \brief Defines the solution mode before the element assembly is started.
  virtual void setMode(SIM::SolutionMode) {}

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
  //! \a X0  and \a nPt as arguments.
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
  //! \details Reimplement this method for problems not requiring the
  //! the element center nor the number of integration points before the
  //! integration loop is started.
  virtual bool initElement(const std::vector<int>& MNPC)
  {
    std::cerr <<" *** Integrand::initElement not implemented."<< std::endl;
    return false;
  }
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  virtual bool initElement(const std::vector<int>& MNPC1,
			   const std::vector<int>& MNPC2, size_t n1)
  {
    std::cerr <<" *** Integrand::initElement not implemented."<< std::endl;
    return false;
  }

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC)
  {
    std::cerr <<" *** Integrand::initElementBou not implemented."<< std::endl;
    return false;
  }
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  virtual bool initElementBou(const std::vector<int>& MNPC1,
			      const std::vector<int>& MNPC2, size_t n1)
  {
    std::cerr <<" *** Integrand::initElementBou not implemented."<< std::endl;
    return false;
  }

  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop is finished, and before the resulting element quantities
  //! are assembled into their system level equivalents.
  //! It can also be used to implement multiple integration point loops within
  //! the same element, provided the necessary integration point values are
  //! stored internally in the object during the first integration loop.
  virtual bool finalizeElement(LocalIntegral*&) { return true; }
  

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used when \a getIntegrandType returns 1.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const
  {
    return this->evalInt(elmInt,detJW,N,dNdX,X);
  }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] d2NdX2 Basis function second derivatives
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] h Characteristic element length
  //!
  //! \details This interface is used when \a getIntegrandType returns 2.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Matrix3D& d2NdX2, const Vec3& X,
		       double h = 0.0) const
  {
    return this->evalInt(elmInt,detJW,N,dNdX,d2NdX2,X,h);
  }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] Navg Volume-averaged basis function values over the element
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used when \a getIntegrandType returns 3.
  //! Use when the integrand requires second-derivatives of the basis functions.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vector& Navg, const Vec3& X) const
  {
    return this->evalInt(elmInt,detJW,N,dNdX,Navg,X);
  }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N1 Basis function values, field 1
  //! \param[in] N2 Basis function values, field 2
  //! \param[in] dN1dX Basis function gradients, field 1
  //! \param[in] dN2dX Basis function gradients, field 2
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used for mixed formulations only.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N1, const Vector& N2,
		       const Matrix& dN1dX, const Matrix& dN2dX,
		       const Vec3& X) const
  {
    return this->evalInt(elmInt,detJW,N1,N2,dN1dX,dN2dX,X);
  }

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalBou(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const Vec3& normal) const
  {
    return this->evalBou(elmInt,detJW,N,dNdX,X,normal);
  }

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N1 Basis function values, field 1
  //! \param[in] N2 Basis function values, field 2
  //! \param[in] dN1dX Basis function gradients, field 1
  //! \param[in] dN2dX Basis function gradients, field 2
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This interface is used for mixed formulations.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalBou(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N1, const Vector& N2,
		       const Matrix& dN1dX, const Matrix& dN2dX,
		       const Vec3& X, const Vec3& normal) const
  {
    return this->evalBou(elmInt,detJW,N1,N2,dN1dX,dN2dX,X,normal);
  }

protected:
  //! \brief Evaluates the integrand at interior points for stationary problems.
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const { return false; }
  //! \brief Evaluates the integrand at interior points for stationary problems.
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Matrix3D& d2NdX2, const Vec3& X,
		       double h = 0.0) const { return false; }
  //! \brief Evaluates the integrand at interior points for stationary problems.
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX, const Vector& Navg,
		       const Vec3& X) const { return false; }
  //! \brief Evaluates the integrand at interior points for stationary problems.
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N1, const Vector& N2,
		       const Matrix& dN1dX, const Matrix& dN2dX,
		       const Vec3& X) const { return false; }

  //! \brief Evaluates the integrand at boundary points for stationary problems.
  virtual bool evalBou(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX, const Vec3& X,
		       const Vec3& normal) const { return false; }
  //! \brief Evaluates the integrand at boundary points for stationary problems.
  virtual bool evalBou(LocalIntegral*& elmInt, double detJW,
		       const Vector& N1, const Vector& N2,
		       const Matrix& dN1dX, const Matrix& dN2dX, const Vec3& X,
		       const Vec3& normal) const { return false; }

public:
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const
  {
    std::cerr <<" *** Integrand::evalSol not implemented"<< std::endl;
    return false;
  }

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSecSol(Vector& s,
                          const Vector& N, const Matrix& dNdX,
                          const Vec3& X, const std::vector<int>& MNPC) const
  {
    return false;
  }

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
		       const std::vector<int>& MNPC2) const
  {
    return this->evalSol(s,N1,dN1dX,X,MNPC1);
  }

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const VecFunc& asol, const Vec3& X) const
  {
    std::cerr <<" *** Integrand::evalSol (exact) not implemented"<< std::endl;
    return false;
  }

  //! \brief Evaluates the analytical secondary solution at an integration point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSecSol(Vector& s, const TensorFunc& asol, const Vec3& X) const
  {
    return false;
  }

  //! \brief Evaluates the analytical scalar solution at an integration point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (vector field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSolScal(real& s, const RealFunc& asol, const Vec3& X) const
  {
    std::cerr <<" *** Integrand::evalSol (exact) not implemented"<< std::endl;
    return false;
  }

  //! \brief Evaluates the analytical secondary scalar solution at an integration point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (vector field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSecSolScal(Vector& s, const VecFunc& asol, const Vec3& X) const
  {
    return false;
  }

  //! \brief Writes surface tractions/fluxes for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param nBlck Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& nBlck) const { return true; }

  //! \brief Returns whether there are any traction/flux values to write to VTF.
  virtual bool hasTractionValues() const { return false; }

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return 1; }

  //! \brief Returns the number of secondary solution field components.
  virtual size_t getNoFields() const { return 0; }

  //! \brief Returns the number of secondary solution field components.
  virtual size_t getNoSecFields() const { return 0; }

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getFieldLabel(size_t i, const char* prefix = 0) const
  {
    return 0;
  }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \param[in] asol Pointer to the analytical solution field (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const { return 0; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \param[in] asol Pointer to the analytical solution field (optional)
  //!
  //! \details This version is used for scalar problems only.
  virtual NormBase* getNormIntegrandScal(AnaSol* asol = 0) const { return 0; }

  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions() const { return primsol.size(); }

  //! \brief Accesses the primary solution vector of current patch.
  Vector& getSolution(size_t n = 0) { return primsol[n]; }

protected:
  Vectors primsol; //!< Primary solution vectors for this patch
};


/*!
  \brief Abstract base class representing a system level norm quantity.
*/

class NormBase : public Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  NormBase() {}

public:
  //! \brief Empty destructor.
  virtual ~NormBase() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return false; }

  //! \brief Returns a 1-based index of the external energy norm.
  virtual size_t indExt() const { return 0; }
};

#endif
