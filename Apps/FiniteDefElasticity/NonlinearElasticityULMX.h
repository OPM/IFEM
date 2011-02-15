// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityULMX.h
//!
//! \date Dec 14 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity mixed problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_UL_MX_H
#define _NONLINEAR_ELASTICITY_UL_MX_H

#include "NonlinearElasticityUL.h"
#include "Tensor.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements a mixed Updated Lagrangian formulation with
  internal pressure and volumetric change modes.
*/

class NonlinearElasticityULMX : public NonlinearElasticityUL
{
  //! \brief A struct with integration point data needed by \a finalizeElement.
  struct ItgPtData
  {
    Tensor F;     //!< Deformation gradient at current step/iteration
    Tensor Fp;    //!< Deformation gradient at previous (converged) step
    Matrix dNdx;  //!< Basis function gradients at current configuration
    Vector Phi;   //!< Internal modes for the pressure/volumetric-change fields
    Vec4   X;     //!< Cartesian coordinates of current integration point
    double detJW; //!< Jacobian determinant times integration point weight
    //! \brief Default constructor.
    ItgPtData() : F(3), Fp(3) { detJW = 0.0; }
  };

public:
  //! \brief The default constructor invokes the parent class constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] mver Material version parameter
  //! \param[in] pp Polynomial order of the pressure/volumetric-change field
  NonlinearElasticityULMX(unsigned short int n = 3, int mver = 0, int pp = 1);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityULMX() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3& Xc, size_t nPt);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This method mainly stores the integration point quantities
  //! depending on the basis function values in the internal member \a myData.
  //! The actual numerical integration of the tangent stiffness and internal
  //! forces is then performed by the \a finalizeElement method.
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const;

  //! \brief Finalizes the element matrices after the numerical integration.
  //! \details This method is used to implement the multiple integration point
  //! loops within an element, required by the mixed formulation with internal
  //! modes. All data needed during the integration has been stored internally
  //! in the \a myData member by the \a evalInt method.
  //! \param elmInt The local integral object to receive the contributions
  virtual bool finalizeElement(LocalIntegral*&);

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return 4; }

private:
  int  p;  //!< Polynomial order of the internal pressure field
  Vec3 X0; //!< Cartesian coordinates of the element center

  Matrix* Hh; //!< Pointer to element Hh-matrix associated with pressure modes

  mutable size_t                 iP;     //!< Local integration point counter
  mutable std::vector<ItgPtData> myData; //!< Local integration point data

  friend class ElasticityNormULMX;
};


/*!
  \brief Class representing the integrand of the mixed elasticity energy norm.
*/

class ElasticityNormULMX : public ElasticityNormUL
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  ElasticityNormULMX(NonlinearElasticityULMX& p) : ElasticityNormUL(p) {}
  //! \brief Empty destructor.
  virtual ~ElasticityNormULMX() {}

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3& Xc, size_t nPt);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This method just forwards to the \a evalInt method of the
  //! corresponding NonlinearElasticityULMX object, for which this class is
  //! declared as friend such that it can access the data members. The actual
  //! norm integration us then performed by the \a finalizeElement method.
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const;

  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to implement the multiple integration point
  //! loops within an element required by the mixed formulation. It performs
  //! a subset of the tasks done by NonlinearElasticityULMX::finalizeElement.
  //! \param elmInt The local integral object to receive the contributions
  virtual bool finalizeElement(LocalIntegral*& elmInt);

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return 4; }
};

#endif
