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


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements a mixed Updated Lagrangian formulation with
  internal pressure and volumetric change modes.
*/

class NonlinearElasticityULMX : public NonlinearElasticityUL
{
public:
  //! \brief The default constructor invokes the parent class constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, and axisymmetric 3D formulation is assumed
  //! \param[in] pp Polynomial order of the pressure/volumetric-change field
  NonlinearElasticityULMX(unsigned short int n = 3,
			  bool axS = false, int pp = 1);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityULMX() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
					  bool neumann) const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3& Xc, size_t nPt,
			   LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This method mainly stores the integration point quantities
  //! depending on the basis function values in the internal member \a myData.
  //! The actual numerical integration of the tangent stiffness and internal
  //! forces is then performed by the \a finalizeElement method.
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  //! \brief Finalizes the element matrices after the numerical integration.
  //! \details This method is used to implement the multiple integration point
  //! loops within an element, required by the mixed formulation with internal
  //! modes. All data needed during the integration has been stored internally
  //! in the \a elmInt object by the \a evalInt method.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] iG Global index of the first integration point in the element
  virtual bool finalizeElement(LocalIntegral& elmInt, const TimeDomain& prm,
			       size_t iG);

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return ELEMENT_CENTER; }
  //! \brief Returns the number of reduced-order integration points.
  //! \return -1 which is used to flag that the integrand needs to know the
  //! number of integration points within each element.
  //! \sa The argument \a nPt of initElement.
  virtual int getReducedIntegration() const { return -1; }

private:
  int p; //!< Polynomial order of the internal pressure field

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

  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl The element number
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
					  bool neumann) const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3& Xc, size_t nPt,
			   LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt The local integral object for current element
  virtual bool initElementBou(const std::vector<int>& MNPC,
			      LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This method just forwards to the \a evalInt method of the
  //! corresponding NonlinearElasticityULMX object, for which this class is
  //! declared as friend such that it can access the data members. The actual
  //! norm integration us then performed by the \a finalizeElement method.
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to implement the multiple integration point
  //! loops within an element required by the mixed formulation. It performs
  //! a subset of the tasks done by NonlinearElasticityULMX::finalizeElement.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] iG Global index of the first integration point in the element
  virtual bool finalizeElement(LocalIntegral& elmInt, const TimeDomain& prm,
			       size_t iG);

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return ELEMENT_CENTER; }
  //! \brief Returns the number of reduced-order integration points.
  //! \return -1 which is used to flag that the integrand needs to know the
  //! number of integration points within each element.
  //! \sa The argument \a nPt of initElement.
  virtual int getReducedIntegration() const { return -1; }
};

#endif
