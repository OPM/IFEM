// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityFbar.h
//!
//! \date Aug 18 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity with Fbar method.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_FBAR_H
#define _NONLINEAR_ELASTICITY_FBAR_H

#include "NonlinearElasticityUL.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements an Updated Lagrangian F-bar formulation.
*/

class NonlinearElasticityFbar : public NonlinearElasticityUL
{
  //! \brief A struct with volumetric sampling point data.
  struct VolPtData
  {
    double J;    //!< Determinant of current deformation gradient
    Matrix dNdx; //!< Basis function gradients at current configuration
    //! \brief Default constructor.
    VolPtData() { J = 1.0; }
  };

public:
  //! \brief The default constructor invokes the parent class constructor.
  //! \param[in] n Number of spatial dimensions
  NonlinearElasticityFbar(unsigned short int n = 3, int nvp = 1);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityFbar() {}

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] nPt Number of volumetric sampling points in this element
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3&, size_t nPt);

  //! \brief Evaluates volumetric sampling point data at an interior point.
  //! \param[in] fe Finite element data of current point
  virtual bool reducedInt(const FiniteElement& fe) const;

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
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return 10+npt1; }

private:
  int    npt1;  //!< Number of volumetric sampling points in each direction
  int    pbar;  //!< Polynomial order of the internal volumetric data field
  double scale; //!< Scaling factor for extrapolation from sampling points

  mutable size_t iP;   //!< Local sampling point counter
  mutable Matrix dMdx; //!< Modified basis function gradients
  mutable std::vector<VolPtData> myVolData; //!< Volumetric sampling point data

  friend class ElasticityNormFbar;
};


/*!
  \brief Class representing the integrand of the Fbar elasticity energy norm.
*/

class ElasticityNormFbar : public ElasticityNormUL
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The elasticity problem to evaluate norms for
  ElasticityNormFbar(NonlinearElasticityFbar& p) : ElasticityNormUL(p) {}
  //! \brief Empty destructor.
  virtual ~ElasticityNormFbar() {}

  //! \brief Evaluates volumetric sampling point data at an interior point.
  //! \param[in] fe Finite element data of current point
  virtual bool reducedInt(const FiniteElement& fe) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return myProblem.getIntegrandType(); }
};

#endif
