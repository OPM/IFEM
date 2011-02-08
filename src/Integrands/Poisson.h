// $Id: Poisson.h,v 1.4 2010-12-27 18:29:02 kmo Exp $
//==============================================================================
//!
//! \file Poisson.h
//!
//! \date Apr 16 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Integrand implementations for Poisson problems.
//!
//==============================================================================

#ifndef _POISSON_H
#define _POISSON_H

#include "IntegrandBase.h"
#include "ElmMats.h"
#include "Vec3.h"
#include <map>

class VTF;


/*!
  \brief Class representing the integrand of the Poisson problem.
  \details This class supports constant isotrophic conductivity properties only.
  Properties with spatial variation has to be implemented as sub-classes.
*/

class Poisson : public Integrand
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  Poisson(unsigned short int n = 3);
  //! \brief Empty destructor.
  virtual ~Poisson() {}

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { tracFld = tf; myMats.rhsOnly = true; }
  //! \brief Clears the integration point traction values.
  void clearTracVal() { tracVal.clear(); }

  //! \brief Defines the heat source.
  void setSource(RealFunc* src) { heatSrc = src; }

  //! \brief Defines the conductivity.
  void setMaterial(double K) { kappa = K; }

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElement(const std::vector<int>& MNPC);
  //! \brief Initializes current element for boundary numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const Vec3& normal) const;

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
  //! \param[out] s The FE solution values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const Matrix& dNdX, const Vec3& X) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical solution values at current point
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSolScal(Vector& s, const VecFunc& asol, const Vec3& X) const;

  //! \brief Writes the surface tractions for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvT(VTF* vtf, int iStep, int& nBlock) const;

  //! \brief Returns whether there are any traction values to write to VTF.
  virtual bool hasTractionValues() const { return !tracVal.empty(); }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution field (optional)
  virtual NormBase* getNormIntegrandScal(VecFunc* asol = 0) const;

  //! \brief Returns the number of secondary solution fields.
  virtual size_t getNoFields() const { return nsd; }

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getFieldLabel(size_t i, const char* prefix = 0) const;

protected:
  //! \brief Returns the conductivety at current point.
  virtual double getConductivety(const Vec3&) const { return kappa; }

public:
  //! \brief Sets up the conductivity matrix at current point.
  //! \param[out] C \f$ nsd\times nsd\f$-matrix
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool formCmatrix(Matrix& C, const Vec3& X, bool invers = false) const;

private:
  // Physical properties (constant)
  double kappa; //!< Conductivity

protected:
  // Finite element quantities
  Matrix* eM; //!< Pointer to element matrix
  Vector* eS; //!< Pointer to element right-hand-side vector
  Vector* eV; //!< Pointer to element solution vector

  mutable ElmMats myMats; //!< Local element matrices

  VecFunc*  tracFld; //!< Pointer to boundary traction field
  RealFunc* heatSrc; //!< Pointer to interior heat source

  mutable std::map<Vec3,Vec3> tracVal; //!< Traction field point values

  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)

  // Work arrays declared as members to avoid frequent re-allocation
  // within the numerical integration loop (for reduced overhead)
  mutable Matrix C;  //!< Conductivity matrix
  mutable Matrix CB; //!< Result of the matrix-matrix product C*dNdX^T
};


/*!
  \brief Class representing the integrand of Poisson energy norms.
*/

class PoissonNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] a The analytical heat flux (optional)
  PoissonNorm(Poisson& p, VecFunc* a) : problem(p), anasol(a) {}
  //! \brief Empty destructor.
  virtual ~PoissonNorm() {}

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElement(const std::vector<int>& MNPC);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const;

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields() const { return anasol ? 3 : 1; }

private:
  Poisson& problem; //!< The problem-specific data
  VecFunc* anasol;  //!< Analytical heat flux
};

#endif
