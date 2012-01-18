// $Id$
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
#include "Vec3.h"
#include <map>

class ElmNorm;


/*!
  \brief Class representing the integrand of the Poisson problem.
  \details This class supports constant isotropic constitutive properties only.
  Properties with spatial variation has to be implemented as sub-classes.
*/

class Poisson : public IntegrandBase
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
  void setTraction(VecFunc* tf) { tracFld = tf; }
  //! \brief Clears the integration point traction values.
  void clearTracVal() { tracVal.clear(); }

  //! \brief Defines the conductivity (constitutive property).
  void setMaterial(double K) { kappa = K; }

  //! \brief Defines the heat source.
  void setSource(RealFunc* src) { heatSrc = src; }

  //! \brief Evaluates the boundary traction field (if any) at specified point.
  double getTraction(const Vec3& X, const Vec3& n) const;
  //! \brief Evaluates the heat source (if any) at specified point.
  double getHeat(const Vec3& X) const;

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
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
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
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
  //! \param[in] eV Element solution vector
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const Vector& eV,
                       const Matrix& dNdX, const Vec3& X) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical solution values at current point
  //! \param[in] asol The analytical solution field (heat flux)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const VecFunc& asol, const Vec3& X) const;

  //! \brief Writes the surface tractions for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& nBlock) const;

  //! \brief Returns whether there are any traction values to write to VTF.
  virtual bool hasTractionValues() const { return !tracVal.empty(); }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld > 1 ? nsd : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual const char* getField1Name(size_t, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

  //! \brief Sets up the constitutive matrix at current point.
  //! \param[out] C \f$ nsd\times nsd\f$-matrix
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool formCmatrix(Matrix& C, const Vec3& X, bool invers = false) const;

private:
  // Physical properties (constant)
  double kappa; //!< Conductivity

protected:
  VecFunc*  tracFld; //!< Pointer to boundary traction field
  RealFunc* heatSrc; //!< Pointer to interior heat source

  mutable std::map<Vec3,Vec3> tracVal; //!< Traction field point values

  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
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
  PoissonNorm(Poisson& p, VecFunc* a = 0) : NormBase(p), anasol(a) {}
  //! \brief Empty destructor.
  virtual ~PoissonNorm() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields() const { return anasol ? 4 : 2; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite Element quantities
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt,  const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

private:
  VecFunc* anasol; //!< Analytical heat flux
};

#endif
