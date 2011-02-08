// $Id: Stokes.h,v 1.4 2011-02-08 12:21:31 rho Exp $
//==============================================================================
//!
//! \file Stokes.h
//!
//! \date Sep 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Base class for Integrand implementation of Stokes problems.
//!
//==============================================================================

#ifndef _STOKES_H
#define _STOKES_H

#include "IntegrandBase.h"
#include "ElmMats.h"
#include "Vec3.h"
#include "AnaSol.h"
#include <map>

class VTF;


/*!
  \brief Base class representing the integrand of Stokes problems.
  \details Implements common features for Stokes/Navier-Stokes problems.
*/

class Stokes : public Integrand
{
public:
  //TODO: Temporary; remove and use SIM-enum in all sub-classes instead (kmo).
  enum ProblemFormulation { LAPLACE = SIM::LAPLACE, STRESS = SIM::STRESS };
  //END TODO:

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] form The solution formulation to use
  //! \param[in] itg_type The integrand type to use
  Stokes(unsigned short int n,
         SIM::Formulation form = SIM::LAPLACE,
         int itg_type = 1);
  //! \brief The destructor frees dynamically allocated objects.
  virtual ~Stokes();

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(TractionFunc* tf) { tracFld = tf; }
  //! \brief Clears the integration point traction values.
  void clearTracVal() { tracVal.clear(); }

  //! \brief Defines the gravitation vector.
  void setGravity(double gx, double gy, double gz)
  { g[0] = gx; g[1] = gy; g[2] = gz; }

  //! \brief Defines fluid properties for current volume patch.
  //! \param[in] Density   Mass density
  //! \param[in] Viscosity Dynamic viscosity
  void setFluidProperties(double Density, double Viscosity)
  { rho = Density; mu = Viscosity; }

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElement(const std::vector<int>& MNPC);
  //! \brief Initializes current element for boundary numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC);

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

  //! \brief Evaluates the secondary solution at current integration point.
  //! \param[out] s The solution field values
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] MNPC Matrix of nodal point correspondance
  virtual bool evalSol(Vector& s,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the analytical secondary solution at the given point.
  //! \param[out] s The solution field values
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalSol(Vector& s, const TensorFunc& asol, const Vec3& X) const;

  //! \brief Writes the surface tractions for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvT(VTF* vtf, int iStep, int& nBlock) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] asol Pointer to analytical solution field (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns problem formulation type.
  SIM::Formulation getFormulation() const { return formulation; }

  //! \brief Returns which integrand to use.
  virtual int getIntegrandType() const { return integrandType; }

  //! \brief Returns the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return rho; }

  //! \brief Returns the body force per unit mass
  //! \param[in]  X Cartesian coordinate of current integration point
  //! \param[out] f Body force vector per unit mass 
  virtual bool getBodyForce(const Vec3& X, Vector& f) const;

  //! \brief Returns viscosity at current point.
  virtual double getViscosity(const Vec3&) const { return mu; }

  //! \brief Returns a pointer to current element solution vector.
  const Vector* getElementSolution() const { return eVs[0]; }

  //! \brief Returns a pointer to current element solution vector.
  virtual const Vector* getElementVelocity() const { return 0; }

  //! \brief Returns a pointer to current element solution vector.
  virtual const Vector* getElementPressure() const { return 0; }

  //! \brief Returns the number of solution fields.
  virtual size_t getNoFields() const { return nsd; }

  //! \brief Returns the number of secondary solution fields.
  virtual size_t getNoSecFields() const { return nsd*nsd; }

  //! \brief Returns the name of the primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getFieldLabel(size_t i, const char* prefix = 0) const;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getSecFieldLabel(size_t i, const char* prefix = 0) const
  { return 0; }

protected:
  //! \brief Utility used by the virtual \a evalInt and \a evalBou methods.
  //! \param elmInt Pointer to the integrated element quantities
  bool getIntegralResult(LocalIntegral*& elmInt) const;

  // Problem parameters
  unsigned short int nsd;           //!< Number of space dimensions (1, 2 or, 3)
  unsigned short int nf;            //!< Number of primary field variables
  SIM::Formulation   formulation;   //!< Problem formulation flag
  int                integrandType; //!< Integrand type flag

  // Physical properties (constant)
  double mu;   //!< Dynamic viscosity
  double rho;  //!< Mass density
  double g[3]; //!< Gravitation vector

  // Finite element quantities
  Matrix* eM;               //!< Pointer to element coefficient matrix
  Vector* eS;               //!< Pointer to element right-hand-side vector
  std::vector<Vector*> eVs; //!< Pointers to element solution vectors

  ElmMats* myMats;       //!< Local element matrices

  TractionFunc* tracFld; //!< Pointer to boundary traction field

  mutable std::map<Vec3,Vec3> tracVal; //!< Traction field point values
};


/*!
  \brief Class representing the integrand of Stokes energy norms.
*/

class StokesNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Stokes problem to evaluate norms for
  //! \param[in] a The analytical velocity and pressure fields (optional)
 StokesNorm(Stokes& p, AnaSol* a) : problem(p), anasol(a) {}
  //! \brief Empty destructor.
  virtual ~StokesNorm() {}

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

  //! \brief Returns the number of primary norm quantities.
  virtual size_t getNoFields() const { return anasol ? 2 : 0; }

  //! \brief Returns the number of secondary norm quantities.
  virtual size_t getNoSecFields() const { return anasol ? 3 : 1; }

protected:
  Stokes&     problem; //!< The problem-specific data
  AnaSol*     anasol;  //!< Analytical solution fields
};


#endif
