// $Id$
//==============================================================================
//!
//! \file Elasticity.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for linear and nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _ELASTICITY_H
#define _ELASTICITY_H

#include "IntegrandBase.h"
#include "Vec3.h"
#include <map>

class LocalSystem;
class ElmMats;
class ElmNorm;
class VTF;


/*!
  \brief Base class representing the integrand of elasticity problems.
  \details Implements common features for linear and nonlinear elasticity
  problems. None of the \a evalInt and \a evalBou methods are implemented here.
  Thus, it is regarded as an abstract base class with a protected constructor.
*/

class Elasticity : public Integrand
{
protected:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] ps If \e true, assume plane stress in 2D
  Elasticity(unsigned short int n = 3, bool ps = true);

public:
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~Elasticity();

  //! \brief Prints out problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(TractionFunc* tf);
  //! \brief Clears the integration point traction values.
  void clearTracVal() { tracVal.clear(); }

  //! \brief Defines the gravitation vector.
  void setGravity(double gx, double gy = 0.0, double gz = 0.0)
  { g[0] = gx; g[1] = gy; g[2] = gz; }

  //! \brief Defines the body force field.
  void setBodyForce(VecFunc* bf) { bodyFld = bf; }

  //! \brief Defines material properties for current volume patch.
  //! \param[in] Young   Young's modulus
  //! \param[in] Poiss   Poisson's ratio
  //! \param[in] Density Mass density
  virtual void setMaterial(double Young, double Poiss, double Density)
  { Emod = Young; nu = Poiss; rho = Density; }

  //! \brief Defines the local coordinate system for stress output.
  void setLocalSystem(LocalSystem* cs) { locSys = cs; }

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElement(const std::vector<int>& MNPC);
  //! \brief Initializes current element for boundary numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC);

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
  //! \param[out] s The FE stress values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const Matrix& dNdX, const Vec3& X) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical stress values at current point
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const STensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the primary solution at a result point.
  //! \param[in] N Basis function values at current point
  //! \return Primary solution vector at current point
  Vec3 evalSol(const Vector& N) const;

  //! \brief Evaluates the boundary traction field (if any) at specified point.
  Vec3 getTraction(const Vec3& X, const Vec3& n) const;
  //! \brief Evaluates the body force field (if any) at specified point.
  virtual Vec3 getBodyforce(const Vec3& X) const;
  //! \brief Returns whether an external load is defined
  virtual bool haveLoads() const;

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
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns the number of secondary solution fields.
  virtual size_t getNoFields() const;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getFieldLabel(size_t i, const char* prefix = 0) const;

protected:
  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return rho; }

  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[out] eps Strain tensor at current point
  //!
  //! \details The strain displacement matrix \b B is established
  //! and stored in the mutable class member \a Bmat.
  virtual bool kinematics(const Matrix& dNdX, SymmTensor& eps) const;

  //! \brief Evaluates the constitutive relation at current point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density
  //! \param[in] eps Strain tensor at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] calcStress If \e false, claculate the C-matrix only
  virtual bool constitutive(Matrix& C, SymmTensor& sigma, double& U,
			    const SymmTensor& eps, const Vec3& X,
			    char calcStress = true) const;

  //! \brief Calculates integration point geometric stiffness contributions.
  //! \param EM Element matrix to receive the stiffness contributions
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] sigma Stress tensor at current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formKG(Matrix& EM, const Matrix& dNdX,
	      const Tensor& sigma, double detJW) const;

  //! \brief Calculates integration point mass matrix contributions.
  //! \param EM Element matrix to receive the mass contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formMassMatrix(Matrix& EM, const Vector& N,
		      const Vec3& X, double detJW) const;

  //! \brief Calculates integration point body force vector contributions.
  //! \brief Evaluates the body force vector at current point.
  //! \param ES Element vector to receive the body force contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formBodyForce(Vector& ES, const Vector& N,
		     const Vec3& X, double detJW) const;

  //! \brief Calculates the strain-displacement matrix \b B at current point.
  //! \param[in] dNdX Basis function gradients at current point
  bool formBmatrix(const Matrix& dNdX) const;

  //! \brief Utility used by the virtual \a evalInt and \a evalBou methods.
  //! \param elmInt Pointer to the integrated element quantities
  bool getIntegralResult(LocalIntegral*& elmInt) const;

public:
  //! \brief Sets up the tangential constitutive matrix at current point.
  //! \param[out] C \f$6\times6\f$-matrix (in 3D) or \f$3\times3\f$-matrix
  //! (in 2D), representing the constitutive tensor
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool formCmatrix(Matrix& C, const Vec3& X, bool invers = false) const;

private:
  // Physical properties (constant)
  double Emod; //!< Young's modulus
  double nu;   //!< Poisson's ratio
  double rho;  //!< Mass density
  double g[3]; //!< Gravitation vector

protected:
  // Finite element quantities
  Matrix* eKm; //!< Pointer to element material stiffness matrix
  Matrix* eKg; //!< Pointer to element geometric stiffness matrix
  Matrix* eM;  //!< Pointer to element mass matrix
  Vector* eS;  //!< Pointer to element load vector
  Vector* iS;  //!< Pointer to element internal force vector
  Vector* eV;  //!< Pointer to element displacement vector

  ElmMats* myMats; //!< Local element matrices, result of numerical integration

  LocalSystem*  locSys;  //!< Local coordinate system for result output
  TractionFunc* tracFld; //!< Pointer to boundary traction field
  VecFunc*      bodyFld; //!< Pointer to body force field

  mutable std::map<Vec3,Vec3> tracVal; //!< Traction field point values

  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
  bool       planeStress; //!< Plane stress/strain option for 2D problems

  // Work arrays declared as members to avoid frequent re-allocation
  // within the numerical integration loop (for reduced overhead)
  mutable Matrix Bmat; //!< Strain-displacement matrix
  mutable Matrix Cmat; //!< Constitutive matrix
};


/*!
  \brief Class representing the integrand of linear elasticity energy norms.
*/

class ElasticityNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  //! \param[in] a The analytical stress field (optional)
  ElasticityNorm(Elasticity& p, STensorFunc* a = 0) : problem(p), anasol(a) {}
  //! \brief Empty destructor.
  virtual ~ElasticityNorm() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }
  //! \brief Returns a 1-based index of the external energy norm.
  virtual size_t indExt() const { return 2; }

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

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields() const;

protected:
  //! \brief Get the element norm object to use in the integration.
  //! \param elmInt The local integral object to receive norm contributions
  //!
  //! \details If \a elmInt is NULL or cannot be casted to a ElmNorm pointer,
  //! a local static buffer is used instead.
  static ElmNorm& getElmNormBuffer(LocalIntegral*& elmInt);

  Elasticity& problem; //!< The problem-specific data

private:
  STensorFunc* anasol; //!< Analytical stress field
};

#endif
