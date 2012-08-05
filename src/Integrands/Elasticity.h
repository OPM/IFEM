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

class LocalSystem;
class Material;


/*!
  \brief Base class representing the integrand of elasticity problems.
  \details Implements common features for linear and nonlinear elasticity
  problems. None of the \a evalInt and \a evalBou methods are implemented here.
  Thus, it is regarded as an abstract base class with a protected constructor.
*/

class Elasticity : public IntegrandBase
{
protected:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] ax \e If \e true, and axisymmetric 3D formulation is assumed
  Elasticity(unsigned short int n = 3, bool ax = false);

public:
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~Elasticity();

  //! \brief Prints out the problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(TractionFunc* tf) { tracFld = tf; }
  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { fluxFld = tf; }
  //! \brief Defines the body force field.
  void setBodyForce(VecFunc* bf) { bodyFld = bf; }

  //! \brief Defines the gravitation vector.
  void setGravity(double gx, double gy = 0.0, double gz = 0.0)
  { grav[0] = gx; grav[1] = gy; grav[2] = gz; }

  //! \brief Defines the material properties.
  virtual void setMaterial(Material* mat) { material = mat; }

  //! \brief Defines the local coordinate system for stress output.
  void setLocalSystem(LocalSystem* cs) { locSys = cs; }

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  virtual bool evalSol(Vector& s, const Vector& eV,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, bool toLocal = false) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical stress values at current point
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const STensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the primary solution at a result point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \return Primary solution vector at current point
  Vec3 evalSol(const Vector& eV, const Vector& N) const;

  //! \brief Evaluates the boundary traction field (if any) at specified point.
  Vec3 getTraction(const Vec3& X, const Vec3& n) const;
  //! \brief Evaluates the body force field (if any) at specified point.
  virtual Vec3 getBodyforce(const Vec3& X) const;
  //! \brief Returns whether an external load is defined.
  virtual bool haveLoads() const;

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
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const;
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField1Name(size_t i, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

  //! \brief Prints out the maximum secondary solution values.
  //! \param os Output stream to write the values to
  //! \param[in] precision Number of digits after the decimal point
  //! \param[in] comp Which component to print (0 means all)
  void printMaxVals(std::ostream& os, std::streamsize precision,
                    size_t comp = 0) const;

protected:
  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[out] B The strain-displacement matrix
  //! \param[out] eps Strain tensor at current point
  virtual bool kinematics(const Vector& eV,
			  const Vector& N, const Matrix& dNdX, double r,
			  Matrix& B, Tensor&, SymmTensor& eps) const;

  //! \brief Calculates integration point geometric stiffness contributions.
  //! \param EM Element matrix to receive the stiffness contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[in] sigma Stress tensor at current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formKG(Matrix& EM, const Vector& N, const Matrix& dNdX,
	      double r, const Tensor& sigma, double detJW) const;

  //! \brief Calculates integration point mass matrix contributions.
  //! \param EM Element matrix to receive the mass contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formMassMatrix(Matrix& EM, const Vector& N,
		      const Vec3& X, double detJW) const;

  //! \brief Calculates integration point body force vector contributions.
  //! \param ES Element vector to receive the body force contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formBodyForce(Vector& ES, const Vector& N,
		     const Vec3& X, double detJW) const;

  //! \brief Calculates the strain-displacement matrix.
  //! \param[in] Bmat The strain-displacement matrix
  //! \param[in] dNdX Basis function gradients at current point
  bool formBmatrix(Matrix& Bmat, const Matrix& dNdX) const;
  //! \brief Calculates the axi-symmetric strain-displacement matrix.
  //! \param[in] Bmat The strain-displacement matrix
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  bool formBmatrix(Matrix& Bmat, const Vector& N, const Matrix& dNdX,
		   double r) const;

public:
  //! \brief Sets up the inverse constitutive matrix at current point.
  //! \param[out] Cinv \f$6\times6\f$-matrix (in 3D) or \f$3\times3\f$-matrix
  //! (in 2D), representing the inverse constitutive tensor
  //! \param[in] X Cartesian coordinates of current point
  bool formCinverse(Matrix& Cinv, const Vec3& X) const;

  //! \brief Returns \e true if this is an axial-symmetric problem.
  bool isAxiSymmetric() const { return axiSymmetry; }

protected:
  // Finite element quantities, i.e., indices into element matrices and vectors.
  // These indices will be identical for all elements in a model and can thus
  // be stored here, even when doing multi-threading. Note that these indices
  // 1-based, since the value zero is used to signal non-existing matrix/vector.
  unsigned short int eKm; //!< Index to element material stiffness matrix
  unsigned short int eKg; //!< Index to element geometric stiffness matrix
  unsigned short int eM;  //!< Index to element mass matrix
  unsigned short int eS;  //!< Index to element load vector
  unsigned short int iS;  //!< Index to element internal force vector

  // Physical properties
  Material* material; //!< Material data and constitutive relation
  double    grav[3];  //!< Gravitation vector

  LocalSystem*  locSys;  //!< Local coordinate system for result output
  TractionFunc* tracFld; //!< Pointer to implicit boundary traction field
  VecFunc*      fluxFld; //!< Pointer to explicit boundary traction field
  VecFunc*      bodyFld; //!< Pointer to body force field

  typedef std::pair<Vec3,double> PointValue; //!< Convenience type

  mutable std::vector<PointValue> maxVal;  //!< Maximum result values
  mutable std::vector<Vec3Pair>   tracVal; //!< Traction field point values

  unsigned short int nsd; //!< Number of space dimensions (1, 2 or 3)
  unsigned short int nDF; //!< Dimension on deformation gradient (2 or 3)
  bool       axiSymmetry; //!< \e true if the problem is axi-symmetric
};


/*!
  \brief Class representing the integrand of elasticity energy norms.
*/

class ElasticityNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  //! \param[in] a The analytical stress field (optional)
  ElasticityNorm(Elasticity& p, STensorFunc* a = 0);
  //! \brief Empty destructor.
  virtual ~ElasticityNorm() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Add external energy terms to relevant norms
  virtual void addBoundaryTerms(Vectors& gNorm, double extEnergy);

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields(int fld=0) const;

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

  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  //! \param elmInt The local integral object to receive the contributions
  virtual bool finalizeElement(LocalIntegral& elmInt, const TimeDomain&,size_t);

  virtual bool hasElementContributions(size_t i, size_t j)
  { 
    return (i == 1 && j < 3) || i > 1;
  }

  virtual const char* getName(size_t i, size_t j, const char* prefix);

private:
  STensorFunc* anasol; //!< Analytical stress field
};

#endif
