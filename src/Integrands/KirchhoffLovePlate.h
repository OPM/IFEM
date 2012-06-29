// $Id$
//==============================================================================
//!
//! \file KirchhoffLovePlate.h
//!
//! \date Sep 13 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for linear Kirchhoff-Love thin plate problems.
//!
//==============================================================================

#ifndef _KIRCHHOFF_LOVE_PLATE_H
#define _KIRCHHOFF_LOVE_PLATE_H

#include "IntegrandBase.h"
#include "Vec3.h"

class LocalSystem;
class Material;


/*!
  \brief Class representing the integrand of thin plate problems.
  \details The formulation is based on Kirchhoff-Love plate theory and therefore
  requires second-derivatives of the basis functions.
*/

class KirchhoffLovePlate : public IntegrandBase
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions (1=beam, 2=plate)
  KirchhoffLovePlate(unsigned short int n = 2);
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~KirchhoffLovePlate();

  //! \brief Prints out the problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Defines the gravitation constant.
  void setGravity(double g) { gravity = g; }

  //! \brief Defines the gravitation constant.
  void setThickness(double t) { thickness = t; }

  //! \brief Defines the pressure field.
  void setPressure(RealFunc* pf) { presFld = pf; }

  //! \brief Defines the material properties.
  void setMaterial(Material* mat) { material = mat; }

  //! \brief Defines the local coordinate system for stress resultant output.
  void setLocalSystem(LocalSystem* cs) { locSys = cs; }

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const
  { 
    return Integrand::SECOND_DERIVATIVES;
  }

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

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
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const Vector&, const Matrix&,
		       const Matrix3D& d2NdX2, const Vec3& X,
		       const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress resultant values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  bool evalSol(Vector& s, const Vector& eV,
	       const Matrix3D& d2NdX2, const Vec3& X,
	       bool toLocal = false) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical stress resultant values at current point
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const STensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the pressure field (if any) at specified point.
  virtual double getPressure(const Vec3& X) const;
  //! \brief Returns whether an external load is defined.
  virtual bool haveLoads() const;

  //! \brief Writes the surface pressure for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the pressure vectors
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& nBlock) const;
  //! \brief Returns whether there are any pressure values to write to VTF.
  virtual bool hasTractionValues() const { return !presVal.empty(); }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const;
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual const char* getField1Name(size_t, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

protected:
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
  //! \param[in] iP Global integration point counter
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formBodyForce(Vector& ES, const Vector& N,
		     size_t iP, const Vec3& X, double detJW) const;

  //! \brief Calculates the strain-displacement matrix \b B at current point.
  //! \param[out] Bmat The strain-displacement matrix
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  bool formBmatrix(Matrix& Bmat, const Matrix3D& d2NdX2) const;

public:
  //! \brief Sets up the constitutive matrix at current point.
  //! \param[out] C \f$3\times3\f$-matrix, representing the constitutive tensor
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, the inverse matrix is establised instead
  bool formCmatrix(Matrix& C, const Vec3& X, bool invers = false) const;

protected:
  // Finite element quantities, i.e., indices into element matrices and vectors.
  // These indices will be identical for all elements in a model and can thus
  // be stored here, even when doing multi-threading. Note that these indices
  // 1-based, since the value zero is used to signal non-existing matrix/vector.
  unsigned short int eK; //!< Index to element stiffness matrix
  unsigned short int eM; //!< Index to element mass matrix
  unsigned short int eS; //!< Index to element load vector

  // Physical properties
  Material* material;  //!< Material data and constitutive relation
  double    thickness; //!< Plate thickness
  double    gravity;   //!< Gravitation constant

  LocalSystem* locSys;  //!< Local coordinate system for result output
  RealFunc*    presFld; //!< Pointer to pressure field

  mutable std::vector<Vec3Pair> presVal; //!< Pressure field point values

  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
};


/*!
  \brief Class representing the integrand of Kirchhoff-Love energy norms.
*/

class KirchhoffLovePlateNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  //! \param[in] a The analytical stress resultant field (optional)
  KirchhoffLovePlateNorm(KirchhoffLovePlate& p, STensorFunc* a = 0);
  //! \brief Empty destructor.
  virtual ~KirchhoffLovePlateNorm() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields(int fld=0) const;

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const
  { 
    return Integrand::SECOND_DERIVATIVES;
  }

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

  //! \brief Returns whether or not the element norm contributions should
  //         be stored for visualization
  virtual bool hasElementContributions(size_t i, size_t j)
  { 
    return (i == 1 && j < 2) || j < 3;
  }

  //! \brief Return the name of a particular norm identified by group and entry
  const char* getName(size_t i, size_t j, const char* prefix);

  //! \brief Add external energy terms to relevant norms
  void addBoundaryTerms(Vectors& gNorm, double extEnergy);

private:
  STensorFunc* anasol; //!< Analytical stress resultant field
};

#endif
