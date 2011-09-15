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
  KirchhoffLovePlate();

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

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElement(const std::vector<int>& MNPC);
  //! \brief Initializes current element for boundary numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC);

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return 2; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const Vector& N, const Matrix& dNdX,
		       const Matrix3D& d2NdX2, const Vec3& X,
		       const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress resultant values at current point
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol(Vector& s, const Matrix3D& d2NdX2, const Vec3& X) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical stress resultant values at current point
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const STensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the primary solution at a result point.
  //! \param[in] N Basis function values at current point
  //! \return Primary solution vector at current point
  double evalSol(const Vector& N) const;

  //! \brief Evaluates the pressure field (if any) at specified point.
  virtual double getPressure(const Vec3& X) const;
  //! \brief Returns whether an external load is defined.
  virtual bool haveLoads() const;

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
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formBodyForce(Vector& ES, const Vector& N,
		     const Vec3& X, double detJW) const;

  //! \brief Calculates the strain-displacement matrix \b B at current point.
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  bool formBmatrix(const Matrix3D& d2NdX2) const;

public:
  //! \brief Sets up the constitutive matrix at current point.
  //! \param[out] C \f$3\times3\f$-matrix, representing the constitutive tensor
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, the inverse matrix is establised instead
  bool formCmatrix(Matrix& C, const Vec3& X, bool invers = false) const;

protected:
  // Finite element quantities
  Matrix* eK; //!< Pointer to element stiffness matrix
  Matrix* eM; //!< Pointer to element mass matrix
  Vector* eS; //!< Pointer to element load vector
  Vector* eV; //!< Pointer to element displacement vector

  // Physical properties
  Material* material;  //!< Material data and constitutive relation
  double    thickness; //!< Plate thickness
  double    gravity;   //!< Gravitation constant

  LocalSystem* locSys;  //!< Local coordinate system for result output
  RealFunc*    presFld; //!< Pointer to pressure field

  // Work arrays declared as members to avoid frequent re-allocation
  // within the numerical integration loop (for reduced overhead)
  mutable Matrix Bmat; //!< Strain-displacement matrix
  mutable Matrix Cmat; //!< Constitutive matrix
  mutable Matrix CB;   //!< Result of the matrix-matrix product C*B
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
  //! \brief Returns a 1-based index of the external energy norm.
  virtual size_t indExt() const { return 2; }

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields() const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElement(const std::vector<int>& MNPC);

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return 2; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

private:
  STensorFunc* anasol; //!< Analytical stress resultant field
};

#endif
