// $Id$
//==============================================================================
//!
//! \file IntegrandBase.h
//!
//! \date Nov 11 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base classes representing FEM integrands.
//!
//==============================================================================

#ifndef _INTEGRAND_BASE_H
#define _INTEGRAND_BASE_H

#include "Integrand.h"
#include "SIMenums.h"
#include "Function.h"
#include "MatVec.h"
#include <map>

class NormBase;
class AnaSol;
class VTF;
class Field;
class Fields;


/*!
  \brief Base class representing a system level integrated quantity.
*/

class IntegrandBase : public Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  IntegrandBase() : npv(1), m_mode(SIM::INIT) {}

public:
  //! \brief Empty destructor.
  virtual ~IntegrandBase() {}

  //! \brief Prints out the problem definition to the given output stream.
  virtual void print(std::ostream&) const {}


  // Global initialization interface
  // ===============================

  //! \brief Defines the solution mode before the element assembly is started.
  virtual void setMode(SIM::SolutionMode mode) { m_mode = mode; }
  //! \brief Initializes the integrand with the number of integration points.
  //! \details This method is invoked only once during the preprocessing stage.
  virtual void initIntegration(size_t, size_t) {}
  //! \brief Initializes the integrand for a new integration loop.
  //! \details This method is invoked once before starting the numerical
  //! integration over the entire spatial domain.
  virtual void initIntegration(const TimeDomain&) {}
  //! \brief Initializes the integrand for a new result point loop.
  //! \details This method is invoked once before starting the evaluation of
  //! the secondary solution at all result sampling points, after the converged
  //! primary solution has been found.
  virtual void initResultPoints(double) {}


  // Element-level initialization interface
  // ======================================

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param[in] elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  //! The default implementation forwards to an overloaded method not taking
  //! \a X0 and \a nPt as arguments.
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3& X0, size_t nPt, LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  //! \param[in] elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC1,
			   const std::vector<int>& MNPC2, size_t n1,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  //! \param[in] elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC1,
			      const std::vector<int>& MNPC2, size_t n1,
                              LocalIntegral& elmInt);


  // Secondary solution field evaluation interface
  // =============================================

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s,
		       const Vector& N, const Matrix& dNdX,
		       const Matrix3D& d2NdX2,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] N1 Basis function values at current point, field 1
  //! \param[in] N2 Basis function values at current point, field 2
  //! \param[in] dN1dX Basis function gradients at current point, field 1
  //! \param[in] dN2dX Basis function gradients at current point, field 2
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  virtual bool evalSol(Vector& s,
		       const Vector& N1, const Vector& N2,
		       const Matrix& dN1dX, const Matrix& dN2dX, const Vec3& X,
		       const std::vector<int>& MNPC1,
		       const std::vector<int>& MNPC2) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s,
		       const TensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (symmetric tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s,
		       const STensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (vector field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s,
		       const VecFunc& asol, const Vec3& X) const;


  // Various service methods
  // =======================

  //! \brief Returns whether a mixed formulation is used.
  virtual bool mixedFormulation() const { return false; }

  //! \brief Writes surface tractions/fluxes for a given time step to VTF-file.
  virtual bool writeGlvT(VTF*, int, int&) const { return true; }

  //! \brief Returns whether there are any traction/flux values to write to VTF.
  virtual bool hasTractionValues() const { return false; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const { return 0; }

  //! \brief Returns the number of primary/secondary solution field components.
  virtual size_t getNoFields(int = 2) const { return 0; }
  //! \brief Returns the name of a primary solution field component.
  virtual const char* getField1Name(size_t, const char* = 0) const { return 0; }
  //! \brief Returns the name of a secondary solution field component.
  virtual const char* getField2Name(size_t, const char* = 0) const { return 0; }

  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions() const { return primsol.size(); }

  //! \brief Accesses the primary solution vector of current patch.
  Vector& getSolution(size_t n = 0) { return primsol[n]; }

  //! \brief Resets the primary solution vectors.
  void resetSolution();

  //! \brief Registers where we can inject a mixed-basis scalar field.
  virtual void setNamedField(const std::string&, Field*) {}
  //! \brief Registers where we can inject a mixed-basis vector field.
  virtual void setNamedFields(const std::string&, Fields*) {}

  //! \brief Returns a vector where we can store a named field.
  Vector* getNamedVector(const std::string& name) const;

protected:
  //! \brief Registers a vector to inject a named field into.
  //! \param[in] name Name of field
  //! \param[in] vec Vector to inject field into
  void registerVector(const std::string& name, Vector* vec);

private:
  std::map<std::string,Vector*> myFields; //!< Named fields of this integrand

protected:
  Vectors            primsol; //!< Primary solution vectors for current patch
  unsigned short int npv;     //!< Number of primary solution variables per node
  SIM::SolutionMode  m_mode;  //!< Current solution mode
};


typedef std::vector<LocalIntegral*> LintegralVec; //!< Local integral container


/*!
  \brief Base class representing a system level norm quantity.
*/

class NormBase : public Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  NormBase(IntegrandBase& p) : myProblem(p), nrcmp(0), lints(0) {}

public:
  //! \brief Empty destructor.
  virtual ~NormBase() {}

  //! \brief Initializes the integrand with the number of integration points.
  virtual void initIntegration(size_t, size_t) {}
  //! \brief Initializes the integrand for a new integration loop.
  virtual void initIntegration(const TimeDomain& time);
  //! \brief Sets the number of projected solutions.
  void initProjection(size_t nproj) { prjsol.resize(nproj); }
  //! \brief Sets a vector of LocalIntegrals to be used during norm integration.
  void setLocalIntegrals(LintegralVec* elementNorms) { lints = elementNorms; }

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl The element number
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                          bool neumann) const;

  //! \brief Initializes current element for numerical integration.
  virtual bool initElement(const std::vector<int>& MNPC,
			   const Vec3& X0, size_t nPt, LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical integration.
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical integration (mixed).
  virtual bool initElement(const std::vector<int>& MNPC1,
			   const std::vector<int>& MNPC2, size_t n1,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration.
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration (mixed).
  virtual bool initElementBou(const std::vector<int>& MNPC1,
			      const std::vector<int>& MNPC2, size_t n1,
                              LocalIntegral& elmInt);

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return false; }

  //! \brief Returns the number of field components.
  virtual size_t getNoFields() const { return 0; }

  //! \brief Returns the name of a norm quantity.
  static const char* getName(size_t&, bool, const char* = 0);

  //! \brief Returns whether a norm component stores element contributions.
  static bool hasElementContributions(size_t i) { return i != 1; }

  //! \brief Accesses a projected secondary solution vector of current patch.
  Vector& getProjection(size_t i);

private:
  //! \brief Initializes projected field for current element.
  bool initProjection(const std::vector<int>& MNPC, LocalIntegral& elmInt);

protected:
  IntegrandBase& myProblem; //!< The problem-specific data

  Vectors prjsol; //!< Projected secondary solution vectors for current patch

  unsigned short int nrcmp; //!< Number of projected solution components
  LintegralVec*      lints; //!< Local integrals used during norm integration
};

#endif
