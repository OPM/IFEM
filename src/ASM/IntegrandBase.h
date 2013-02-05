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
#include "ASMenums.h"
#include "Function.h"
#include "MatVec.h"
#include <map>

class NormBase;
class ForceBase;
class GlobalIntegral;
class GlbNorm;
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
  virtual void initIntegration(const TimeDomain&, const Vector&, bool=false) {}
  //! \brief Initializes the integrand for a new result point loop.
  //! \details This method is invoked once before starting the evaluation of
  //! the secondary solution at all result sampling points, after the converged
  //! primary solution has been found.
  virtual void initResultPoints(double) {}
  //! \brief Initializes the global node number mapping for current patch.
  virtual void initNodeMap(const std::vector<int>&) {}
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const { return *gq; }


  // Element-level initialization interface
  // ======================================

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const Vec3& X0, size_t nPt, LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on current patch
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC1,
                           const std::vector<int>& MNPC2, size_t n1,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on current patch
  //! \param elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC1,
                              const std::vector<int>& MNPC2, size_t n1,
                              LocalIntegral& elmInt);

  //! \brief Returns whether this integrand has explicit interior contributions.
  virtual bool hasInteriorTerms() const { return true; }
  //! \brief Returns whether this integrand has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }


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

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
		       const std::vector<int>& MNPC) const;

  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Mixed finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  virtual bool evalSol(Vector& s, const MxFiniteElement& fe, const Vec3& X,
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

  //! \brief Returns \e true if simulation diverged on integration point level.
  virtual bool diverged(size_t = 0) const { return false; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const { return 0; }
  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  virtual ForceBase* getForceIntegrand(const Vec3* = 0, AnaSol* = 0) const
  { return 0; }

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
  NormBase(IntegrandBase& p) : myProblem(p), nrcmp(0), lints(0),
                               finalOp(ASM::SQRT) {}

public:
  //! \brief Empty destructor.
  virtual ~NormBase() {}

  //! \brief Initializes the integrand with the number of integration points.
  virtual void initIntegration(size_t, size_t) {}
  //! \brief Sets the number of projected solutions.
  void initProjection(size_t nproj) { prjsol.resize(nproj); }
  //! \brief Sets a vector of LocalIntegrals to be used during norm integration.
  void setLocalIntegrals(LintegralVec* elementNorms) { lints = elementNorms; }

  //! \brief Returns a local integral container for the element \a iEl.
  virtual LocalIntegral* getLocalIntegral(size_t, size_t iEl, bool) const;

  //! \brief Initializes current element for numerical integration.
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration.
  virtual bool initElement(const std::vector<int>& MNPC,
                           const Vec3& X0, size_t nPt, LocalIntegral& elmInt);
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

  //! \brief Adds external energy terms to relevant norms.
  virtual void addBoundaryTerms(Vectors&, double) const {}

  //! \brief Returns the number of norm groups or size of a specified group.
  virtual size_t getNoFields(int group = 0) const { return 0; }

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group
  //! \param[in] j The norm number
  //! \param[in] prefix Common prefix for all norm names
  virtual const char* getName(size_t i, size_t j, const char* prefix = 0) const;

  //! \brief Returns whether a norm quantity stores element contributions.
  virtual bool hasElementContributions(size_t, size_t) const { return true; }

  //! \brief Accesses a projected secondary solution vector of current patch.
  Vector& getProjection(size_t i);

  //! \brief Sets the final operation to apply to norms.
  void setFinalOperation(ASM::FinalNormOp op) { finalOp = op; }
  //! \brief Returns the final operation applied to norms.
  ASM::FinalNormOp getFinalOperation() { return finalOp; }

protected:
  //! \brief Initializes the projected fields for current element.
  bool initProjection(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  IntegrandBase& myProblem; //!< The problem-specific data

  Vectors prjsol; //!< Projected secondary solution vectors for current patch

  unsigned short int nrcmp;   //!< Number of projected solution components
  LintegralVec*      lints;   //!< Local integrals used during norm integration
  ASM::FinalNormOp   finalOp; //!< The final operation to apply to norms
};


/*!
  \brief Base class representing a system level boundary force quantity.
*/

class ForceBase : public Integrand
{
protected:
  //! \brief The constructor is protected to allow sub-classes only.
  ForceBase(IntegrandBase& p) : myProblem(p), eBuffer(NULL) {}

public:
  //! \brief The destructor frees the internally allocated objects.
  virtual ~ForceBase();

  //! \brief Allocates internal element force buffers.
  bool initBuffer(size_t nel);

  //! \brief Assembles the global forces.
  void assemble(RealArray& force) const;

  //! \brief Initializes the integrand with the number of integration points.
  virtual void initIntegration(size_t, size_t) {}

  //! \brief Returns a local integral container for the element \a iEl.
  virtual LocalIntegral* getLocalIntegral(size_t, size_t iEl, bool) const;

  //! \brief Dummy implementation (only boundary integration is relevant).
  virtual bool initElement(const std::vector<int>&, LocalIntegral&)
  { return false; }

  //! \brief Dummy implementation (only boundary integration is relevant).
  virtual bool initElement(const std::vector<int>&, const Vec3&,
                           size_t, LocalIntegral&)
  { return false; }

  //! \brief Dummy implementation (only boundary integration is relevant).
  virtual bool initElement(const std::vector<int>&, const std::vector<int>&,
                           size_t, LocalIntegral&)
  { return false; }

  //! \brief Initializes current element for boundary integration.
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);
  //! \brief Initializes current element for boundary integration (mixed).
  virtual bool initElementBou(const std::vector<int>& MNPC1,
                              const std::vector<int>& MNPC2, size_t n1,
                              LocalIntegral& elmInt);

  //! \brief Returns the number of force components.
  virtual size_t getNoComps() const = 0;

protected:
  IntegrandBase& myProblem; //!< The problem-specific data
  LintegralVec      eForce; //!< Local integrals used during force integration
  double*          eBuffer; //!< Element force buffer used during integration
};

#endif
