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
#include "ASMenums.h"
#include "LinAlgenums.h"
#include "MatVec.h"
#include <map>

class ASMbase;
class NormBase;
class ForceBase;
class GlobalIntegral;
class AnaSol;
class VTF;
class Field;
class Fields;
class VecFunc;
class TensorFunc;
class STensorFunc;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Base class representing a system level integrated quantity.
*/

class IntegrandBase : public Integrand
{
protected:
  //! \brief The constructor is protected to allow sub-classes only.
  explicit IntegrandBase(unsigned short int n) : nsd(n), npv(1),
                                                 m_mode(SIM::INIT) {}

public:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement*) { return false; }

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const {}


  // Global initialization interface
  // ===============================

  //! \brief Defines the solution mode before the element assembly is started.
  virtual void setMode(SIM::SolutionMode mode);
  //! \brief Returns current solution mode.
  SIM::SolutionMode getMode(bool = false) const override { return m_mode; }
  //! \brief Initializes an integration parameter for the integrand.
  virtual void setIntegrationPrm(unsigned short int, double) {}
  //! \brief Returns an integration parameter for the integrand.
  virtual double getIntegrationPrm(unsigned short int) const { return 0.0; }
  //! \brief Initializes the integrand with the number of integration points.
  //! \details This method is invoked only once during the preprocessing stage.
  virtual void initIntegration(size_t, size_t) {}
  //! \brief Initializes the integrand for a new integration loop.
  //! \details This method is invoked once before starting the numerical
  //! integration over the entire spatial domain.
  virtual void initIntegration(const TimeDomain&, const Vector&, bool = false){}
  //! \brief Initializes and toggles the use of left-hand-side matrix buffers.
  virtual void initLHSbuffers(size_t) {}
  //! \brief Initializes the integrand for a new result point loop.
  //! \details This method is invoked once before starting the evaluation of
  //! the secondary solution at all result sampling points, after the converged
  //! primary solution has been found. It is reimplemented for integrands
  //! containing internal result buffers that need to be (re-)initialized.
  virtual void initResultPoints(double, bool = false) {}
  //! \brief Initializes the global node number mapping for current patch.
  virtual void initNodeMap(const std::vector<int>&) {}
  //! \brief Assigns a secondary integral to be computed (for reaction forces).
  virtual void setSecondaryInt(GlobalIntegral* = nullptr) {}
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const;
  //! \brief Interface for initialization of integrand with patch-specific data.
  virtual void initForPatch(const ASMbase* pch);


  // Element-level initialization interface
  // ======================================

  using Integrand::getLocalIntegral;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl Global element number (1-based)
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                  bool neumann) const override;

protected:
  //! \brief Initializes the first primary solution vector for current element.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[out] elmVec Primary element solution vector
  //! \param[in] nskip If nonzero, skip the last \a nskip nodes of \a MNPC
  bool initElement1(const std::vector<int>& MNPC, Vectors& elmVec,
                    size_t nskip = 0) const;
  //! \brief Initializes all primary solution vectors for current element.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[out] elmVec Primary element solution vectors
  //! \param[in] nskip If nonzero, skip the last \a nskip nodes of \a MNPC
  bool initElement2(const std::vector<int>& MNPC, Vectors& elmVec,
                    size_t nskip = 0) const;

public:
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  bool initElement(const std::vector<int>& MNPC,
                   LocalIntegral& elmInt) override;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC,
                   const FiniteElement& fe,
                   const Vec3& X0, size_t nPt, LocalIntegral& elmInt) override;
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC,
                   const std::vector<size_t>& elem_sizes,
                   const std::vector<size_t>& basis_sizes,
                   LocalIntegral& elmInt) override;
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC,
                   const MxFiniteElement& fe,
                   const std::vector<size_t>& elem_sizes,
                   const std::vector<size_t>& basis_sizes,
                   LocalIntegral& elmInt) override;

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  bool initElementBou(const std::vector<int>& MNPC,
                      LocalIntegral& elmInt) override;
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt Local integral for element
  bool initElementBou(const std::vector<int>& MNPC,
                      const std::vector<size_t>& elem_sizes,
                      const std::vector<size_t>& basis_sizes,
                      LocalIntegral& elmInt) override;

  //! \brief Returns whether this integrand has explicit interior contributions.
  virtual bool hasInteriorTerms() const { return true; }
  //! \brief Returns whether this integrand has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Assigns the group of active elements.
  void activateElmGroup(const std::vector<int>& elms = {}) { elmGrp = elms; }
  //! \brief Returns \e true, if the element \a iel is deactivated.
  bool inActive(int iel) const;


  // Secondary solution field evaluation interface
  // =============================================

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  //! \param[in] nskip If nonzero, skip the last \a nskip nodes of \a MNPC
  bool evalSol1(Vector& s, const FiniteElement& fe, const Vec3& X,
                const std::vector<int>& MNPC, size_t nskip = 0) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] elmVec Element-level primary solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol2(Vector& s, const Vectors& elmVec,
                        const FiniteElement& fe, const Vec3& X) const;

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
  //! \param[in] MNPC Nodal point correspondance for the bases
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  virtual bool evalSol(Vector& s, const MxFiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC,
                       const std::vector<size_t>& elem_sizes,
                       const std::vector<size_t>& basis_sizes) const;

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

  //! \brief Returns an evaluated principal direction vector field for plotting.
  virtual bool getPrincipalDir(Matrix&, size_t, size_t) const { return false; }

  //! \brief Returns max number of 2ndary solution components to print per line.
  virtual size_t getNo2ndSolPerLine() const { return 999; }

  //! \brief Computes some derived primary solution quantities.
  virtual void primaryScalarFields(Matrix&) {}


  // Various service methods
  // =======================

  //! \brief Returns the derivative order of the differential operator.
  virtual int derivativeOrder() const { return 1; }

  //! \brief Writes surface tractions/fluxes for a given time step to VTF-file.
  virtual bool writeGlvT(VTF*, int, int&, int&) const { return true; }

  //! \brief Returns whether there are any traction/flux values to write to VTF.
  virtual bool hasTractionValues() const { return false; }

  //! \brief Returns \e true if simulation diverged on integration point level.
  virtual bool diverged(size_t = 0) const { return false; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  virtual NormBase* getNormIntegrand(AnaSol* = nullptr) const
  { return nullptr; }
  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  virtual ForceBase* getForceIntegrand(const Vec3*, AnaSol* = nullptr) const
  { return nullptr; }
  //! \brief Returns a pointer to an Integrand for nodal force evaluation.
  virtual ForceBase* getForceIntegrand() const { return nullptr; }

  //! \brief Returns the number of spatial dimensions.
  size_t getNoSpaceDim() const { return nsd; }
  //! \brief Returns the number of primary/secondary solution field components.
  virtual size_t getNoFields(int = 2) const { return 0; }
  //! \brief Returns the number of global %Lagrange multipliers in the model.
  virtual size_t getNoGLMs() const { return 0; }

  //! \brief Returns the name of a primary solution field component.
  //! \param[in] idx Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t idx, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] idx Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t idx, const char* prefix = 0) const;
  //! \brief Filters a result components for output.
  virtual bool suppressOutput(size_t, ASM::ResultClass) const { return false; }

  //! \brief Returns the number of solution vectors.
  virtual size_t getNoSolutions(bool = true) const { return primsol.size(); }

  //! \brief Returns the patch-wise extraction function field, if any.
  virtual Vector* getExtractionField(size_t = 1) { return nullptr; }
  //! \brief Accesses the primary solution vector of current patch.
  Vector& getSolution(size_t n = 0) { return primsol[n]; }
  //! \brief Accesses the primary solution vectors of current patch.
  Vectors& getSolutions() { return primsol; }

  //! \brief Resets the primary solution vectors.
  void resetSolution();

  //! \brief Prints out the patch-wise solution vectors.
  void printSolution(std::ostream& os, int pindx);

  //! \brief Registers where we can inject a mixed-basis scalar field.
  virtual void setNamedField(const std::string&, Field*);
  //! \brief Registers where we can inject a mixed-basis vector field.
  virtual void setNamedFields(const std::string&, Fields*);

  //! \brief Returns a vector where we can store a named field.
  Vector* getNamedVector(const std::string& name) const;

  //! \brief Defines the properties of the resulting linear system.
  //! \details This method is used by PETSc to optimize assembly and
  //! matrix-vector products. For maximum speed always override this method
  //! to reflect symmetry/definiteness of your operator.
  virtual LinAlg::LinearSystemType getLinearSystemType() const
  {
    return LinAlg::GENERAL_MATRIX;
  }

  //! \brief Registers a vector to inject a named field into.
  //! \param[in] name Name of field
  //! \param[in] vec Vector to inject field into
  void registerVector(const std::string& name, Vector* vec);

  //! \brief Returns nodal DOF flags for monolithic coupled integrands.
  virtual void getNodalDofTypes(std::vector<char>&) const {}

private:
  std::map<std::string,Vector*> myFields; //!< Named fields of this integrand

protected:
  unsigned short int nsd;     //!< Number of spatial dimensions (1, 2 or 3)
  unsigned short int npv;     //!< Number of primary solution variables per node
  SIM::SolutionMode  m_mode;  //!< Current solution mode
  std::vector<int>   elmGrp;  //!< List of currently active elements
  Vectors            primsol; //!< Primary solution vectors for current patch
};


using LintegralVec = std::vector<LocalIntegral*>; //!< Local integral container


/*!
  \brief Base class representing a system level norm quantity.
*/

class NormBase : public Integrand
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  explicit NormBase(IntegrandBase& p) : myProblem(p), projBou(false), nrcmp(0),
                                        lints(nullptr), finalOp(ASM::SQRT) {}

public:
  //! \brief The destructor deletes the projected secondary solution fields.
  virtual ~NormBase();

  //! \brief Initializes the integrand with the number of integration points.
  virtual void initIntegration(size_t, size_t) {}
  //! \brief Sets the number of projected solutions.
  void initProjection(size_t nproj);
  //! \brief Sets a vector of LocalIntegrals to be used during norm integration.
  void setLocalIntegrals(LintegralVec* elementNorms) { lints = elementNorms; }

  using Integrand::getLocalIntegral;
  //! \brief Returns a local integral container for the element \a iEl.
  LocalIntegral* getLocalIntegral(size_t, size_t iEl, bool) const override;

  //! \brief Initializes current element for numerical integration.
  bool initElement(const std::vector<int>& MNPC,
                   LocalIntegral& elmInt) override;
  //! \brief Initializes current element for numerical integration.
  bool initElement(const std::vector<int>& MNPC,
                   const FiniteElement& fe, const Vec3& X0, size_t nPt,
                   LocalIntegral& elmInt) override;
  //! \brief Initializes current element for numerical integration (mixed).
  bool initElement(const std::vector<int>& MNPC,
                   const std::vector<size_t>& elem_sizes,
                   const std::vector<size_t>& basis_sizes,
                   LocalIntegral& elmInt) override;
  //! \brief Initializes current element for numerical integration (mixed).
  bool initElement(const std::vector<int>& MNPC,
                   const MxFiniteElement&,
                   const std::vector<size_t>& elem_sizes,
                   const std::vector<size_t>& basis_sizes,
                   LocalIntegral& elmInt) override
  { return this->initElement(MNPC,elem_sizes,basis_sizes,elmInt); }


  //! \brief Initializes current element for boundary integration.
  bool initElementBou(const std::vector<int>& MNPC,
                      LocalIntegral& elmInt) override;
  //! \brief Initializes current element for boundary integration (mixed).
  bool initElementBou(const std::vector<int>& MNPC,
                      const std::vector<size_t>& elem_sizes,
                      const std::vector<size_t>& basis_sizes,
                      LocalIntegral& elmInt) override;

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return false; }

  //! \brief Assigns the group of active elements.
  void activateElmGroup(const std::vector<int>& elms = {})
  {
    myProblem.activateElmGroup(elms);
  }

  //! \brief Adds external energy terms to relevant norms.
  //! \param gNorm Global norm quantities
  //! \param[in] energy Global external energy
  void addBoundaryTerms(Vectors& gNorm, double energy) const;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \details If \a group is zero, the number of norm groups is returned.
  //! If \a group is greater than zero, the size of that groups is returned.
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group
  //! \param[in] j The norm number
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix = 0) const;

  //! \brief Returns whether a norm quantity stores element contributions.
  virtual bool hasElementContributions(size_t, size_t) const { return true; }

  //! \brief Accesses a projected secondary solution vector of current patch.
  Vector& getProjection(size_t i);

  //! \brief Sets the final operation to apply to norms.
  void setFinalOperation(ASM::FinalNormOp op) { finalOp = op; }
  //! \brief Returns the final operation applied to norms.
  ASM::FinalNormOp getFinalOperation() { return finalOp; }

  //! \brief Defines which FE quantities are needed by the integrand.
  int getIntegrandType() const override;
  //! \brief Returns the number of reduced-order integration points.
  int getReducedIntegration(int n) const override;

  //! \brief Evaluates reduced integration terms at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool reducedInt(LocalIntegral& elmInt,
                  const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns whether projections are fed through external means.
  virtual bool hasExternalProjections() const { return false; }

  //! \brief Sets a projected secondary solution as a field quantity.
  //! \param[in] f The field defining the projected secondary solution
  //! \param[in] idx Projection index
  virtual void setProjectedFields(Fields* f, size_t idx);

  //! \brief Assigns a parameter value to property functions of the integrand.
  //! \param[in] name Parameter name
  //! \param[in] value Parameter value
  void setParam(const std::string& name, double value) override
  {
    myProblem.setParam(name,value);
  }
  //! \brief Assigns parameter values to property functions of the integrand.
  //! \param[in] name Parameter name
  //! \param[in] value Parameter value
  void setParam(const std::string& name, const Vec3& value) override
  {
    myProblem.setParam(name,value);
  }

  //! \brief Returns current solution mode.
  SIM::SolutionMode getMode(bool) const override { return myProblem.getMode(); }

protected:
  //! \brief Initializes the projected fields for current element.
  bool initProjection(const std::vector<int>& MNPC, LocalIntegral& elmInt,
                      size_t nExtraNodes = 0);

  //! \brief Applies the operation \a finalOp on the given \a value.
  double applyFinalOp(double value) const;

  IntegrandBase& myProblem; //!< The problem-specific data

  std::vector<Fields*> prjFld; //!< Projected secondary solution fields

  Vectors prjsol; //!< Projected secondary solution vectors for current patch
  bool   projBou; //!< If \e true, the boundary integrand needs prjsol too

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
  explicit ForceBase(IntegrandBase& p) : myProblem(p), eBuffer(nullptr) {}

public:
  //! \brief The destructor frees the internally allocated objects.
  virtual ~ForceBase();

  //! \brief Allocates internal element force buffers.
  void initBuffer(size_t nel);

  //! \brief Assembles the global forces.
  void assemble(RealArray& force) const;

  //! \brief Initializes the integrand with the number of integration points.
  virtual void initIntegration(size_t, size_t) {}

  using Integrand::getLocalIntegral;
  //! \brief Returns a local integral container for the element \a iEl.
  LocalIntegral* getLocalIntegral(size_t, size_t iEl,
                                  bool = false) const override;

  //! \brief Dummy implementation (only boundary integration is relevant).
  bool initElement(const std::vector<int>&,
                   LocalIntegral&) override { return false; }

  //! \brief Dummy implementation (only boundary integration is relevant).
  bool initElement(const std::vector<int>&,
                   const FiniteElement&, const Vec3&, size_t,
                   LocalIntegral&) override { return false; }

  //! \brief Dummy implementation (only boundary integration is relevant).
  bool initElement(const std::vector<int>&,
                   const std::vector<size_t>&,
                   const std::vector<size_t>&,
                   LocalIntegral&) override { return false; }

  //! \brief Dummy implementation (only boundary integration is relevant).
  bool initElement(const std::vector<int>&,
                   const MxFiniteElement&,
                   const std::vector<size_t>&,
                   const std::vector<size_t>&,
                   LocalIntegral&) override { return false; }

  //! \brief Initializes current element for boundary integration.
  bool initElementBou(const std::vector<int>& MNPC,
                      LocalIntegral& elmInt) override;
  //! \brief Initializes current element for boundary integration (mixed).
  bool initElementBou(const std::vector<int>& MNPC,
                      const std::vector<size_t>& elem_sizes,
                      const std::vector<size_t>& basis_sizes,
                      LocalIntegral& elmInt) override;

  //! \brief Returns the number of force components.
  virtual size_t getNoComps() const = 0;

  //! \brief Returns whether this integrand has explicit interior contributions.
  virtual bool hasInteriorTerms() const { return false; }
  //! \brief Returns whether this integrand has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Assigns a parameter value to property functions of the integrand.
  //! \param[in] name Parameter name
  //! \param[in] value Parameter value
  void setParam(const std::string& name, double value) override
  {
    myProblem.setParam(name,value);
  }
  //! \brief Assigns parameter values to property functions of the integrand.
  //! \param[in] name Parameter name
  //! \param[in] value Parameter value
  void setParam(const std::string& name, const Vec3& value) override
  {
    myProblem.setParam(name,value);
  }

  //! \brief Returns current solution mode.
  SIM::SolutionMode getMode(bool) const override { return myProblem.getMode(); }

protected:
  //! \brief Clears out internal buffers.
  void clearBuffer();

  IntegrandBase& myProblem; //!< The problem-specific data
  LintegralVec      eForce; //!< Local integrals used during force integration
  double*          eBuffer; //!< Element force buffer used during integration
};

#endif
