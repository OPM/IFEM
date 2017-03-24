// $Id$
//==============================================================================
//!
//! \file SIMbase.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for NURBS-based FEM simulators.
//!
//==============================================================================

#ifndef _SIM_BASE_H
#define _SIM_BASE_H

#include "SIMadmin.h"
#include "SIMdependency.h"
#include "TimeDomain.h"
#include "Property.h"
#include "Function.h"
#include "MatVec.h"

class IntegrandBase;
class NormBase;
class ForceBase;
class AnaSol;
class SAM;
class AlgEqSystem;
class LinSolParams;
class TimeStep;
class SystemVector;
class Vec4;

//! Property code to integrand map
typedef std::multimap<int,IntegrandBase*> IntegrandMap;


/*!
  \brief Struct for storage of data associated with one mode shape.
*/

struct Mode
{
  int    eigNo;  //!< Eigenvalue identifier
  double eigVal; //!< Eigenvalue associated with this mode
  Vector eigVec; //!< Eigenvector associated with this mode
  // \brief Default constructor.
  Mode() : eigNo(0), eigVal(0.0) {}
};


/*!
  \brief Base class for NURBS-based FEM simulators.
  \details This class incapsulates data and methods need for solving
  partial differential equations using NURBS-based finite elements.
  It only contains the problem-independent data and methods.
  Sub-classes are derived with additional info regarding the problem to solve.
*/

class SIMbase : public SIMadmin, public SIMdependency
{
protected:
  //! \brief The constructor initializes the pointers to dynamic data members.
  SIMbase(IntegrandBase* itg = nullptr);

public:
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~SIMbase();


  // Model input and pre-processing methods
  // ======================================

  //! \brief Creates the computational FEM model from the spline patches.
  //! \param[in] resetNumb If \e 'y', start element and node numbers from zero
  virtual bool createFEMmodel(char resetNumb = 'y') = 0;

  //! \brief Initializes the property containers of the model.
  //! \details Use this method to clear the model before re-reading
  //! the input file in the refinement step of an adaptive simulation.
  virtual void clearProperties();

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \param[in] ignored Indices of patches to ignore in the analysis
  //! \param[in] fixDup Merge duplicated FE nodes on patch interfaces?
  bool preprocess(const std::vector<int>& ignored = std::vector<int>(),
                  bool fixDup = false);

  //! \brief Allocates the system matrices of the FE problem to be solved.
  //! \param[in] mType The matrix format to use
  //! \param[in] nMats Number of system matrices
  //! \param[in] nVec Number of system right-hand-side vectors
  //! \param[in] nScl Number of global scalar quantities
  //! \param[in] withRF Whether nodal reaction forces should be computed or not
  bool initSystem(int mType, size_t nMats = 1, size_t nVec = 1, size_t nScl = 0,
                  bool withRF = false);

  //! \brief Associates a system vector to a system matrix.
  //! \sa AlgEqSystem::setAssociatedVector
  //! \param[in] iMat Index of a coefficient matrix
  //! \param[in] iVec Index of the system vector to associate with the matrix
  bool setAssociatedRHS(size_t iMat, size_t iVec);

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  //! \param[in] resetSol If \e true, the internal solution vectors are cleared
  bool setMode(int mode, bool resetSol = false);

  //! \brief Initializes an integration parameter for the integrand.
  //! \param[in] i Index of the integration parameter to define
  //! \param[in] prm The parameter value to assign
  void setIntegrationPrm(unsigned short int i, double prm);

  //! \brief Defines the spatial numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  //! \param[in] redimBuffers Toggle initialization of internal buffer arrays
  //! \param[in] printQP If \e true, print out total number of quadrature points
  void setQuadratureRule(size_t ng, bool redimBuffers = false,
                         bool printQP = false);

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const;

  //! \brief Returns a pointer to the problem-specific data object.
  const IntegrandBase* getProblem() const { return myProblem; }

  //! \brief Clears the reference to the problem-specific data object.
  //! \details This method is used when the same IntegrandBase object is shared
  //! by several SIMbase objects, to avoid that it is deleted more than once.
  void clearProblem() { myProblem = nullptr; }

  //! \brief Returns the name of this simulator.
  //! \details This method is typically reimplemented in sub-classes that are
  //! parts of a partitioned solution method and are used to identify the basis
  //! for the result fields associated with each simulator in the HDF5 output.
  virtual std::string getName() const { return "SIMbase"; }
  //! \brief Returns whether a mixed formulation is used (used by HDF5 output).
  virtual bool mixedProblem() const { return false; }

  //! \brief Obtain the linear solver parameters.
  const LinSolParams* getSolParams() const { return mySolParams; }

  //! \brief Returns the number of parameter dimensions in the model.
  virtual unsigned short int getNoParamDim() const = 0;
  //! \brief Returns the number of spatial dimensions in the model.
  virtual size_t getNoSpaceDim() const { return nsd; }
  //! \brief Returns the number of primary solution fields.
  //! \param[in] basis Which basis to consider when mixed methods (0 = all)
  size_t getNoFields(int basis = 0) const;
  //! \brief Returns the model size in terms of number of DOFs.
  size_t getNoDOFs() const;
  //! \brief Returns the model size in terms of number of (unique) nodes.
  size_t getNoNodes(bool unique = false, int basis = 0) const;
  //! \brief Returns the model size in terms of number of elements.
  //! \param[in] includeXElms If \e true, include any extra-ordinary elements
  size_t getNoElms(bool includeXElms = false) const;
  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions() const;
  //! \brief Returns the total number of patches in the model.
  int getNoPatches() const { return nGlPatches; }
  //! \brief Returns the number of right-hand-side vectors.
  virtual size_t getNoRHS() const;
  //! \brief Returns the number of bases in the model.
  unsigned char getNoBasis() const;

  //! \brief Returns the type (DOF classification) of the specified global node.
  char getNodeType(int inod) const;
  //! \brief Returns the spatial coordinates of the specified global node.
  Vec4 getNodeCoord(int inod) const;
  //! \brief Returns \e true if all DOFs in the specified global node are fixed.
  bool isFixed(int inod, int dof = 123) const;
  //! \brief Returns the global node number from a process-local node number.
  int getGlobalNode(int node) const;
  //! \brief Returns the process-local node number from a global node number.
  int getLocalNode(int node) const;

  //! \brief Finds the list of global nodes associated with a boundary.
  //! \param[in] pcode Property code identifying the boundary
  //! \param[out] glbNodes Global node numbers on the boundary
  //! \param[out] XYZ Spatial coordinates of the boundary nodes (optional)
  void getBoundaryNodes(int pcode, std::vector<int>& glbNodes,
                        std::vector<Vec3>* XYZ = nullptr) const;

  //! \brief Finds the node that is closest to the given point \b X.
  int findClosestNode(const Vec3&) const;

  //! \brief Initializes time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  bool initDirichlet(double time = 0.0);
  //! \brief Checks for time-dependent in-homogeneous Dirichlet conditions.
  bool hasTimeDependentDirichlet() const;

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  //! \param[in] prevSol Pointer to previous primary solution in DOF-order
  //!
  //! \details If \a prevSol is null, the Dirichlet coefficients are set to
  //! zero (used when doing equilibrium iterations on a fixed time level).
  //! If \a prevSol points to an empty vector, the coefficients are set to the
  //! current (updated) values of the functions defining the Dirichlet condition
  //! (used for the initial time step), otherwise they are set to the difference
  //! between the new values from the Dirichlet functions, and the previous
  //! values stored in the provided \a prevSol vector.
  virtual bool updateDirichlet(double time = 0.0,
                               const Vector* prevSol = nullptr);

  //! \brief Updates problem-dependent state based on the current solution.
  virtual bool updateConfiguration(const Vector&) { return true; }
  //! \brief Updates the nodal rotations for problems with rotational DOFs.
  virtual bool updateRotations(const Vector&, double = 0.0) { return true; }

  //! \brief Updates the grid coordinates.
  //! \param[in] displ The displacement increment to update the grid with
  bool updateGrid(const Vector& displ);
  //! \brief Updates the grid coordinates.
  //! \param[in] field Name of the displacement increment field to update with
  bool updateGrid(const std::string& field);


  // Computational methods
  // =====================

  //! \brief Administers assembly of the linear equation system.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] pSol Previous primary solution vectors in DOF-order
  //! \param[in] newLHSmatrix If \e false, only integrate the RHS vector
  //! \param[in] poorConvg If \e true, the nonlinear driver is converging poorly
  virtual bool assembleSystem(const TimeDomain& time, const Vectors& pSol,
                              bool newLHSmatrix = true, bool poorConvg = false);

  //! \brief Administers assembly of the linear equation system.
  //! \param[in] pSol Primary solution vectors in DOF-order
  //!
  //! \details Use this version for linear/stationary problems only.
  bool assembleSystem(const Vectors& pSol = Vectors())
  { return this->assembleSystem(TimeDomain(),pSol); }

  //! \brief Extracts the assembled load vector for inspection/visualization.
  //! \param[out] loadVec Global load vector in DOF-order
  bool extractLoadVec(Vector& loadVec) const;
  //! \brief Extracts an assembled global scalar quantity.
  double extractScalar(size_t i = 0) const;

  //! \brief Applies the Dirichlet conditions to given vector.
  //! \param[out] glbVec Global vector in DOF-order
  bool applyDirichlet(Vector& glbVec) const;

  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[out] rCond Reciprocal condition number
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] newLHS If \e false, reuse the LHS-matrix from previous call.
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  virtual bool solveSystem(Vector& solution, int printSol, double* rCond,
                           const char* compName = "displacement",
                           bool newLHS = true, size_t idxRHS = 0);

  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[in] compName Solution name to be used in norm output
  bool solveSystem(Vector& solution, int printSol = 0,
                   const char* compName = "displacement")
  {
    return this->solveSystem(solution,printSol,nullptr,compName);
  }

  //! \brief Solves a linear system of equations with multiple right-hand-sides.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[in] compName Solution name to be used in norm output
  bool solveMatrixSystem(Vectors& solution, int printSol = 0,
                         const char* compName = "displacement");

  //! \brief Finds the DOFs showing the worst convergence behavior.
  //! \param[in] x Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[in] nWorst How many bad DOFs to detect
  //! \param[in] eps Only record values larger than this tolerance
  //! \param[in] iteNorm Which norm to consider (1=res, 2=dis, 3=energy)
  //! \param[out] worst Node and local DOF number and values of the worst DOFs
  void getWorstDofs(const Vector& x, const Vector& r,
                    size_t nWorst, double eps, int iteNorm,
                    std::map<std::pair<int,int>,RealArray>& worst) const;

  //! \brief Evaluates some iteration norms for convergence assessment.
  //! \param[in] x Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[out] eNorm Energy norm of solution increment
  //! \param[out] rNorm Residual norm of solution increment
  //! \param[out] dNorm Displacement norm of solution increment
  void iterationNorms(const Vector& x, const Vector& r,
                      double& eNorm, double& rNorm, double& dNorm) const;

  //! \brief Evaluates some norms of the primary solution vector
  //! \param[in] x Global primary solution vector
  //! \param[out] inf Infinity norms in each spatial direction
  //! \param[out] ind Global index of the node corresponding to the inf-value
  //! \param[in] nf Number of components in the primary solution field
  //! \param[in] type Only consider nodes of this DOF type (for mixed methods)
  //! \return L2-norm of the solution vector
  double solutionNorms(const Vector& x, double* inf = nullptr,
                       size_t* ind = nullptr, size_t nf = 0,
                       char type = 'D') const;

  //! \brief Integrates some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations.
  //! \param[in] psol Primary solution vectors
  //! \param[in] ssol Secondary solution vectors
  //! \param[out] gNorm Global norm quantities
  //! \param[out] eNorm Element-wise norm quantities
  //!
  //! \details If an analytical solution is provided, norms of the exact
  //! error in the solution are computed as well. If projected secondary
  //! solutions are provided (i.e., \a ssol is not empty), norms of the
  //! difference between these solutions and the directly evaluated secondary
  //! solution are computed as well.
  bool solutionNorms(const TimeDomain& time,
                     const Vectors& psol, const Vectors& ssol,
                     Vectors& gNorm, Matrix* eNorm = nullptr);
  //! \brief Integrates some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations.
  //! \param[in] psol Primary solution vectors
  //! \param[out] gNorm Global norm quantities
  //! \param[out] eNorm Element-wise norm quantities
  //!
  //! \details Use this version if no projected solutions are needed/available.
  bool solutionNorms(const TimeDomain& time, const Vectors& psol,
                     Vectors& gNorm, Matrix* eNorm = nullptr)
  { return this->solutionNorms(time,psol,Vectors(),gNorm,eNorm); }
  //! \brief Integrates some solution norm quantities.
  //! \param[in] psol Primary solution vectors
  //! \param[in] ssol Secondary solution vectors
  //! \param[out] eNorm Element-wise norm quantities
  //! \param[out] gNorm Global norm quantities
  //!
  //! \details Use this version for linear/stationary problems only.
  bool solutionNorms(const Vectors& psol, const Vectors& ssol,
                     Matrix& eNorm, Vectors& gNorm)
  { return this->solutionNorms(TimeDomain(),psol,ssol,gNorm,&eNorm); }
  //! \brief Integrates some solution norm quantities.
  //! \param[in] psol Primary solution vectors
  //! \param[out] eNorm Element-wise norm quantities
  //! \param[out] gNorm Global norm quantities
  //!
  //! \details Use this version for linear/stationary problems,
  //! and when no projected solutions are needed/available.
  bool solutionNorms(const Vectors& psol, Matrix& eNorm, Vectors& gNorm)
  { return this->solutionNorms(TimeDomain(),psol,Vectors(),gNorm,&eNorm); }

  //! \brief Prints integrated solution norms to the log stream.
  //! \param[in] norms The norm values
  //! \param[in] w Total number of characters in the norm labels
  virtual void printNorms(const Vectors& norms, size_t w = 36) const;

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solution The solution vector
  //! \param[in] printSol Print solution only if size is less than this value
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual void printSolutionSummary(const Vector& solution, int printSol = 0,
                                    const char* compName = nullptr,
                                    std::streamsize outPrec = 0);

  //! \brief Computes the total reaction forces in the model.
  //! \param[out] RF Reaction force in each spatial direction + energy
  //! \param[in] psol Primary solution vector
  bool getCurrentReactions(RealArray& RF, const Vector& psol) const;
  //! \brief Computes the total reaction forces associated with a boundary.
  //! \param[out] RF Reaction force in each spatial direction
  //! \param[in] pcode Property code identifying the boundary
  bool getCurrentReactions(RealArray& RF, int pcode) const;
  //! \brief Returns current reaction force vector.
  const Vector* getReactionForces() const;

  //! \brief Performs a generalized eigenvalue analysis of the assembled system.
  //! \param[in] iop Which eigensolver method to use
  //! \param[in] nev Number of eigenvalues/vector (see ARPack documentation)
  //! \param[in] ncv Number of Arnoldi vectors (see ARPack documentation)
  //! \param[in] shift Eigenvalue shift
  //! \param[out] solution Computed eigenvalues and associated eigenvectors
  //! \param[in] iA Index of system matrix \b A in \a myEqSys->A
  //! \param[in] iB Index of system matrix \b B in \a myEqSys->A
  bool systemModes(std::vector<Mode>& solution,
                   int nev, int ncv, int iop, double shift,
                   size_t iA = 0, size_t iB = 1);
  //! \brief Performs a generalized eigenvalue analysis of the assembled system.
  //! \param[out] solution Computed eigenvalues and associated eigenvectors
  //! \param[in] iA Index of system matrix \b A in \a myEqSys->A
  //! \param[in] iB Index of system matrix \b B in \a myEqSys->A
  bool systemModes(std::vector<Mode>& solution, size_t iA = 0, size_t iB = 1)
  {
    return this->systemModes(solution,opt.nev,opt.ncv,opt.eig,opt.shift,iA,iB);
  }


  // Post-processing methods
  // =======================

  //! \brief Projects the secondary solution associated with a primary solution.
  //! \param[out] ssol Control point values of the secondary solution
  //! \param[in] psol Control point values of the primary solution
  //! \param[in] pMethod Projection method to use
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //!
  //! \details The secondary solution, defined through the Integrand object,
  //! corresponding to the primary solution \a psol is projected onto the
  //! spline basis to obtain the control point values of the secondary solution.
  bool project(Matrix& ssol, const Vector& psol,
               SIMoptions::ProjectionMethod pMethod = SIMoptions::GLOBAL,
               const TimeDomain& time = TimeDomain()) const;

  //! \brief Projects a scalar function onto the specified basis.
  //! \param[out] values Resulting control point values
  //! \param[in] f The function to evaluate
  //! \param[in] basis Which basis to consider
  //! \param[in] iField Field component offset in values vector
  //! \param[in] nFields Number of field components in values vector
  //! \param[in] time Current time
  bool project(Vector& values, const RealFunc* f,
               int basis = 1, int iField = 0, int nFields = 1,
               double time = 0.0) const;

  //! \brief Evaluates the secondary solution field for specified patch.
  //! \param[out] field Control point values of the secondary solution field
  //! \param[in] pindx Local patch index to evaluate solution field for
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the Integrand object.
  //! The solution is evaluated at the Greville points and then projected onto
  //! the spline basis to obtain the control point values.
  bool evalSecondarySolution(Matrix& field, int pindx) const;

  //! \brief Returns whether an analytical solution is available or not.
  bool haveAnaSol() const { return mySol ? true : false; }

  //! \brief Returns a pointer to a norm integrand object for this simulator.
  //! \note The object is allocated dynamically and has therefore to be
  //! manually deleted before the variable receiving the pointer value goes
  //! out of scope.
  NormBase* getNormIntegrand() const;
  //! \brief Returns a pointer to a force integrand object for this simulator.
  //! \note The object is allocated dynamically and has therefore to be
  //! manually deleted before the variable receiving the pointer value goes
  //! out of scope.
  ForceBase* getBoundaryForceIntegrand(const Vec3* X0 = nullptr) const;
  //! \brief Returns a pointer to a force integrand object for this simulator.
  //! \note The object is allocated dynamically and has therefore to be
  //! manually deleted before the variable receiving the pointer value goes
  //! out of scope.
  ForceBase* getNodalForceIntegrand() const;

  //! \brief Returns the SAM object for this SIM
  const SAM* getSAM() const { return mySam; }

protected:
  //! \brief Returns a vector function associated with given patch and property.
  //! \param[in] patch 1-based index of the patch to retrieve property for
  //! \param[in] ptype The property type associated with the vector function
  VecFunc* getVecFunc(size_t patch, Property::Type ptype) const;

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homegeneous Dirichlet condition property code
  //! \param ngnod Total number of global nodes in the model (might be updated)
  //! \param[in] basis Which basis to apply the constraint to (mixed methods)
  virtual bool addConstraint(int patch, int lndx, int ldim, int dirs, int code,
                             int& ngnod, char basis = 1) = 0;

  //! \brief Preprocessing performed before the FEM model generation.
  virtual void preprocessA() {}
  //! \brief Specialized preprocessing performed before assembly initialization.
  virtual bool preprocessBeforeAsmInit(int&) { return true; }
  //! \brief Preprocessing performed after the system assembly initialization.
  virtual bool preprocessB() { return true; }
  //! \brief Preprocesses the result sampling points.
  virtual void preprocessResultPoints() = 0;

  //! \brief Extracts all local solution vector(s) for a specified patch.
  //! \param[in] problem Integrand to receive the patch-level solution vectors
  //! \param[in] sol Global primary solution vectors in DOF-order
  //! \param[in] pindx Local patch index to extract solution vectors for
  //!
  //! \details This method is typically invoked before the \a integrate method
  //! on the the specified path, in order to extract all patch-level vector
  //! quantities needed by the Integrand. This also includes any dependent
  //! vectors from other simulator classes that have been registered.
  //! All patch-level vectors are stored within the provided integrand.
  virtual bool extractPatchSolution(IntegrandBase* problem,
                                    const Vectors& sol, size_t pindx) const;

public:
  //! \brief Extracts all local solution vector(s) for a specified patch.
  //! \param[in] sol Global primary solution vectors in DOF-order
  //! \param[in] pindx Local patch index to extract solution vectors for
  bool extractPatchSolution(const Vectors& sol, size_t pindx) const
  { return this->extractPatchSolution(myProblem,sol,pindx); }

  //! \brief Extracts a local solution vector for a specified patch.
  //! \param[in] sol Global primary solution vector in DOF-order
  //! \param[out] vec Local solution vector associated with specified patch
  //! \param[in] pindx Local patch index to extract solution vector for
  //! \param[in] nndof Number of DOFs per node (optional)
  //! \param[in] basis Basis to extract for (optional)
  //! \return Total number of DOFs in the patch (first basis only if mixed)
  size_t extractPatchSolution(const Vector& sol, Vector& vec, int pindx,
                              unsigned char nndof = 0,
                              unsigned char basis = 0) const;

  //! \brief Injects a patch-wise solution vector into the global vector.
  //! \param sol Global primary solution vector in DOF-order
  //! \param[in] vec Local solution vector associated with specified patch
  //! \param[in] pindx Local patch index to inject solution vector for
  //! \param[in] nndof Number of DOFs per node (optional)
  //! \param[in] basis Basis to inject for (optional)
  bool injectPatchSolution(Vector& sol, const Vector& vec, int pindx,
                           unsigned char nndof = 0,
                           unsigned char basis = 0) const;

  //! \brief Extracts element results for a specified patch.
  //! \param[in] glbRes Global element result array
  //! \param[out] elRes Patch-level element result array
  //! \param[in] pindx Local patch index to extract element results for
  bool extractPatchElmRes(const Matrix& glbRes, Matrix& elRes, int pindx) const;

  //! \brief Returns the local patch index for the given global patch number.
  //! \details For serial applications this is an identity mapping only, whereas
  //! for parallel applications the local (1-based) patch index on the current
  //! processor is returned. If \a patchNo is out of range, -1 is returned.
  //! If \a patchNo is not on current processor, 0 is returned.
  int getLocalPatchIndex(int patchNo) const;

  //! \brief Returns a const reference to our FEM model.
  const PatchVec& getFEModel() const { return myModel; }
  //! \brief Returns a pointer to a specified patch of our FEM model.
  //! \param[in] idx 1-based patch index
  //! \param[in] glbIndex If \e true, the patch index is assumed to be global
  //! for the whole model, otherwise it is assumed local within current process.
  //! For serial applications this option has no effect.
  ASMbase* getPatch(int idx, bool glbIndex = false) const;

  //! \brief Initializes material properties for the given patch.
  bool setPatchMaterial(size_t patch);

  //! \brief Returns a const reference to our global-to-local node mapping.
  const std::map<int,int>& getGlob2LocMap() const { return myGlb2Loc; }

  //! \brief Returns the beginning of the property array.
  PropertyVec::const_iterator begin_prop() const { return myProps.begin(); }
  //! \brief Returns the end of the property array.
  PropertyVec::const_iterator end_prop() const { return myProps.end(); }

  //! \brief Returns current system light-hand-side vector.
  SystemVector* getRHSvector(size_t idx = 0, bool copy = false) const;
  //! \brief Adds a system vector to the given right-hand-side vector.
  void addToRHSvector(size_t idx, const SystemVector& vec, double scale = 1.0);

  //! \brief Returns a scalar function associated with \a code.
  RealFunc* getSclFunc(int code) const;

protected:
  //! \brief Initializes material properties for integration of interior terms.
  virtual bool initMaterial(size_t) { return true; }
  //! \brief Initializes the body load properties for current patch.
  virtual bool initBodyLoad(size_t) { return true; }
  //! \brief Initializes for integration of Neumann terms for a given property.
  virtual bool initNeumann(size_t) { return true; }

  //! \brief Assembles problem-dependent discrete terms, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase*,
                                     const TimeDomain&) { return true; }

  //! \brief Computes (possibly problem-dependent) external energy contribution.
  virtual double externalEnergy(const Vectors& psol) const;

  //! \brief Generates element groups for multi-threading of boundary integrals.
  //! \param[in] p Property object identifying a patch boundary
  //! \param[in] silence If \e true, suppress threading group outprint
  void generateThreadGroups(const Property& p, bool silence = false);

  //! \brief Adds a MADOF with an extraordinary number of DOFs on a given basis.
  //! \param[in] basis The basis to specify number of DOFs for
  //! \param[in] nndof Number of nodal DOFs on the given basis
  bool addMADOF(unsigned char basis, unsigned char nndof);

private:
  //! \brief Returns an extraordinary MADOF array.
  //! \param[in] basis The basis to specify number of DOFs for
  //! \param[in] nndof Number of nodal DOFs on the given basis
  const std::vector<int>& getMADOF(unsigned char basis,
                                   unsigned char nndof) const;

public:
  static bool ignoreDirichlet; //!< Set to \e true for free vibration analysis
  static bool preserveNOrder;  //!< Set to \e true to preserve node ordering

  typedef std::vector<unsigned char> CharVec; //!< Convenience declaration

protected:
  //! \brief Scalar field container
  typedef std::map<int,RealFunc*>     SclFuncMap;
  //! \brief Vector field container
  typedef std::map<int,VecFunc*>      VecFuncMap;
  //! \brief Traction field container
  typedef std::map<int,TractionFunc*> TracFuncMap;

  // Model attributes
  bool           isRefined; //!< Indicates if the model is adaptively refined
  bool           lagMTOK;   //!< Indicates that global multipliers is okay with multithreading (app specific).
  unsigned char  nsd;       //!< Number of spatial dimensions
  PatchVec       myModel;   //!< The actual NURBS/spline model
  PropertyVec    myProps;   //!< Physical property mapping
  SclFuncMap     myScalars; //!< Scalar property fields
  VecFuncMap     myVectors; //!< Vector property fields
  TracFuncMap    myTracs;   //!< Traction property fields
  IntegrandBase* myProblem; //!< The main integrand of this simulator
  IntegrandMap   myInts;    //!< Set of all integrands involved
  AnaSol*        mySol;     //!< Analytical/Exact solution

  //! \brief A struct with data for system matrix/vector dumps.
  struct DumpData
  {
    std::string fname;  //!< File name
    char        format; //!< File format flag
    int         step;   //!< Dump step identifier
    int         count;  //!< Internal step counter, dump only when step==count
    //! \brief Default constructor.
    DumpData() : format('P'), step(1), count(0) {}
    //! \brief Checks if the matrix or vector should be dumped now.
    bool doDump() { return !fname.empty() && ++count == step; }
  };

  // Post-processing attributes
  std::vector<DumpData> lhsDump; //!< Coefficient matrix dump specifications
  std::vector<DumpData> rhsDump; //!< Right-hand-side vector dump specifications
  std::vector<DumpData> solDump; //!< Solution vector dump specifications

  // Parallel computing attributes
  int               nGlPatches; //!< Number of global patches
  std::vector<int>  myPatches;  //!< Global patch numbers for current processor
  std::map<int,int> myGlb2Loc;  //!< Global-to-local node number mapping
  const std::map<int,int>* g2l; //!< Pointer to global-to-local node mapping

  // Equation solver attributes
  AlgEqSystem*  myEqSys;     //!< The actual linear equation system
  SAM*          mySam;       //!< Auxiliary data for FE assembly management
  LinSolParams* mySolParams; //!< Input parameters for PETSc

private:
  size_t nIntGP; //!< Number of interior integration points in the whole model
  size_t nBouGP; //!< Number of boundary integration points in the whole model

  //! Additional MADOF arrays for mixed problems (extraordinary DOF counts)
  std::map<int, std::vector<int> > mixedMADOFs;
};

#endif
