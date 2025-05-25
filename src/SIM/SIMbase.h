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
#include "MatVec.h"
#include <set>

class IntegrandBase;
class NormBase;
class ForceBase;
class AnaSol;
class SAM;
class AlgEqSystem;
class LinSolParams;
class SystemMatrix;
class SystemVector;
class FunctionBase;
class RealFunc;
class VecFunc;
class TractionFunc;
class ScalarFunc;
class Vec4;
class Vec3;


/*!
  \brief Struct for storage of data associated with one mode shape.
*/

struct Mode
{
  int    eigNo;   //!< Eigenvalue identifier
  double eigVal;  //!< Eigenvalue associated with this mode
  double damping; //!< Modal damping coefficient
  Vector eigVec;  //!< Eigenvector associated with this mode
  Vector eqnVec;  //!< Eigenvector associated with this mode in equation order

  //! \brief Default constructor.
  Mode() : eigNo(0), eigVal(0.0), damping(0.0) {}
  //! \brief Orthonormalize the eigenvector w.r.t. the given matrix.
  bool orthonormalize(const SystemMatrix& mat);
  //! \brief Compute modal damping based on the given damping matrix.
  bool computeDamping(const SystemMatrix& mat);
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
  explicit SIMbase(IntegrandBase* itg);

public:
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~SIMbase();


  // Model input and pre-processing methods
  // ======================================

  //! \brief Reads model data from the specified input file \a *fileName.
  bool readModel(const char* fileName);

  //! \brief Creates the computational FEM model from the spline patches.
  virtual bool createFEMmodel(char resetNumb) = 0;

  //! \brief Initializes the property containers of the model.
  //! \details Use this method to clear the model before re-reading
  //! the input file in the refinement step of an adaptive simulation.
  virtual void clearProperties();

  //! \brief Interface for app-specific single-step simulation initialisation.
  virtual void initForSingleStep() {}
  //! \brief Interface for app-specific multi-step simulation initialisation.
  virtual void initForMultiStep() {}

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \param[in] ignored Indices of patches to ignore in the analysis
  //! \param[in] fixDup Merge duplicated FE nodes on patch interfaces?
  //! \return Total number of unique FE nodes in the model
  virtual bool preprocess(const std::vector<int>& ignored = std::vector<int>(),
                          bool fixDup = false);

  //! \brief Merges the global equation system of \a that simulator with this.
  //! \param that The simulator whose equation system is to be merged
  //! \param[in] old2new Global node number mapping
  //! \param[in] poff Global patch index offset
  virtual bool merge(SIMbase* that,
                     const std::map<int,int>* old2new = nullptr, int poff = 0);

  //! \brief Allocates the system matrices of the FE problem to be solved.
  //! \param[in] mType The matrix format to use
  //! \param[in] nMats Number of system matrices
  //! \param[in] nVec Number of system right-hand-side vectors
  //! \param[in] nScl Number of global scalar quantities
  //! \param[in] withRF Whether nodal reaction forces should be computed or not
  bool initSystem(LinAlg::MatrixType mType,
                  size_t nMats = 1, size_t nVec = 1, size_t nScl = 0,
                  bool withRF = false);

  //! \brief Lets this simulator share equation system with \a that simulator.
  bool initSystem(const SIMbase* that);

  //! \brief Initializes left-hand-side element matrix buffers for integrand.
  void initLHSbuffers();

  //! \brief Associates a system vector to a system matrix.
  //! \sa AlgEqSystem::setAssociatedVector
  //! \param[in] iMat Index of a coefficient matrix
  //! \param[in] iVec Index of the system vector to associate with the matrix
  bool setAssociatedRHS(size_t iMat, size_t iVec);

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  //! \param[in] needIntegr If \e false, silently ignore non-existing integrand
  //! \param[in] resetSol If \e true, the internal solution vectors are cleared
  bool setMode(int mode, bool needIntegr = true, bool resetSol = false);

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
  virtual bool printProblem() const;

  //! \brief Returns a pointer to the problem-specific data object.
  const IntegrandBase* getProblem() const { return myProblem; }

  //! \brief Returns interface checker type for model.
  virtual ASM::InterfaceChecker* getInterfaceChecker(size_t) const
  { return nullptr; }

  //! \brief Clears the reference to the problem-specific data object.
  //! \details This method is used when the same IntegrandBase object is shared
  //! by several SIMbase objects, to avoid that it is deleted more than once.
  void clearProblem() { myProblem = nullptr; }

  //! \brief Returns the name of this simulator.
  //! \details This method is typically overridden by sub-classes that are
  //! parts of a partitioned solution method and are used to identify the basis
  //! for the result fields associated with each simulator in the HDF5 output.
  virtual std::string getName() const { return "SIMbase"; }
  //! \brief Returns whether a mixed formulation is used (used by HDF5 output).
  virtual bool mixedProblem() const { return false; }

  //! \brief Returns the linear equation solver parameters (for PETSc).
  const LinSolParams* getSolParams() const { return mySolParams; }

  //! \brief Returns the number of parameter dimensions in the model.
  virtual unsigned short int getNoParamDim() const = 0;
  //! \brief Returns the number of spatial dimensions in the model.
  virtual size_t getNoSpaceDim() const { return nsd; }
  //! \brief Returns the number of primary solution fields.
  //! \param[in] basis Which basis to consider when mixed methods (0 = all)
  size_t getNoFields(int basis = 0) const;
  //! \brief Returns the model size in terms of number of DOFs.
  size_t getNoDOFs(bool subSim = false) const;
  //! \brief Returns the model size in terms of number of unique nodes.
  //! \param[in] basis Which basis to return the number of nodes for (0 = all)
  size_t getNoNodes(int basis = 0) const;
  //! \brief Returns the model size in terms of number of elements.
  //! \param[in] includeXelms If \e true, include any extra-ordinary elements
  //! \param[in] includeZelms If \e true, count all regular elements,
  //! including zero-volume elements due to multiple knots, etc.
  size_t getNoElms(bool includeXelms = false, bool includeZelms = false) const;
  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions(bool allocated = true) const;
  //! \brief Returns the total number of patches in the model.
  int getNoPatches() const { return nGlPatches; }
  //! \brief Returns the number of unknowns in the linear equation system.
  size_t getNoEquations() const;
  //! \brief Returns the number of constraint equations in the model.
  size_t getNoConstraints() const;
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
  //! \brief Finds the Matrix of Nodal Point Correspondance for element \a iel.
  bool getElmNodes(std::vector<int>& mnpc, int iel) const;
  //! \brief Obtain element-element connectivities
  virtual std::vector<std::vector<int>> getElmConnectivities() const = 0;

  //! \brief Finds the list of global nodes associated with a boundary.
  //! \param[in] pcode Property code identifying the boundary
  //! \param[out] glbNodes Global node numbers on the boundary
  //! \param[out] XYZ Spatial coordinates of the boundary nodes (optional)
  void getBoundaryNodes(int pcode, std::vector<int>& glbNodes,
                        std::vector<Vec3>* XYZ = nullptr) const;

  //! \brief Finds the node that is closest to the given point \b X.
  int findClosestNode(const Vec3&) const;

  //! \brief Returns a predefined node set.
  std::vector<int> getNodeSet(const std::string& setName) const;

  //! \brief Initializes time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  bool initDirichlet(double time = 0.0);
  //! \brief Checks for time-dependent in-homogeneous Dirichlet conditions.
  bool hasTimeDependentDirichlet() const;

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  //! \param[in] prevSol Pointer to previous primary solution in DOF-order
  virtual bool updateDirichlet(double time = 0.0,
                               const Vector* prevSol = nullptr);

  //! \brief Updates problem-dependent state based on the current solution.
  virtual bool updateConfiguration(const Vector&) { return true; }
  //! \brief Updates the nodal rotations for problems with rotational DOFs.
  virtual bool updateRotations(const RealArray&, double = 0.0) { return true; }

  //! \brief Updates the grid coordinates.
  //! \param[in] displ The displacement increment to update the grid with
  bool updateGrid(const RealArray& displ);
  //! \brief Updates the grid coordinates.
  //! \param[in] field Name of the displacement increment field to update with
  bool updateGrid(const std::string& field);

  //! \brief Sets the refinement status (for restart of adaptive simulations).
  //! \param[in] nref Number of refinement levels so far
  void setRefined(int nref) { isRefined = nref; }
  //! \brief Returns current refinement status.
  int getRefined() const { return isRefined; }

  //! \brief Returns \e true if an element activation function is specified.
  bool hasElementActivator() const;

  //! \brief Modifies the current solution vector when activating elements.
  //! \param solution Current primary solution vector
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  void updateForNewElements(Vector& solution, const TimeDomain& time) const;


  // Computational methods
  // =====================

  //! \brief Administers assembly of the linear equation system.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] prevSol Previous primary solution vectors in DOF-order
  //! \param[in] newLHSmatrix If \e false, only integrate the RHS vector
  //! \param[in] poorConvg If \e true, the nonlinear driver is converging poorly
  virtual bool assembleSystem(const TimeDomain& time, const Vectors& prevSol,
                              bool newLHSmatrix = true, bool poorConvg = false);

  //! \brief Administers assembly of the linear equation system.
  //! \param[in] t0 Time for evaluation of time-dependent property functions
  //! \param[in] pSol Primary solution vectors in DOF-order
  //!
  //! \details Use this version for linear/stationary problems only.
  bool assembleSystem(double t0 = 0.0, const Vectors& pSol = Vectors())
  { return this->assembleSystem(TimeDomain(t0),pSol); }

  //! \brief Extracts the assembled load vector for inspection/visualization.
  //! \param[out] loadVec Global load vector in DOF-order
  //! \param[in] idx Index to the system vector to extract
  //! \param[in] hd Header for outprint of resultant
  bool extractLoadVec(Vector& loadVec, size_t idx = 0,
                      const char* hd = nullptr) const;
  //! \brief Extracts the assembled global scalar quantities.
  bool extractScalars(RealArray& values) const;
  //! \brief Extracts an assembled global scalar quantity.
  double extractScalar(size_t idx = 0) const;

  //! \brief Applies the Dirichlet conditions to given vector.
  //! \param[out] glbVec Global vector in DOF-order
  bool applyDirichlet(Vector& glbVec) const;

  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  //! \param[out] rCond Reciprocal condition number
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[in] dumpEqSys If \e true, activate dump of equation system to file
  bool solveEqSystem(Vector& solution, size_t idxRHS, double* rCond,
                     int printSol = 0, bool dumpEqSys = false,
                     const char* compName = "displacement");

  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[out] rCond Reciprocal condition number
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  virtual bool solveSystem(Vector& solution, int printSol, double* rCond,
                           const char* compName = "displacement",
                           size_t idxRHS = 0)
  {
    return this->solveEqSystem(solution,idxRHS,rCond,printSol,true,compName);
  }

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
  //! \param[out] solution Global primary solution vectors
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[in] cmpName Solution name to be used in norm output
  bool solveSystem(Vectors& solution, int printSol = 0,
                   const char* cmpName = "displacement");

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
  virtual void iterationNorms(const Vector& x, const Vector& r, double& eNorm,
                              double& rNorm, double& dNorm) const;

  //! \brief Evaluates some norms of the primary solution vector.
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
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] psol Primary solution vectors
  //! \param[in] ssol Secondary solution vectors
  //! \param[out] gNorm Global norm quantities
  //! \param[out] eNorm Element-wise norm quantities
  //! \param[in] name Name of solution being the source of calculation
  bool solutionNorms(const TimeDomain& time,
                     const Vectors& psol, const Vectors& ssol,
                     Vectors& gNorm, Matrix* eNorm = nullptr,
                     const char* name = nullptr);
  //! \brief Integrates some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
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
  //! \param[in] name Name of solution being the source of calculation
  //!
  //! \details Use this version for linear/stationary problems only.
  bool solutionNorms(const Vector& psol, const Vectors& ssol,
                     Matrix& eNorm, Vectors& gNorm,
                     const char* name = nullptr)
  {
    return this->solutionNorms(TimeDomain(),Vectors(1,psol),ssol,
                               gNorm,&eNorm,name);
  }
  //! \brief Integrates some solution norm quantities.
  //! \param[in] psol Primary solution vectors
  //! \param[out] eNorm Element-wise norm quantities
  //! \param[out] gNorm Global norm quantities
  //!
  //! \details Use this version for linear/stationary problems,
  //! and when no projected solutions are needed/available.
  bool solutionNorms(const Vector& psol, Matrix& eNorm, Vectors& gNorm)
  {
    return this->solutionNorms(TimeDomain(),Vectors(1,psol),Vectors(),
                               gNorm,&eNorm);
  }
  //! \brief Integrates some solution norm quantities.
  //! \param[in] psol Primary solution vectors
  //! \param[in] ssol Secondary solution vectors
  //! \param[out] gNorm Global norm quantities
  //! \param[in] name Name of solution being the source of calculation
  //!
  //! \details Use this version for linear/stationary problems,
  //! and when the element-wise norms are not needed.
  bool solutionNorms(const Vector& psol, const Vectors& ssol, Vectors& gNorm,
                     const char* name = nullptr)
  {
    return this->solutionNorms(TimeDomain(),Vectors(1,psol),ssol,gNorm,
                               nullptr,name);
  }

  //! \brief Prints out load/time step identification.
  //! \param[in] istep Load- or time step counter
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //!
  //! \details This method is used by multi-step simulators to print out
  //! a heading when starting a new load/time increment. Override this method
  //! if your simulator have some additional data to be printed.
  virtual void printStep(int istep, const TimeDomain& time) const;

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
  //! \param[in] pcode Property code identifying the boundary for which the
  //! reaction forces are computed (0 means the entire model, or all boundaries)
  bool getCurrentReactions(RealArray& RF, const Vector& psol,
                           int pcode = 0) const;
  //! \brief Checks for total reaction forces associated with a boundary.
  //! \param[in] pcode Property code identifying the boundary (0 = all)
  bool haveReactions(int pcode = 0) const;
  //! \brief Returns current reaction force container.
  virtual const RealArray* getReactionForces() const;

  //! \brief Performs a generalized eigenvalue analysis of the assembled system.
  //! \param[in] iop Which eigensolver method to use
  //! \param[in] nev Number of eigenvalues/vector (see ARPack documentation)
  //! \param[in] ncv Number of Arnoldi vectors (see ARPack documentation)
  //! \param[in] shift Eigenvalue shift
  //! \param[out] solution Computed eigenvalues and associated eigenvectors
  //! \param[in] iA Index of system matrix \b A in \ref myEqSys
  //! \param[in] iB Index of system matrix \b B in \ref myEqSys
  bool systemModes(std::vector<Mode>& solution,
                   int nev, int ncv, int iop, double shift,
                   size_t iA = 0, size_t iB = 1);
  //! \brief Performs a generalized eigenvalue analysis of the assembled system.
  //! \param[out] solution Computed eigenvalues and associated eigenvectors
  //! \param[in] iA Index of system matrix \b A in \ref myEqSys
  //! \param[in] iB Index of system matrix \b B in \ref myEqSys
  bool systemModes(std::vector<Mode>& solution, size_t iA = 0, size_t iB = 1)
  {
    return this->systemModes(solution,opt.nev,opt.ncv,opt.eig,opt.shift,iA,iB);
  }

  //! \brief Returns whether reaction forces are to be computed or not.
  virtual bool haveBoundaryReactions(bool = false) const { return false; }

  //! \brief Assembles reaction and interface forces for specified boundaries.
  //! \param[in] solution Current primary solution vector
  //! \param[in] t0 Current time (for time-dependent loads)
  //! \param[out] R Nodal reaction force container
  //! \param[out] S Nodal interface force container
  bool assembleForces(const Vector& solution, double t0,
                      RealArray* R, Vector* S = nullptr);


  // Post-processing methods
  // =======================

  //! \brief Projects the secondary solution associated with a primary solution.
  //! \param[out] ssol Control point values of the secondary solution
  //! \param[in] psol Control point values of the primary solution
  //! \param[in] method Projection method to use
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //!
  //! \details The secondary solution, defined through the Integrand object,
  //! corresponding to the primary solution \a psol is projected onto the
  //! spline basis to obtain the control point values of the secondary solution.
  virtual bool project(Matrix& ssol, const Vector& psol,
                       SIMoptions::ProjectionMethod method = SIMoptions::GLOBAL,
                       const TimeDomain& time = TimeDomain()) const;
  //! \brief Projects the secondary solution associated with a primary solution.
  //! \param[out] ssol Vector of control point values of the secondary solution
  //! \param[in] psol Vector of control point values of the primary solution
  //! \param[in] method Projection method to use
  //! \param[in] iComp One-based index of the component to return (0 = all)
  //!
  //! \details Convenience overload, for stationary problems only.
  bool project(Vector& ssol, const Vector& psol,
               SIMoptions::ProjectionMethod method = SIMoptions::GLOBAL,
               size_t iComp = 0) const;

  //! \brief Projects the analytical secondary solution, if any.
  //! \param[out] ssol Vector of control point values of the secondary solution
  //! \param[in] method Projection method to use
  bool projectAnaSol(Vector& ssol, SIMoptions::ProjectionMethod method) const;

  //! \brief Projects a function onto the specified basis.
  //! \param[out] values Resulting control point values
  //! \param[in] f The function to evaluate
  //! \param[in] basis Which basis to consider
  //! \param[in] iField Field component offset in values vector
  //! \param[in] nFields Number of field components in values vector
  //! \param[in] method Projection method to use
  //! \param[in] time Current time
  bool project(RealArray& values, const FunctionBase* f,
               int basis = 1, int iField = 0, int nFields = 1,
               SIMoptions::ProjectionMethod method = SIMoptions::GLOBAL,
               double time = 0.0) const;

  //! \brief Evaluates the secondary solution field for specified patch.
  //! \param[out] field Control point values of the secondary solution field
  //! \param[in] pindx Local patch index to evaluate solution field for
  bool evalSecondarySolution(Matrix& field, int pindx) const;

  //! \brief Returns whether projections must be handled through fields or not.
  virtual bool fieldProjections() const;

  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const { return mySol ? true : false; }
  //! \brief Returns whether a dual solution is available or not.
  virtual bool haveDualSol() const { return dualField ? true : false; }

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

  //! \brief Returns a const pointer to the SAM object of this simulator.
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
  //! \param[in] code In-homogeneous Dirichlet condition property code
  //! \param ngnod Total number of global nodes in the model (might be updated)
  //! \param[in] basis Which basis to apply the constraint to (mixed methods)
  //! \param[in] ovrD If \e true, override conflicting Dirichlet conditions
  virtual bool addConstraint(int patch, int lndx, int ldim, int dirs, int code,
                             int& ngnod, char basis = 1, bool ovrD = false);

  //! \brief Preprocessing performed before the FEM model generation.
  virtual void preprocessA() {}
  //! \brief Specialized preprocessing performed before assembly initialization.
  virtual bool preprocessBeforeAsmInit(int&) { return true; }
  //! \brief Preprocessing performed after the system assembly initialization.
  virtual bool preprocessB() { return true; }
  //! \brief Preprocesses the result sampling points.
  virtual void preprocessResultPoints() = 0;

  //! \brief Renumbers the global node numbers after resolving patch topology.
  //! \param[in] renumMNPC If \e true, also update element connectivity tables
  //! \return Total number of unique nodes in the model, negative on error
  int renumberNodes(bool renumMNPC = false);
  //! \brief Interface for renumbering of app-specific node number tables.
  virtual bool renumberNodes(const std::map<int,int>&) { return true; }

  //! \brief Extracts all local solution vector(s) for a specified patch.
  //! \param[in] problem The integrand to receive patch-level solution vectors
  //! \param[in] sol Global primary solution vectors in DOF-order
  //! \param[in] pindx Local patch index to extract solution vectors for
  //!
  //! \details This method is typically invoked before ASMbase::integrate()
  //! on the specified patch, in order to extract all patch-level vector
  //! quantities needed by the Integrand. This also includes any dependent
  //! vectors from other simulator classes that have been registered.
  //! All patch-level vectors are stored within the provided integrand.
  virtual bool extractPatchSolution(IntegrandBase* problem,
                                    const Vectors& sol, size_t pindx) const;

public:
  using SIMdependency::registerDependency;
  //! \brief Registers a dependency on a field from another %SIM object.
  //! \param[in] name Name of field we depend on
  //! \param[in] sim The %SIM object holding the field we depend on
  //! \param[in] nvc Number of components in field
  //! \param[in] basis The basis the dependent field should be defined on
  void registerDependency(const std::string& name, SIMdependency* sim,
                          short int nvc = 1, unsigned char basis = 1);

  //! \brief Extracts all local solution vector(s) for a specified patch.
  //! \param[in] sol Global primary solution vectors in DOF-order
  //! \param[in] pindx Local patch index to extract solution vectors for
  bool extractPatchSolution(const Vectors& sol, size_t pindx) const
  { return this->extractPatchSolution(myProblem,sol,pindx); }

  //! \brief Extracts a local solution vector for a specified patch.
  //! \param[in] sol Global primary solution vector in DOF-order
  //! \param[out] vec Local solution vector associated with specified patch
  //! \param[in] pch The patch to extract solution vector for
  //! \param[in] nndof Number of DOFs per node (optional)
  //! \param[in] basis Basis to extract for (optional)
  //! \return Total number of DOFs in the patch (first basis only if mixed)
  size_t extractPatchSolution(const RealArray& sol, RealArray& vec,
                              const ASMbase* pch,
                              unsigned char nndof = 0,
                              unsigned char basis = 0) const;

  //! \brief Injects a patch-wise solution vector into the global vector.
  //! \param sol Global primary solution vector in DOF-order
  //! \param[in] vec Local solution vector associated with specified patch
  //! \param[in] pch The patch to inject solution vector for
  //! \param[in] nndof Number of DOFs per node (optional)
  //! \param[in] basis Basis to inject for (optional)
  bool injectPatchSolution(RealArray& sol, const RealArray& vec,
                           const ASMbase* pch,
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
  bool setPatchMaterial(size_t patch) const;

  //! \brief Returns a const reference to our global-to-local node mapping.
  const std::map<int,int>& getGlob2LocMap() const { return myGlb2Loc; }

  //! \brief Returns the beginning of the property array.
  PropertyVec::const_iterator begin_prop() const { return myProps.begin(); }
  //! \brief Returns the end of the property array.
  PropertyVec::const_iterator end_prop() const { return myProps.end(); }

  //! \brief Returns current Rayleigh system damping matrix.
  //! \param[in] iM Index of the system mass matrix
  //! \param[in] iK Index of the system stiffness matrix
  SystemMatrix* getRayleighDampingMatrix(size_t iM = 1, size_t iK = 0) const;
  //! \brief Returns current system left-hand-side matrix.
  SystemMatrix* getLHSmatrix(size_t idx = 0, bool copy = false) const;
  //! \brief Returns current system right-hand-side vector.
  SystemVector* getRHSvector(size_t idx = 0, bool copy = false) const;
  //! \brief Adds a system vector to the given right-hand-side vector.
  void addToRHSvector(size_t idx, const SystemVector& vec, double scale = 1.0);

  //! \brief Returns a scalar function associated with \a code.
  RealFunc* getSclFunc(int code) const;

  //! \brief Sets the multi-dimension simulator sequence flag.
  void setMDflag(char flag) { mdFlag = flag; }
  //! \brief Returns \e true, if this is the equation system owner.
  bool isFirst() const { return mdFlag <= 1; }

  //! \brief Dumps left-hand-side matrix and right-hand-side vector to file.
  void dumpEqSys(bool initialBlankLine = false);
  //! \brief Dumps a solution vector to file.
  void dumpSolVec(const Vector& x,bool isExpanded = true, bool expOnly = false);

  //! \brief Prints out nodal reaction forces to the log stream.
  virtual int printNRforces(const std::vector<int>& = {}) const { return 0; }

protected:
  //! \brief Returns the multi-dimension simulator sequence flag.
  char getMDflag() const { return mdFlag; }
  //! \brief Shifts global node and element numbers by constant offsets.
  virtual void shiftGlobalNums(int, int) {}

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
  virtual double externalEnergy(const Vectors& psol, const TimeDomain&) const;

  //! \brief Applies app-specific post-processing on norms.
  virtual bool postProcessNorms(Vectors&, Matrix*) { return true; }

  //! \brief Generates element groups for multi-threading of boundary integrals.
  //! \param[in] p Property object identifying a patch boundary
  //! \param[in] silence If \e true, suppress threading group outprint
  void generateThreadGroups(const Property& p, bool silence = false);

  //! \brief Called if number of threads changes.
  void changeNumThreads();

  //! \brief Adds a MADOF with an extraordinary number of DOFs on a given basis.
  //! \param[in] basis The basis to specify number of DOFs for
  //! \param[in] nndof Number of nodal DOFs on the given basis
  //! \param[in] other If \e true, include other bases in MADOF as well
  bool addMADOF(unsigned char basis, unsigned char nndof, bool other = true);

  //! \brief Returns a pointer to the external energy path integral value.
  double* theExtEnerg() { return &extEnergy; }
  //! \brief Returns a const pointer to the external energy path integral value.
  const double* getExtEnerg() const { return &extEnergy; }

private:
  //! \brief Returns an extraordinary MADOF array.
  //! \param[in] basis The basis to specify number of DOFs for
  //! \param[in] nndof Number of nodal DOFs on the given basis
  const std::vector<int>& getMADOF(unsigned char basis,
                                   unsigned char nndof) const;

public:
  static bool ignoreDirichlet; //!< Set to \e true for free vibration analysis
  static bool preserveNOrder;  //!< Set to \e true to preserve node ordering

protected:
  //! \brief Scalar field container
  typedef std::map<int,RealFunc*>     SclFuncMap;
  //! \brief Vector field container
  typedef std::map<int,VecFunc*>      VecFuncMap;
  //! \brief Traction field container
  typedef std::map<int,TractionFunc*> TracFuncMap;

  //! \brief Property code to integrand map
  typedef std::multimap<int,IntegrandBase*> IntegrandMap;

  // Model attributes
  unsigned char  nsd;       //!< Number of spatial dimensions
  PatchVec       myModel;   //!< The actual NURBS/spline model
  PropertyVec    myProps;   //!< Physical property mapping
  SclFuncMap     myScalars; //!< Scalar property fields
  VecFuncMap     myVectors; //!< Vector property fields
  TracFuncMap    myTracs;   //!< Traction property fields
  IntegrandBase* myProblem; //!< The main integrand of this simulator
  IntegrandMap   myInts;    //!< Set of all integrands involved
  AnaSol*        mySol;     //!< Analytical/Exact solution
  FunctionBase*  dualField; //!< Dual solution field (extraction function)

  std::vector<FunctionBase*> extrFunc; //!< Extraction functions for VCP

  //! \brief Total number of unique nodes in the model
  //! \details This variable is only used to carry the number of nodes in
  //! the model, in the case that renumberNodes() is called before preprocess()
  //! and until the SAMpatch::init() method is invoked.
  int  nGlbNodes; //!< This value will be equal to SAM::getNoNodes('A').
  int  isRefined; //!< Indicates if the model is adaptively refined
  bool lagMTOK;   //!< Indicates if global multipliers is OK with multithreading
  bool fixZeros;  //!< If \e true, constrain zero pivots before solving

  //! \brief A struct with data for system matrix/vector dumps.
  struct DumpData
  {
    std::string   fname;  //!< File name
    LinAlg::StorageFormat format; //!< File format flag
    std::set<int> step;   //!< Dump step identifiers
    bool          expand; //!< If \e true, dump expanded solution vectors
    int           count;  //!< Internal step counter, dump only when step==count
    double        eps;    //!< Zero tolerance for printing small values

    //! \brief Default constructor.
    DumpData() : format(LinAlg::FLAT), expand(false), count(0), eps(1.0e-6) {}

    //! \brief Checks if the matrix or vector should be dumped now.
    bool doDump() { return !fname.empty() && step.find(++count) != step.end(); }
  };

  // Post-processing attributes
  std::vector<DumpData> lhsDump; //!< Coefficient matrix dump specifications
  std::vector<DumpData> rhsDump; //!< Right-hand-side vector dump specifications
  std::vector<DumpData> solDump; //!< Solution vector dump specifications

  // Parallel computing attributes
  int               nGlPatches; //!< Number of global patches
  std::vector<int>  myPatches;  //!< Global patch numbers for current processor
  std::vector<int>  myLoc2Glb;  //!< Local-to-global node number mapping
  std::map<int,int> myGlb2Loc;  //!< Global-to-local node number mapping
  const std::map<int,int>* g2l; //!< Pointer to global-to-local node mapping
  std::map<int,int> myDegenElm; //!< Degenerated elements mapping

  // Equation solver attributes
  AlgEqSystem*  myEqSys;     //!< The actual linear equation system
  SAM*          mySam;       //!< Auxiliary data for FE assembly management
  LinSolParams* mySolParams; //!< Input parameters for PETSc
  LinSolParams* myGl2Params; //!< Input parameters for PETSc, for L2 projection

private:
  size_t nIntGP; //!< Number of interior integration points in the whole model
  size_t nBouGP; //!< Number of boundary integration points in the whole model
  size_t nDofS;  //!< Number of degrees of freedom in this sub-simulator
  char   mdFlag; //!< Sequence flag for multi-dimensional simulators

  //! Additional MADOF arrays with varying DOF counts per node
  mutable std::map<int,std::vector<int>> extraMADOFs;

  mutable double extEnergy;  //!< Path integral of external forces
  mutable Vector prevForces; //!< Reaction forces of previous time step
};

#endif
