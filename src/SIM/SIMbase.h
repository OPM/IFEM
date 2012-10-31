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

#include "SIMinput.h"
#include "SIMdependency.h"
#include "TimeDomain.h"
#include "TopologySet.h"
#include "Property.h"
#include "Function.h"
#include "MatVec.h"
#include "Vec3.h"

class IntegrandBase;
class NormBase;
class ForceBase;
class AnaSol;
class VTF;
class SAM;
class AlgEqSystem;
class LinSolParams;
class TimeStep;

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
  \brief Struct defining a result sampling point.
*/

struct ResultPoint
{
  unsigned char npar;   //!< Number of parameters
  size_t        patch;  //!< Patch index [0,nPatch>
  int           inod;   //!< Local node number of the closest node
  double        par[3]; //!< Parameters of the point (u,v,w)
  Vec3          X;      //!< Spatial coordinates of the point
  // \brief Default constructor.
  ResultPoint() : npar(0), patch(0), inod(0) { par[0] = par[1] = par[2] = 0.0; }
};

typedef std::vector<ResultPoint> ResPointVec; //!< Result point container


/*!
  \brief Base class for NURBS-based FEM simulators.
  \details This class incapsulates data and methods need for solving
  partial differential equations using NURBS-based finite elements.
  It only contains the problem-independent data and methods.
  Sub-classes are derived with additional info regarding the problem to solve.
*/

class SIMbase : public SIMinput, public SIMdependency
{
protected:
  //! \brief The constructor initializes the pointers to dynamic data members.
  SIMbase();

public:
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~SIMbase();


  // Model input and pre-processing methods
  // ======================================

  //! \brief Initializes the property containers of the model.
  //! \details Use this method to clear the model before re-reading
  //! the input file in the refinement step of an adaptive simulation.
  virtual void clearProperties();

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an xml document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const;

protected:
  //! \brief Parses the "set" attribute of a material XML-tag.
  //! \param[in] elem The XML element extract the set name from
  //! \param[in] mindex Index into problem-dependent material property container
  //! \return The property code to be associated with the material
  int parseMaterialSet(const TiXmlElement* elem, int mindex);

private:
  //! \brief Parses a subelement of the \a geometry XML-tag.
  bool parseGeometryTag(const TiXmlElement* elem);
  //! \brief Parses a subelement of the \a boundaryconditions XML-tag.
  bool parseBCTag(const TiXmlElement* elem);
  //! \brief Parses a subelement of the \a resultoutput XML-tag.
  bool parseOutputTag(const TiXmlElement* elem);
  //! \brief Parses a subelement of the \a linearsolver XML-tag.
  bool parseLinSolTag(const TiXmlElement* elem);

public:
  //! \brief Creates the FE model by copying the given patches.
  virtual void clonePatches(const PatchVec&, const std::map<int,int>&) {}

  //! \brief Refines a list of elements.
  virtual bool refine(const std::vector<int>&, const std::vector<int>&,
		      const char* = 0) { return false; }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \param[in] ignored Indices of patches to ignore in the analysis
  //! \param[in] fixDup Merge duplicated FE nodes on patch interfaces?
  virtual bool preprocess(const std::vector<int>& ignored = std::vector<int>(),
			  bool fixDup = false);

  //! \brief Defines a vector field property.
  //! \param[in] code The property code to be associated with the property
  //! \param[in] ptype The property type to be associated with the given code
  //! \param[in] field The vector field representing the physical property
  //! \param[in] pflag Flag for local axis directions (see setPropertyType)
  size_t setVecProperty(int code, Property::Type ptype, VecFunc* field = NULL,
			int pflag = -1);

  //! \brief Defines a traction field property.
  //! \param[in] code The property code to be associated with the property
  //! \param[in] ptype The property type to be associated with the given code
  //! \param[in] field The traction field representing the physical property
  bool setTracProperty(int code, Property::Type ptype, TractionFunc* field = 0);

  //! \brief Allocates the system matrices of the FE problem to be solved.
  //! \param[in] mType The matrix format to use
  //! \param[in] nMats Number of system matrices
  //! \param[in] nVec Number of system right-hand-side vectors
  //! \param[in] withRF Whether nodal reaction forces should be computed or not
  bool initSystem(int mType, size_t nMats, size_t nVec, bool withRF = true);

  //! \brief Associates a system vector to a system matrix.
  //! \sa AlgEqSystem::setAssociatedVector
  //! \param[in] iMat Index of a coefficient matrix
  //! \param[in] iVec Index of the system vector to associate with the matrix
  bool setAssociatedRHS(size_t iMat, size_t iVec);

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  //! \param[in] resetSol If \e true, the internal solution vectors are cleared
  bool setMode(int mode, bool resetSol = false);

  //! \brief Defines the spatial numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  //! \param[in] redimBuffers Toggle initialization of internal buffer arrays
  void setQuadratureRule(size_t ng, bool redimBuffers = false);

  //! \brief Prints out problem-specific data to the given stream.
  void printProblem(std::ostream& os) const;

  //! \brief Returns a pointer to the problem-specific data object.
  const IntegrandBase* getProblem() const { return myProblem; }

  //! \brief Returns the name of this simulator.
  //! \details This method is typically reimplemented in sub-classes that are
  //! parts of a partitioned solution method and are used to identify the basis
  //! for the result fields associated with each simulator in the HDF5 output.
  virtual std::string getName() const { return "SIMbase"; }

  //! \brief Returns the number of primary solution fields.
  //! \param[in] basis Which basis to condsider when mixed methods (0 = both)
  size_t getNoFields(int basis = 0) const;
  //! \brief Returns the number of spatial dimensions in the model.
  size_t getNoSpaceDim() const;
  //! \brief Returns the model size in terms of number of DOFs.
  size_t getNoDOFs() const;
  //! \brief Returns the model size in terms of number of (unique) nodes.
  size_t getNoNodes(bool unique = false) const;
  //! \brief Returns the model size in terms of number of elements.
  size_t getNoElms() const;
  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions() const;
  //! \brief Returns the total number of patches in the model.
  int getNoPatches() const { return nGlPatches; }
  //! \brief Returns the number of registered result points.
  size_t getNoResultPoints() const { return myPoints.size(); }
  //! \brief Returns the visualization dump interval.
  int getDumpInterval() const { return opt.saveInc; }

  //! \brief Returns the type (DOF classification) of the specified node.
  char getNodeType(int inod) const;

  //! \brief Initializes time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  bool initDirichlet(double time = 0.0);

  //! \brief Initializes for time-dependent simulation.
  virtual void init(TimeStep&) {}
  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep&) { return false; }
  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep&) { return false; }
  //! \brief Saves the converged results to VTF file of a given time step.
  virtual bool saveStep(int, double, int&) { return false; }

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
  virtual bool updateDirichlet(double time = 0.0, const Vector* prevSol = 0);

  //! \brief Updates problem-dependent state based on the current solution.
  //! \param[in] solution Current primary solution vector in DOF-order
  virtual bool updateConfiguration(const Vector& solution) { return true; }

  //! \brief Updates the grid coordinates.
  //! \param[in] displ The displacement increment to update the grid with
  bool updateGrid(const Vector& displ);


  // Computational methods
  // =====================

  //! \brief Administers assembly of the linear equation system.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] pSol Previous primary solution vectors in DOF-order
  //! \param[in] newLHSmatrix If \e false, only integrate the RHS vector
  //! \param[in] poorConvg If \e true, the nonlinear driver is converging poorly
  bool assembleSystem(const TimeDomain& time, const Vectors& pSol,
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

  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] newLHS If \e false, reuse the LHS-matrix from previous call.
  bool solveSystem(Vector& solution, int printSol = 0,
                   const char* compName = "displacement", bool newLHS = true);

  //! \brief Finds the worst energy DOFs in the residual.
  //! \param[in] x Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[in] nWorst How many bad DOFs to detect
  //! \param[in] eps Only record the energies larger than this tolerance
  //! \param[out] worst Node and local DOF number and values of the worst DOFs
  void getWorstDofs(const Vector& x, const Vector& r, size_t nWorst, double eps,
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
  //! \return L2-norm of the solution vector
  double solutionNorms(const Vector& x, double* inf,
                       size_t* ind, size_t nf = 0) const;

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
		     Vectors& gNorm, Matrix* eNorm = 0);
  //! \brief Integrates some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations.
  //! \param[in] psol Primary solution vectors
  //! \param[out] gNorm Global norm quantities
  //! \param[out] eNorm Element-wise norm quantities
  //!
  //! \details Use this version if no projected solutions are needed/available.
  bool solutionNorms(const TimeDomain& time, const Vectors& psol,
		     Vectors& gNorm, Matrix* eNorm = 0)
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

  //! \brief Prints integrated solution norms to the specified output stream.
  virtual std::ostream& printNorms(const Vectors&, std::ostream& os)
  { return os; }

  //! \brief Computes the total reaction forces in the model.
  //! \param[out] RF Reaction force in each spatial direction + energy
  //! \param[in] psol Primary solution vector
  bool getCurrentReactions(RealArray& RF, const Vector& psol) const;

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

  //! \brief Evaluates the secondary solution field for specified patch.
  //! \param[out] field Control point values of the secondary solution field
  //! \param[in] pindx Local patch index to evaluate solution field for
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the Integrand object.
  //! The solution is evaluated at the Greville points and then projected onto
  //! the spline basis to obtain the control point values.
  bool evalSecondarySolution(Matrix& field, int pindx) const;


  // Post-processing methods
  // =======================

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fnam File name used to construct the VTF-file name from
  bool writeGlv(const char* fnam) { int n = 0; return this->writeGlvG(n,fnam); }

  //! \brief Writes current model geometry to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] inpFile File name used to construct the VTF-file name from
  //!
  //! \details The spline patches are tesselated into linear finite elements
  //! with a fixed number of elements within each knot-span of non-zero length.
  //! The solution fields are then evaluated at the nodal points of the
  //! generated FE mesh and written to the VTF-file as vector and scalar fields
  //! by the other \a writeGlv* methods.
  virtual bool writeGlvG(int& nBlock, const char* inpFile = 0);

  //! \brief Writes additional, problem-specific, results to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] iStep Load/time step identifier
  virtual bool writeGlvA(int& nBlock, int iStep = 1) const { return true; }

  //! \brief Writes boundary conditions as scalar fields to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] iStep Load/time step identifier
  bool writeGlvBC(int& nBlock, int iStep = 1) const;

  //! \brief Writes boundary tractions for a given time step to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvT(int iStep, int& nBlock) const;

  //! \brief Writes a vector field for a given load/time step to the VTF-file.
  //! \param[in] vec The vector field to output (nodal values in DOF-order)
  //! \param[in] fieldName Name identifying the vector field
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  bool writeGlvV(const Vector& vec, const char* fieldName,
		 int iStep, int& nBlock, int idBlock = 2) const;

  //! \brief Writes solution fields for a given load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] psolOnly If \e true, skip secondary solution field evaluation
  //! \param[in] pvecName Optional name of the primary vector field solution
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  bool writeGlvS(const Vector& psol, int iStep, int& nBlock, double time = 0.0,
		 char psolOnly = 0, const char* pvecName = 0,
		 int idBlock = 10, int psolComps = 0);

  //! \brief Writes projected solutions for a given time step to the VTF-file.
  //! \param[in] ssol Secondary solution vector (control point values)
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] prefix Common prefix for the field components
  bool writeGlvP(const Vector& ssol, int iStep, int& nBlock,
		 int idBlock = 100, const char* prefix = "Global projected");

  //! \brief Writes a mode shape to the VTF-file.
  //! \param[in] mode The mode shape eigenvector and associated eigenvalue
  //! \param[in] freq \e true if the eigenvalue is a frequency
  //! \param nBlock Running result block counter
  //!
  //! \details The eigenvalue is used as a label on the step state info.
  bool writeGlvM(const Mode& mode, bool freq, int& nBlock);

  //! \brief Writes element norms for a given load/time step to the VTF-file.
  //! \param[in] norms The element norms to output
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] prefix Prefices for projected solutions
  bool writeGlvN(const Matrix& norms, int iStep, int& nBlock,
		 const char** prefix = 0);
  //! \brief Writes a scalar function to the VTF-file.
  //! \param[in] f The function to output
  //! \param[in] fname Name of the function
  //! \param[in] iStep Load/time step identifier
  //! \param[in] idBlock Starting value of result block numbering
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] time Load/time step parameter
  bool writeGlvF(const RealFunc& f, const char* fname,
		 int iStep, int& nBlock, int idBlock = 50, double time = 0.0);

  //! \brief Writes time/load step info to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] value Time or load parameter of the step
  //! \param[in] itype Type identifier of the step
  bool writeGlvStep(int iStep, double value = 0.0, int itype = 0);

  //! \brief Closes the current VTF-file.
  void closeGlv();

  //! \brief Returns the current VTF-file object.
  VTF* getVTF() const { return myVtf; }
  //! \brief Defines the VTF-file for subsequent results output.
  void setVTF(VTF* vtf) { myVtf = vtf; }

  //! \brief Dumps the (possibly refined) geometry in g2-format.
  //! \param os Output stream to write the geometry data to
  bool dumpGeometry(std::ostream& os) const;
  //! \brief Dumps the (possibly refined) spline basis in g2-format.
  //! \param os Output stream to write the spline data to
  //! \param[in] basis Which basis to dump for mixed methods (0 = geometry)
  //! \param[in] patch Which patch to dump for (0 = all)
  bool dumpBasis(std::ostream& os, int basis = 0, size_t patch = 0) const;
  //! \brief Dumps the primary solution in ASCII format for inspection.
  //! \param[in] psol Primary solution vector
  //! \param os Output stream to write the solution data to
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpPrimSol(const Vector& psol, std::ostream& os,
                   bool withID = true) const;
  //! \brief Dumps the entire solution in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param os Output stream to write the solution data to
  bool dumpSolution(const Vector& psol, std::ostream& os) const;
  //! \brief Dumps solution results at specified points in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] formatted If \e false, write all result points on a single line
  //!            without point identifications, but with time as first column
  //! \param[in] precision Number of digits after the decimal point
  bool dumpResults(const Vector& psol, double time, std::ostream& os,
		   bool formatted = true, std::streamsize precision = 3) const;
  //! \brief Dumps coordinate at specified points in ASCII format.
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] formatted If \e false, write all result points on a single line
  //!            without point identifications, but with time as first column
  //! \param[in] precision Number of digits after the decimal point
  bool dumpResultCoords(double time, std::ostream& os, bool formatted = true,
			std::streamsize precision = 3) const;
  //! \brief Dumps additional problem-specific results in ASCII format.
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] precision Number of digits after the decimal point
  virtual void dumpMoreResults(double time, std::ostream& os,
                               std::streamsize precision = 3) const {}

  //! \brief Evaluates the secondary solution for a given load/time step.
  //! \param[in] psol Primary solution vector
  //! \param[in] time Load/time step parameter
  //!
  //! \details This method only evaluates the solutions fields, but does not
  //! return any data. The method is used only for load/time steps that are not
  //! not be saved, but the solution has to be evaluated at every increment in
  //! any case to ensure consistency (like, when constitutive models with
  //! history variables are in use).
  bool eval2ndSolution(const Vector& psol, double time);

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
  ForceBase* getBoundaryForceIntegrand(const Vec3* X0 = NULL) const;

  //! \brief Returns a unique integer code for a Property set.
  //! \param[in] setName Name of the topology set the property is defined on
  //! \param[in] comp The solution components on which the property is applied
  //!
  //! \details The actual Property objects are also created (one for each entity
  //! in the topology set) and their type is set to UNDEFINED. The method
  //! setPropertyType must be used to assign the actual Property type.
  int getUniquePropertyCode(const std::string& setName, int comp = 0);
  //! \brief Creates a set of Property objects.
  //! \param[in] setName Name of the topology set the property is defined on
  //! \param[in] pc The property code to be associated with this set
  bool createPropertySet(const std::string& setName, int pc);

protected:
  //! \brief Defines the type of a property set.
  //! \param[in] code The property code to be associated with the property type
  //! \param[in] ptype The property type to be associated with the given code
  //! \param[in] pindex 0-based index into problem-dependent property container
  size_t setPropertyType(int code, Property::Type ptype, int pindex = -1);

  //! \brief Defines a Neumann boundary condition property by parsing a string.
  //! \param[in] prop The string to parse the property definition from
  //! \param[in] type Additional option defining the type of property definition
  //! \param[in] ndir Direction of the surface traction on the Neumann boundary
  //! \param[in] code The property code to be associated with this property
  bool setNeumann(const std::string& prop, const std::string& type,
		  int ndir, int code);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homegeneous Dirichlet condition property code
  //! \param ngnod Total number of global nodes in the model (might be updated)
  virtual bool addConstraint(int patch, int lndx, int ldim,
			     int dirs, int code, int& ngnod) = 0;

  //! \brief Specialized preprocessing performed before assembly initialization.
  virtual void preprocessBeforeAsmInit(int&) {}

  //! \brief Creates the computational FEM model from the spline patches.
  //! \param[in] resetNumb If \e true, start element and node numbers from zero
  bool createFEMmodel(bool resetNumb = true);

private:
  //! \brief Extracts all local solution vector(s) for a specified patch.
  //! \param[in] problem Integrand to receive the patch-level solution vectors
  //! \param[in] sol Global primary solution vectors in DOF-order
  //! \param[in] pindx Local patch index to extract solution vectors for
  //!
  //! \details This method is typically invoked before the \a integrate method
  //! on the the specified path, in order to extract all patch-level vector
  //! quantities needed by the Integrand. This also includes any dependent
  //! vectors from other simulator classes that have been registered.
  //! All patch-level vectors are stored within the provided Integrand \a *p.
  bool extractPatchSolution(IntegrandBase* problem,
                            const Vectors& sol, size_t pindx);

public:
  //! \brief Extracts all local solution vector(s) for a specified patch.
  //! \param[in] sol Global primary solution vectors in DOF-order
  //! \param[in] pindx Local patch index to extract solution vectors for
  bool extractPatchSolution(const Vectors& sol, size_t pindx)
  { return this->extractPatchSolution(myProblem,sol,pindx); }

  //! \brief Extracts a local solution vector for a specified patch.
  //! \param[in] sol Global primary solution vector in DOF-order
  //! \param[in] pindx Local patch index to extract solution vector for
  //! \return Total number of DOFs in the patch (first basis only if mixed)
  //!
  //! \details The extracted patch-level solution vector is stored within the
  //! Integrand \a *myProblem such that \a evalSolution can be invoked to get
  //! the secondary solution field values within the same patch afterwards.
  size_t extractPatchSolution(const Vector& sol, int pindx);

  //! \brief Injects a patch-wise solution vector into the global vector.
  //! \param sol Global primary solution vector in DOF-order
  //! \param[in] locSol Local solution vector associated with specified patch
  //! \param[in] pindx Local patch index to inject solution vector for
  //! \param[in] nndof Number of DOFs per node (optional)
  bool injectPatchSolution(Vector& sol, const Vector& locSol,
			   int pindx, unsigned char nndof = 0);

  //! \brief Extracts element results for a specified patch.
  //! \param[in] globRes Global element result array
  //! \param[out] elmRes Patch-level element result array
  //! \param[in] pindx Local patch index to extract element results for
  bool extractPatchElmRes(const Matrix& globRes, Matrix& elmRes, int pindx);

  //! \brief Returns the local patch index for the given global patch number.
  //! \details For serial applications this is an identity mapping only, whereas
  //! for parallel applications the local (1-based) patch index on the current
  //! processor is returned. If \a patchNo is out of range, -1 is returned.
  //! If \a patchNo is not on current processor, 0 is returned.
  int getLocalPatchIndex(int patchNo) const;

  //! \brief Returns a const reference to our FEM model.
  const PatchVec& getFEModel() const { return myModel; }
  //! \brief Returns a const reference to our global-to-local node mapping.
  const std::map<int,int>& getGlob2LocMap() const { return myGlb2Loc; }

  //! \brief Returns the beginning of the property array.
  PropertyVec::const_iterator begin_prop() const { return myProps.begin(); }
  //! \brief Returns the end of the property array.
  PropertyVec::const_iterator end_prop() const { return myProps.end(); }

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  virtual ASMbase* readPatch(std::istream& isp, int pchInd) const = 0;

protected:
  //! \brief Initializes material properties for integration of interior terms.
  virtual bool initMaterial(size_t) { return true; }
  //! \brief Initializes the body load properties for current patch.
  virtual bool initBodyLoad(size_t) { return true; }
  //! \brief Initializes for integration of Neumann terms for a given property.
  virtual bool initNeumann(size_t) { return true; }

  //! \brief Reads patches from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[out] vec Array of spline patches that were read
  //! \param[in] whiteSpace For message formatting
  virtual bool readPatches(std::istream& isp, PatchVec& vec,
                           const char* whiteSpace = "") = 0;
  //! \brief Reads global node data for a patch from given input stream.
  //! \param[in] isn The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read node data for
  //! \param[in] basis The basis to read node data for (when mixed FEM)
  //! \param[in] oneBased If \e true the read node numbers are assumed
  //! one-based. If \e false they are assumed to be zero-based.
  virtual bool readNodes(std::istream& isn, int pchInd, int basis = 0,
			 bool oneBased = false) { return false; }
  //! \brief Reads node numbers from given input stream.
  //! \param[in] isn The input stream to read from
  virtual void readNodes(std::istream& isn) {}
  //! \brief Reads a LinSolParams object from the given stream.
  void readLinSolParams(std::istream& is, int npar);

  //! \brief Assembles problem-dependent discrete terms, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase*) { return true; }
  //! \brief Computes (possibly problem-dependent) external energy contribution.
  virtual double externalEnergy(const Vectors& psol) const;

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

  // Model attributes
  PatchVec       myModel;   //!< The actual NURBS/spline model
  PropertyVec    myProps;   //!< Physical property mapping
  TopologySet    myEntitys; //!< Set of named topological entities
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
  ResPointVec myPoints; //!< User-defined result sampling points
  VTF*        myVtf;    //!< VTF-file for result visualization

  std::vector<DumpData> lhsDump; //!< Coefficient matrix dump specifications
  std::vector<DumpData> rhsDump; //!< Right-hand-side vector dump specifications

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
};

#endif
