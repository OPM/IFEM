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
#include "SystemMatrix.h"
#include "TimeDomain.h"
#include "Property.h"
#include "Function.h"
#include "Vec3.h"
#include <map>

class ASMbase;
class IntegrandBase;
class AnaSol;
class VTF;
class SAMpatch;
class AlgEqSystem;
class LinSolParams;


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

class SIMbase : public SIMinput
{
protected:
  //! \brief The constructor initializes the pointers to dynamic data members.
  SIMbase();

public:
  //! \brief The destructor frees the dynamically allocated objects.
  virtual ~SIMbase();


  // Model input and pre-processing methods
  // ======================================

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \param[in] ignoredPatches Indices of patches to ignore in the analysis
  //! \param[in] fixDup Merge duplicated FE nodes on patch interfaces?
  bool preprocess(const std::vector<int>& ignoredPatches = std::vector<int>(),
		  bool fixDup = false);

  //! \brief Read patches from given stream
  //! \param[in] isp The stream to read from
  virtual void readPatches(std::istream& isp) {};

  //! \brief Allocates the system matrices of the FE problem to be solved.
  //! \param[in] mType The matrix format to use
  //! \param[in] nMats Number of system matrices
  //! \param[in] nVec Number of system right-hand-side vectors
  bool initSystem(SystemMatrix::Type mType, size_t nMats, size_t nVec);

  //! \brief Associates a system vector to a system matrix.
  //! \sa AlgEqSystem::setAssociatedVector
  //! \param[in] iMat Index of a coefficient matrix
  //! \param[in] iVec Index of the system vector to associate with the matrix
  bool setAssociatedRHS(size_t iMat, size_t iVec);

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  bool setMode(int mode);

  //! \brief Defines the spatial numerical integration scheme to use.
  virtual void setQuadratureRule(size_t) {}

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
  //! \brief Returns the number of solution vectors.
  size_t getNoSolutions() const;
  //! \brief Returns the total number of patches in the model.
  int getNoPatches() const { return nGlPatches; }
  //! \brief Returns the number of registered result points.
  size_t getNoResultPoints() const { return myPoints.size(); }

  //! \brief Initializes time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  bool initDirichlet(double time = 0.0);

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
  bool updateDirichlet(double time = 0.0, const Vector* prevSol = 0);


  // Computational methods
  // =====================

  //! \brief Administers assembly of the linear equation system.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] pSol Previous primary solution vectors in DOF-order
  //! \param[in] newLHSmatrix If \e false, only integrate the RHS vector
  bool assembleSystem(const TimeDomain& time, const Vectors& pSol = Vectors(),
		      bool newLHSmatrix = true);

  //! \brief Administers assembly of the linear equation system.
  //! \param[in] pSol Previous primary solution vectors in DOF-order
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
  virtual bool solveSystem(Vector& solution, int printSol = 0,
			   const char* compName = "displacement",
			   bool newLHS = true);

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
  virtual double solutionNorms(const Vector& x, double* inf,
			       size_t* ind, size_t nf = 0) const;

  //! \brief Integrates some solution norm quantities.
  //! \details If an analytical solution is provided, norms of the exact
  //! error in the solution are computed as well.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations.
  //! \param[in] psol Primary solution vectors
  //! \param[out] gNorm Global norm quantities
  //! \param[out] eNorm Element-wise norm quantities
  virtual bool solutionNorms(const TimeDomain& time, const Vectors& psol,
			     Vector& gNorm, Matrix* eNorm = 0);

  //! \brief Integrates some solution norm quantities.
  //! \details If an analytical solution is provided, norms of the exact
  //! error in the solution are computed as well.
  //! \param[in] psol Primary solution vectors
  //! \param[out] eNorm Element-wise norm quantities
  //! \param[out] gNorm Global norm quantities
  //!
  //! \details Use this version for linear/stationary problems only.
  virtual bool solutionNorms(const Vectors& psol, Matrix& eNorm, Vector& gNorm)
  { return this->solutionNorms(TimeDomain(),psol,gNorm,&eNorm); }

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

  //! \brief Projects the secondary solution associated with \a psol.
  //! \param psol Control point values of primary(in)/secondary(out) solution
  //!
  //! \details The secondary solution, defined through the Integrand object,
  //! corresponding to the primary solution \a psol is projected onto the
  //! spline basis to obtain the control point values of the secondary solution.
  bool project(Vector& psol);

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
  //! \param[in] inpFile File name used to construct the VTF-file name from
  //! \param[in] nViz    Number of visualization points over a knot-span
  //! \param[in] format  Format of VTF-file (0=ASCII, 1=BINARY)
  //!
  //! \details The spline patches are tesselated into linear finite elements
  //! with a fixed number of elements within each knot-span of non-zero length.
  //! The solution fields are then evaluated at the nodal points of the
  //! generated FE mesh and written to the VTF-file as vector and scalar fields
  //! by the other \a writeGlv* methods.
  bool writeGlv(const char* inpFile, const int* nViz, int format);

  //! \brief Writes boundary conditions as scalar fields to the VTF-file.
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param nBlock Running result block counter
  bool writeGlvBC(const int* nViz, int& nBlock) const;

  //! \brief Writes boundary tractions for a given time step to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvT(int iStep, int& nBlock) const;

  //! \brief Writes a vector field for a given load/time step to the VTF-file.
  //! \param[in] vec The vector field to output (nodal values in DOF-order)
  //! \param[in] fieldName Name identifying the vector field
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvV(const Vector& vec, const char* fieldName,
		 const int* nViz, int iStep, int& nBlock) const;

  //! \brief Writes solution fields for a given load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] psolOnly If \e true, skip secondary solution field evaluation
  //! \param[in] pvecName Optional name of the primary vector field solution
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  //!
  //! \details The primary vector field solution is written as a deformation
  //! plot if \a pvecName is NULL. Otherwise, it is written as a named vector
  //! field. If an analytical solution is provided, the exact secondary solution
  //! fields are written to the VTF-file as well.
  virtual bool writeGlvS(const Vector& psol, const int* nViz, int iStep,
			 int& nBlock, double time = 0.0, char psolOnly = 0,
			 const char* pvecName = 0, int idBlock = 10,
			 int psolComps = 0);

  //! \brief Writes a mode shape and associated eigenvalue to the VTF-file.
  //! \details The eigenvalue is used as a label on the step state info
  //! that is associated with the eigenvector.
  //! \param[in] mode The mode shape eigenvector to output
  //! \param[in] freq \e true if the eigenvalue is a frequency
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param nBlock Running result block counter
  bool writeGlvM(const Mode& mode, bool freq, const int* nViz, int& nBlock);

  //! \brief Writes element norms for a given load/time step to the VTF-file.
  //! \details This method can be used only when the number of visualization
  //! points over each knot-span equals 2 (that is, no additonal points).
  //! \param[in] norms The element norms to output
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvN(const Matrix& norms, int iStep, int& nBlock);

  //! \brief Writes time/load step info to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] value Time or load parameter of the step
  //! \param[in] itype Type indentifier of the step
  bool writeGlvStep(int iStep, double value = 0.0, int itype = 0);

  //! \brief Closes the current VTF-file.
  void closeGlv();

  //! \brief Dumps the (possibly refined) geometry in g2-format.
  //! \param os Output stream to write the geometry data to
  bool dumpGeometry(std::ostream& os) const;
  //! \brief Dumps the (possibly refined) spline basis in g2-format.
  //! \param os Output stream to write the spline data to
  //! \param[in] basis Which basis to dump for mixed methods (0 = geometry)
  //! \param[in] patch Which patch to dump for (0 = all)
  bool dumpBasis(std::ostream& os, int basis = 0, size_t patch=0) const;
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
  //! \brief Dumps the primary solution in ASCII format for inspection.
  //! \param[in] psol Primary solution vector
  //! \param os Output stream to write the solution data to
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpPrimSol(const Vector& psol, std::ostream& os,
		   bool withID = true) const;

protected:
  //! \brief Defines the type of a property set.
  //! \param[in] code The property code to be associated with the property type
  //! \param[in] ptype The property type to be associated with the given code
  //! \param[in] pindex 0-based index into problem-dependent property container
  bool setPropertyType(int code, Property::Type ptype, int pindex = -1);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homegeneous Dirichlet condition property code
  virtual bool addConstraint(int patch, int lndx, int ldim,
			     int dirs, int code = 0) = 0;

  //! \brief Creates the computational FEM model from the spline patches.
  bool createFEMmodel();

public:
  //! \brief Returns the local patch index for the given global patch number.
  //! \details For serial applications this is an identity mapping only, whereas
  //! for parallel applications the local (1-based) patch index on the current
  //! processor is returned. If \a patchNo is out of range, -1 is returned.
  //! If \a patchNo is not on current processor, 0 is returned.
  int getLocalPatchIndex(int patchNo) const;

  //! \brief Extracts the local solution vector for a specified patch.
  //! \param[in] sol Global primary solution vectors in DOF-order
  //! \param[in] pindx Local patch index to extract solution vectors for
  //! \return Total number of DOFs in the patch (first basis only if mixed)
  //!
  //! \details The extracted patch-level solution vector is stored within the
  //! Integrand \a *myProblem such that \a evalSolution can be invoked to get
  //! the secondary solution field values within the same patch afterwards,
  size_t extractPatchSolution(const Vector& sol, int pindx);

  //! \brief Injects a patch-wise solution vector into the global vector.
  //! \param sol Global primary solution vector in DOF-order
  //! \param[in] locSol Local solution vector associated with specified patch
  //! \param[in] pindx Local patch index to inject solution vector for
  //! \param[in] nndof Number of DOFs per node (optional)
  bool injectPatchSolution(Vector& sol, const Vector& locSol,
			   int pindx, unsigned char nndof = 0);

protected:
  //! \brief Extracts local solution vector(s) for the given patch.
  //! \param[in] patch Pointer to the patch to extract solution vectors for
  //! \param[in] sol Global primary solution vectors in DOF-order
  void extractPatchSolution(const ASMbase* patch, const Vectors& sol);

  //! \brief Initializes material properties for integration of interior terms.
  virtual bool initMaterial(size_t) { return true; }

  //! \brief Initializes for integration of Neumann terms for a given property.
  virtual bool initNeumann(size_t) { return true; }

  //! \brief Reads a LinSolParams object from the given stream.
  //! \details This method helps with encapsulating PETSc in libIFEM.
  void readLinSolParams(std::istream& is, int npar);

  //! \brief Finalizes the system assembly.
  bool finalizeAssembly(bool newLHSmatrix);

public:
  //! \brief Enum defining the available discretization methods.
  enum Discretization { Spline, Lagrange, Spectral };

  static Discretization discretization; //!< Spatial discretization option

  static bool ignoreDirichlet; //!< Set to \e true for free vibration analysis

  static int num_threads_SLU; //!< Number of threads for SuperLU_MT

protected:
  //! \brief Spline patch container
  typedef std::vector<ASMbase*>       FEModelVec;
  //! \brief Scalar field container
  typedef std::map<int,RealFunc*>     SclFuncMap;
  //! \brief Vector field container
  typedef std::map<int,VecFunc*>      VecFuncMap;
  //! \brief Traction field container
  typedef std::map<int,TractionFunc*> TracFuncMap;

  // Model attributes
  FEModelVec     myModel;   //!< The actual NURBS/spline model
  PropertyVec    myProps;   //!< Physical property mapping
  SclFuncMap     myScalars; //!< Scalar property fields
  VecFuncMap     myVectors; //!< Vector property fields
  TracFuncMap    myTracs;   //!< Traction property fields
  IntegrandBase* myProblem; //!< Problem-specific data and methods
  AnaSol*        mySol;     //!< Analytical/Exact solution
  VTF*           myVtf;     //!< VTF-file for result visualization
  ResPointVec    myPoints;  //!< User-defined result sampling points

  // Parallel computing attributes
  int              nGlPatches; //!< Number of global patches
  std::vector<int> myPatches;  //!< Global patch numbers for current processor

  // Equation solver attributes
  AlgEqSystem*  myEqSys;     //!< The actual linear equation system
  SAMpatch*     mySam;       //!< Auxiliary data for FE assembly management
  LinSolParams* mySolParams; //!< Input parameters for PETSc
};

#endif
