// $Id$
//==============================================================================
//!
//! \file ASMbase.h
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for spline-based finite element (FE) assembly drivers.
//!
//==============================================================================

#ifndef _ASM_BASE_H
#define _ASM_BASE_H

#include "MatVec.h"
#include "MPCLess.h"
#include "Function.h"
#include <map>
#include <set>

typedef std::vector<int>       IntVec;  //!< General integer vector
typedef std::vector<IntVec>    IntMat;  //!< General 2D integer matrix
typedef std::map<MPC*,int>     MPCMap;  //!< MPC to function code mapping
typedef std::set<MPC*,MPCLess> MPCSet;  //!< Sorted set of MPC equations
typedef MPCSet::const_iterator MPCIter; //!< Iterator over an MPC equation set

struct TimeDomain;
struct ElementBlock;
class GlobalIntegral;
class IntegrandBase;
class Integrand;
class ASMbase;
class Vec3;

typedef std::vector<ASMbase*> ASMVec; //!< Spline patch container


/*!
  \brief Base class for spline-based finite element (FE) assembly drivers.

  \details This class incapsulates the data and methods needed for assembling
  the algebraic equation system resulting from a finite element discretization
  of a set of partial differential equations using splines as basis functions.

  The class does not contain any problem-specific data or methods.
  The methods that need access to problem information are given that through
  Integrand objects that are passed as arguments to those methods.
*/

class ASMbase
{
public:
  //! \brief Struct for boundary condition codes.
  struct BC
  {
    int node; //!< Global node number of the constrained node
    char CX;  //!< Boundary condition code for X-translation
    char CY;  //!< Boundary condition code for Y-translation
    char CZ;  //!< Boundary condition code for Z-translation
    //! \brief Constructor initializing a BC instance.
    BC(int n) : node(n), CX(1), CY(1), CZ(1) {}
  };

  typedef std::vector<BC> BCVec; //!< Nodal boundary condition container

protected:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMbase(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Copy constructor.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  //!
  //! \details The default copy constructor is overridden, such that this patch
  //! shares the FE data (including the basis functions) with the provided
  //! \a patch. The boundary conditions and constraint equations are however
  //! not copied, as the copy constructor is typically used for multi-stage
  //! simulators with different sub-problems discretized on the same grid.
  ASMbase(const ASMbase& patch, unsigned char n_f = 0);

public:
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~ASMbase();

  //! \brief Checks if this patch is empty.
  virtual bool empty() const = 0;

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is) = 0;
  //! \brief Writes the geometry/basis of the patch to the given stream.
  virtual bool write(std::ostream& os, int basis = 0) const = 0;

  //! \brief Generates the finite element topology data for this patch.
  virtual bool generateFEMTopology() = 0;

  //! \brief Clears the contents of this patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Adds extraordinary elements associated with a patch boundary.
  //! \param[in] dim Dimension of the boundary (should be \a nsd - 1)
  //! \param[in] item Local index of the boundary face/edge
  //! \param[in] nXnod Number of extraordinary nodes
  //! \param[out] nodes Global numbers assigned to the extraordinary nodes
  virtual bool addXElms(short int dim, short int item, size_t nXnod,
                        std::vector<int>& nodes);

  //! \brief Adds a set of Lagrange multipliers to the specified element.
  //! \param[in] iel 1-based element index local to current patch
  //! \param[in] mGLag Global node numbers of the Lagrange multipliers
  //! \param[in] nnLag Number of Lagrange multipliers per node
  bool addLagrangeMultipliers(size_t iel, const IntVec& mGLag,
                              unsigned char nnLag = 1);

  //! \brief Defines the numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  void setGauss(int ng) { nGauss = ng; }


  // Service methods for query of various model data
  // ===============================================

  //! \brief Returns the number of spatial dimensions.
  unsigned char getNoSpaceDim() const { return nsd; }
  //! \brief Returns the number of parameter dimensions.
  unsigned char getNoParamDim() const { return ndim; }
  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int b = 0) const { return b > 1 ? 0 : nf; }
  //! \brief Returns the number of Lagrange multipliers per node.
  unsigned char getNoLagPerNode() const { return nLag; }

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction
  virtual bool getOrder(int& p1, int& p2, int& p3) const { return false; }

  //! \brief Returns local 1-based index of the node with given global number.
  //! \details If the given node number is not present, 0 is returned.
  //! \param[in] globalNum Global node number
  virtual size_t getNodeIndex(int globalNum, bool = false) const;
  //! \brief Returns the global node number for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual int getNodeID(size_t inod, bool = false) const;
  //! \brief Returns the global element number for the given element
  //! \param[in] iel 1-based element index local to current patch
  int getElmID(size_t iel) const;
  //! \brief Returns the number of DOFs per node.
  //! \param[in] inod 1-based node index local to current patch
  virtual unsigned char getNodalDOFs(size_t inod) const;
  //! \brief Returns the classification of a node.
  //! \param[in] inod 1-based node index local to current patch
  virtual char getNodeType(size_t inod) const;
  //! \brief Returns \e true if node \a n is a Lagrange multiplier node.
  bool isLMn(size_t n) const { return n >= myLMs.first && n <= myLMs.second; }
  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const = 0;
  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X nsd\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const = 0;
  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel 1-based element index local to current patch
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const = 0;

  //! \brief Prints out the nodal coordinates of this patch to the given stream.
  void printNodes(std::ostream& os, const char* heading = 0) const;

  //! \brief Returns the global node numbers of this patch.
  const IntVec& getGlobalNodeNums() const { return MLGN; }
  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;
  //! \brief Returns the total number of elements in this patch.
  size_t getNoElms(bool includeZeroVolumeElms = false) const;
  //! \brief Returns the total number of MPC equations in this patch.
  size_t getNoMPCs() const { return mpcs.size(); }

  //! \brief Returns the beginning of the BC array.
  BCVec::const_iterator begin_BC() const { return BCode.begin(); }
  //! \brief Returns the end of the BC array.
  BCVec::const_iterator end_BC() const { return BCode.end(); }

  //! \brief Returns the beginning of the MNPC array.
  IntMat::const_iterator begin_elm() const { return MNPC.begin(); }
  //! \brief Returns the end of the MNPC array.
  IntMat::const_iterator end_elm() const { return MNPC.end(); }

  //! \brief Returns the beginning of the MPC set.
  MPCIter begin_MPC() const { return mpcs.begin(); }
  //! \brief Returns the end of the MPC set.
  MPCIter end_MPC() const { return mpcs.end(); }
  //! \brief Returns the MPC equation for a specified slave, if any.
  //! \param[in] node Global node number of the slave node
  //! \param[in] dof Which local DOF which is constrained (1, 2, 3)
  MPC* findMPC(int node, int dof) const;


  // Various preprocessing methods
  // =============================

  //! \brief Merges a given node in this patch with a given global node.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] globalNum Global number of the node to merge \a node with
  //! \param[in] silence If \e true, suppress message on the merged nodes
  bool mergeNodes(size_t inod, int globalNum, bool silence = false);

  //! \brief Renumbers all global node numbers in the entire model.
  //! \param[in] model All spline patches in the model
  //! \param[out] old2new Old-to-new node number mapping
  //! \return The number of unique nodes in the model
  //!
  //! \details After the renumbering, the global node numbers are in the range
  //! [1,\a nnod ], where \a nnod is the number of unique nodes in the model.
  //! The new node numbers computed by this methid preserve the relative
  //! ordering of the nodes. That is not the case when the non-static version
  //! is used.
  static int renumberNodes(const ASMVec& model, std::map<int,int>& old2new);

  //! \brief Renumbers the global node numbers in this patch.
  //! \param old2new Old-to-new node number mapping
  //! \param nnod Number of unique nodes found so far
  //! \return The number of renumbered nodes in this patch
  //!
  //! \details After the renumbering, the global node numbers are in the range
  //! [1,\a nnod ], where \a nnod is the number of unique nodes in the model.
  int renumberNodes(std::map<int,int>& old2new, int& nnod);

  //! \brief Renumbers the global node numbers referred by this patch.
  //! \param[in] old2new Old-to-new node number mapping
  //!
  //! \details The node numbers referred by boundary condition and
  //! multi-point constraint objects in the patch are renumbered.
  //! The nodes themselves are assumed already to be up to date.
  //! \param[in] old2new Old-to-new node number mapping
  bool renumberNodes(const std::map<int,int>& old2new);

  //! \brief Computes the set of all MPCs over the whole model.
  //! \param[in] model All spline patches in the model
  //! \param[out] allMPCs All multi-point constraint equations in the model
  //!
  //! \details The MPC equations are stored distributed over the patches,
  //! based on the slave DOF of the constraint. Therefore, constraints on the
  //! interface nodes between patches (to enforce higher order regularity) may
  //! be defined on both neighboring patches. This method wil also merge such
  //! multiply defined equations into a single equations.
  static void mergeAndGetAllMPCs(const ASMVec& model, MPCSet& allMPCs);

  //! \brief Resolves (possibly multi-level) chaining in MPC equations.
  //! \param[in] allMPCs All multi-point constraint equations in the model
  //!
  //! \details If a master DOF in one MPC (multi-point constraint) equation
  //! is specified as slave by another MPC, it is replaced by the master(s) of
  //! that other equation. Since an MPC-equation may couple nodes belonging to
  //! different patches, this method must have access to all patches.
  static void resolveMPCchains(const MPCSet& allMPCs);

  //! \brief Initializes the multi-point constraint coefficients.
  virtual bool initConstraints() { return true;}

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  //! \param[in] time Current time
  virtual bool updateDirichlet(const std::map<int,RealFunc*>& func,
			       const std::map<int,VecFunc*>& vfunc,
			       double time = 0.0);

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ) = 0;

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int*) {}

  //! \brief Computes the total number of integration points in this patch.
  virtual void getNoIntPoints(size_t& nPt) = 0;
  //! \brief Computes the number of boundary integration points in this patch.
  virtual void getNoBouPoints(size_t& nPt, char ldim, char lindx) = 0;

  //! \brief Generates element groups for multi-threading of interior integrals.
  virtual void generateThreadGroups(const Integrand&, bool = false) {}
  //! \brief Generates element groups for multi-threading of boundary integrals.
  virtual void generateThreadGroups(char, bool = false) {}


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
			 GlobalIntegral& glbInt, const TimeDomain& time) = 0;

  //! \brief Evaluates a boundary integral over a patch face/edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary face/edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time) = 0;

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lEdge Local index of the patch edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrateEdge(Integrand& integrand, int lEdge,
 			     GlobalIntegral& glbInt,
			     const TimeDomain& time) { return false; }


  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameters in range [0.0,1.0] of the point
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that matches the point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const = 0;

  //! \brief Creates a standard FE model of this patch for visualization.
  //! \param[out] grid The generated finite element grid
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int* npe) const = 0;

  //! \brief Extract the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  virtual bool getSolution(Matrix& sField, const Vector& locSol,
			   const IntVec& nodes) const;

  //! \brief Evaluates the primary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
			    const int* npe) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size() \a X \a gpar[2].size().
  //! Otherwise, we assume that it contains the \a u, \a v and \a w parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
			    const RealArray* gpar, bool regular = true) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] project Flag indicating result recovery method
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! If \a npe is NULL, the solution is recovered or evaluated at the Greville
  //! points and then projected onto the spline basis to obtain the control
  //! point values, which then are returned through \a sField.
  //! If \a npe is not NULL and \a project is defined, the solution is also
  //! projected onto the spline basis, and then evaluated at the \a npe points.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
			    const int* npe = 0, char project = '\0') const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! When \a regular is \e true, it is assumed that the parameter value array
  //! \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size() \a X \a gpar[2].size().
  //! Otherwise, we assume that it contains the \a u, \a v and \a w parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
			    const RealArray* gpar, bool regular = true) const;

  //! \brief Projects the secondary solution using a (discrete) global L2-fit.
  //! \param[out] sField Secondary solution field control point values
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e true, a continuous L2-projection is used
  virtual bool globalL2projection(Matrix& sField,
				  const IntegrandBase& integrand,
				  bool continuous = false) const;

  //! \brief Projects the secondary solution using a continuous global L2-fit.
  //! \param[out] sField Secondary solution field control point values
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //!
  //! \details This method uses the \a integrate interface to perform numerical
  //! integration of the projection matrices. It should use the same integration
  //! scheme as \a SIMbase::assembleSystem and is therefore suitable for
  //! problems using internal integration point buffers (with history-dependent
  //! data, etc.) and that must be traveresed in the same sequence each time.
  //! \note The implementation of this method is placed in GlbL2projector.C
  bool L2projection(Matrix& sField,
                    const IntegrandBase& integrand,
                    const TimeDomain& time);


  // Methods for result extraction
  // =============================

  //! \brief Extracts element results for this patch from a global vector.
  //! \param[in] globRes Global matrix of element results
  //! \param[out] elmRes Element results for this patch
  void extractElmRes(const Matrix& globRes, Matrix& elmRes) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] madof Global Matrix of Accumulated DOFs
  void extractNodeVec(const Vector& globVec, Vector& nodeVec,
                      const int* madof) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] nndof Number of DOFs per node (the default is \a nf)
  //! \param[in] basis Which basis to extract nodal values for (mixed methods)
  virtual void extractNodeVec(const Vector& globVec, Vector& nodeVec,
			      unsigned char nndof = 0, int basis = 1) const;

  //! \brief Injects nodal results for this patch into the global vector.
  //! \param[in] nodeVec Nodal result vector for this patch
  //! \param globVec Global solution vector in DOF-order
  //! \param[in] nndof Number of DOFs per node (the default is \a nf)
  virtual bool injectNodeVec(const Vector& nodeVec, Vector& globVec,
			     unsigned char nndof = 0) const;

protected:

  // Internal methods for preprocessing of boundary conditions
  // =========================================================

  //! \brief Checks wether a given DOF is fixed or not.
  //! \param[in] node Global node number of the DOF to check
  //! \param[in] dof Local index of the DOF to check
  bool isFixed(int node, int dof) const;
  //! \brief Adds a general multi-point-constraint (MPC) equation to this patch.
  //! \param mpc Pointer to an MPC-object
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  //! \param[in] silence If \e true, suppress debug print
  bool addMPC(MPC*& mpc, int code = 0, bool silence = false);
  //! \brief Creates and adds a two-point constraint to this patch.
  //! \param[in] slave Global node number of the node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  //! \param[in] master Global node number of the master node of the constraint
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  bool add2PC(int slave, int dir, int master, int code = 0);
  //! \brief Creates and adds a three-point constraint to this patch.
  //! \param[in] slave Global node number of the node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  //! \param[in] master1 Global node number of 1st master node of the constraint
  //! \param[in] master2 Global node number of 2nd master node of the constraint
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  bool add3PC(int slave, int dir, int master1, int master2, int code = 0);
  //! \brief Creates and adds a single-point constraint to this patch.
  //! \param[in] node Global node number of the node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  bool addSPC(int node, int dir, int code);
  //! \brief Creates and adds a periodicity constraint to this patch.
  //! \param[in] master 1-based local index of the master node
  //! \param[in] slave 1-based local index of the slave node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  bool addPeriodicity(size_t master, size_t slave, int dir);
  //! \brief Creates periodicity constraints between two nodes in this patch.
  //! \param[in] master 1-based local index of the master node
  //! \param[in] slave 1-based local index of the slave node to constrain
  //! \param[in] dirs Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  void makePeriodic(size_t master, size_t slave, int dirs = 123);

public:
  //! \brief Constrains DOFs in the given node to the given value.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] dirs Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  void prescribe(size_t inod, int dirs, int code);
  //! \brief Constrains DOFs in the given node to zero.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] dirs Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  void fix(size_t inod, int dirs = 123);

private:
  //! \brief Recursive method used by \a resolveMPCchains.
  //! \param[in] allMPCs All multi-point constraint equations in the model
  //! \param mpc Pointer to the multi-point constraint equation to resolve
  static bool resolveMPCchain(const MPCSet& allMPCs, MPC* mpc);

protected:
  //! \brief Collapses the given two nodes into one.
  //! \details The global node number of the node with the highest number
  //! is changed into the number of the other node.
  static bool collapseNodes(ASMbase& pch1, int node1, ASMbase& pch2, int node2);

public:
  static bool fixHomogeneousDirichlet; //!< If \e true, pre-eliminate fixed DOFs

  size_t idx; //!< Index of this patch in the multi-patch model

protected:
  // Standard finite element data structures
  unsigned char ndim;   //!< Number of parametric dimensions (1, 2 or 3)
  unsigned char nsd;    //!< Number of space dimensions (ndim <= nsd <= 3)
  unsigned char nf;     //!< Number of primary solution fields (1 or larger)
  unsigned char nLag;   //!< Number of Lagrange multipliers per node
  int           nGauss; //!< Numerical integration scheme

  const IntVec& MLGE; //!< Matrix of Local to Global Element numbers
  const IntVec& MLGN; //!< Matrix of Local to Global Node numbers
  const IntMat& MNPC; //!< Matrix of Nodal Point Correspondance
  const bool shareFE; //!< If \e true, this patch uses FE data of another patch

  BCVec  BCode; //!< Array of Boundary condition codes
  MPCMap dCode; //!< Inhomogeneous Dirichlet condition codes for the MPCs
  MPCSet mpcs;  //!< All multi-point constraints with the slave in this patch

  IntVec myMLGE; //!< The actual Matrix of Local to Global Element numbers
  IntVec myMLGN; //!< The actual Matrix of Local to Global Node numbers
  IntMat myMNPC; //!< The actual Matrix of Nodal Point Correspondance

  size_t nXelm;  //!< Number of extra-ordinary elements

private:
  std::pair<size_t,size_t> myLMs; //!< Nodal range of the Lagrange multipliers
};

#endif
