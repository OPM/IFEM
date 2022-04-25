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
#include <map>
#include <set>
#include <array>
#include <string>

typedef std::vector<int>       IntVec;  //!< General integer vector
typedef std::vector<IntVec>    IntMat;  //!< General 2D integer matrix
typedef std::map<MPC*,int,MPCLess> MPCMap; //!< MPC to function code mapping
typedef std::set<MPC*,MPCLess> MPCSet;  //!< Sorted set of MPC equations
typedef MPCSet::const_iterator MPCIter; //!< Iterator over an MPC equation set

struct TimeDomain;
class ElementBlock;
class Field;
class Fields;
class GlobalIntegral;
class IntegrandBase;
class Integrand;
class L2Integrand;
class ASMbase;
class SparseMatrix;
class StdVector;
class FunctionBase;
class RealFunc;
class VecFunc;
class Vec3;
class Tensor;
namespace ASM { class InterfaceChecker; }

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
    char RX;  //!< Boundary condition code for X-rotation
    char RY;  //!< Boundary condition code for Y-rotation
    char RZ;  //!< Boundary condition code for Z-rotation

    //! \brief Constructor initializing a BC instance.
    explicit BC(int n) : node(n), CX(1), CY(1), CZ(1), RX(1), RY(1), RZ(1) {}
  };

  typedef std::vector<BC> BCVec; //!< Nodal boundary condition container

protected:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMbase(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Special copy constructor for sharing of FE data.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  //!
  //! \details This copy constructor makes this patch sharing the FE data
  //! (including the basis functions) with the provided \a patch.
  //! The boundary conditions and constraint equations are however not copied,
  //! as this copy constructor is typically used for multi-stage simulators
  //! with different sub-problems discretized on the same grid.
  ASMbase(const ASMbase& patch, unsigned char n_f);
  //! \brief Default copy constructor, copying everything except \a neighbors.
  //! \param[in] patch The patch to copy
  ASMbase(const ASMbase& patch);

public:
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~ASMbase();

  //! \brief Returns a copy of this patch with identical FE discretization.
  //! \note The copied patch shares the spline data structures with the copy,
  //! in order to save memory. Thus, the copy cannot be read from file, refined,
  //! or changed in other ways that affect the FE geometry and/or topology.
  //! The other properties of the patch (FE topology, boundary conditions,
  //! constraints, loads, etc.) are however not copied.
  ASMbase* cloneUnShared() const;

  //! \brief Checks if this patch is empty.
  virtual bool empty() const = 0;

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is, int basis = 0) = 0;
  //! \brief Writes the geometry/basis of the patch to the given stream.
  virtual bool write(std::ostream& os, int basis = 0) const = 0;

  //! \brief Adds a circular immersed boundary in the physical geometry.
  virtual void addHole(double, double, double) {}
  //! \brief Adds an oval immersed boundary in the physical geometry.
  virtual void addHole(double, double, double, double, double) {}
  //! \brief Defines the immersed geometry from a scalar function.
  virtual bool setGeometry(RealFunc*, double, double) { return false; }

  //! \brief Generates the finite element topology data for this patch.
  virtual bool generateFEMTopology() = 0;

  //! \brief Clears the contents of this patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Adds extraordinary elements associated with a patch boundary.
  //! \param[in] dim Dimension of the boundary (should be \a nsd - 1)
  //! \param[in] item Local index of the boundary face/edge
  //! \param[in] nXn Number of extraordinary nodes
  //! \param[out] nodes Global numbers assigned to the extraordinary nodes
  virtual bool addXElms(short int dim, short int item, size_t nXn,
                        IntVec& nodes);

  //! \brief Adds a set of Lagrange multipliers to the specified element.
  //! \param[in] iel 1-based element index local to current patch
  //! \param[in] mGLag Global node numbers of the Lagrange multipliers
  //! \param[in] nnLag Number of Lagrange multipliers per node
  bool addLagrangeMultipliers(size_t iel, const IntVec& mGLag,
                              unsigned char nnLag = 1);

  //! \brief Adds global Lagrange multipliers to the system.
  //! \param[in] mGLag Global node numbers of the Lagrange multipliers
  //! \param[in] nnLag Number of Lagrange multipliers to add
  bool addGlobalLagrangeMultipliers(const IntVec& mGLag,
                                    unsigned char nnLag = 1);

  //! \brief Defines the numerical integration scheme \a nGauss in the patch.
  void setGauss(int ng) { nGauss = ng; }

  //! \brief Defines the number of solution fields \a nf in the patch.
  //! \details This method is to be used by simulators where \a nf is not known
  //! when the patch is constructed, e.g., it depends on the input file content.
  //! It must be invoked only before SIMbase::preprocess() is invoked.
  void setNoFields(unsigned char n) { nf = n; }

  //! \brief Sets the minimum element size for adaptive refinement.
  virtual void setMinimumSize(double) {}
  //! \brief Defines the minimum element size for adaptive refinement.
  virtual double getMinimumSize(int = 0) const { return 0.0; }
  //! \brief Checks if the specified element is larger than the minimum size.
  virtual bool checkElementSize(int, bool = true) const { return false; }

  //! \brief Resets the global element and node counters.
  static void resetNumbering(int n = 0);


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
  virtual bool getOrder(int&, int&, int&) const { return false; }

  //! \brief Returns local 1-based index of the node with given global number.
  //! \details If the given node number is not present, 0 is returned.
  //! \param[in] globalNum Global node number
  virtual size_t getNodeIndex(int globalNum, bool = false) const;
  //! \brief Returns the global node number for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual int getNodeID(size_t inod, bool = false) const;
  //! \brief Returns local 1-based index of element with given global number.
  //! \details If the given node number is not present, 0 is returned.
  //! \param[in] globalNum Global element number
  size_t getElmIndex(int globalNum) const;
  //! \brief Returns the global element number for the given element.
  //! \param[in] iel 1-based element index local to current patch
  int getElmID(size_t iel) const;
  //! \brief Returns the number of DOFs per node.
  //! \param[in] inod 1-based node index local to current patch
  virtual unsigned char getNodalDOFs(size_t inod) const;
  //! \brief Returns the classification of a node.
  //! \param[in] inod 1-based node index local to current patch
  virtual char getNodeType(size_t inod) const;
  //! \brief Returns \e true if node \a n is a Lagrange multiplier node.
  bool isLMn(size_t n) const { return myLMs.find(n) != myLMs.end(); }
  //! \brief Returns \e true if node \a n is a master node of a rigid coupling.
  bool isRMn(size_t n) const { return myRmaster.find(n) != myRmaster.end(); }
  //! \brief Returns the type of a Lagrange multiplier node.
  //! \param[in] inod 1-based node index for the Lagrange multiplier
  char getLMType(size_t inod) const;
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

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face/edge
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (for mixed methods)
  //! \param[in] thick Thickness of connection
  //! \param[in] orient Local orientation of the boundary face/edge
  //! \param[in] local If \e true, return patch-local numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int basis = 0, int thick = 1,
                                int orient = -1, bool local = false) const = 0;

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary edge/vertex
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (for mixed methods)
  //! \param[in] orient Local orientation flag (for LR splines only)
  //! \param[in] local If \e true, return patch-local numbers
  //! \param[in] open If \e true, exclude edge end points
  virtual void getBoundary1Nodes(int lIndex, IntVec& nodes,
                                 int basis = 0, int orient = -1,
                                 bool local = false, bool open = false) const {}

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face/edge
  //! \param[in] orient Local orientation of the boundary face/edge
  //! \param[out] elms Array of element numbers
  virtual void getBoundaryElms(int lIndex, int orient, IntVec& elms) const = 0;

  //! \brief Returns (1-based) index of a predefined node set in the patch.
  virtual int getNodeSetIdx(const std::string&) const { return 0; }
  //! \brief Returns an indexed predefined node set.
  virtual const IntVec& getNodeSet(int) const { static IntVec v; return v; }
  //! \brief Returns a named node set for update.
  virtual IntVec& getNodeSet(const std::string&, int&)
  { static IntVec v; return v; }

  //! \brief Returns (1-based) index of a predefined element set in the patch.
  virtual int getElementSetIdx(const std::string&) const { return 0; }
  //! \brief Returns an indexed predefined element set.
  virtual const IntVec& getElementSet(int) const { static IntVec v; return v; }
  //! \brief Returns a named element set for update.
  virtual IntVec& getElementSet(const std::string&, int&)
  { static IntVec v; return v; }

  //! \brief Finds the node that is closest to the given point.
  virtual std::pair<size_t,double> findClosestNode(const Vec3&) const
  { return std::make_pair(0,-1.0); }

  //! \brief Prints out the nodal coordinates of this patch to the given stream.
  void printNodes(std::ostream& os) const;
  //! \brief Prints out element connections of this patch to the given stream.
  void printElements(std::ostream& os) const;

  //! \brief Increase all global node numbers by \a nshift.
  virtual void shiftGlobalNodeNums(int nshift);
  //! \brief Sets the global node numbers for this patch.
  void setGlobalNodeNums(const IntVec& nodes) { myMLGN = nodes; }
  //! \brief Returns the actual global node numbers of this patch.
  const IntVec& getMyNodeNums() const { return myMLGN; }
  //! \brief Returns the global node numbers of this patch.
  const IntVec& getGlobalNodeNums() const { return MLGN; }

  //! \brief Increase all global element numbers by \a eshift.
  virtual void shiftGlobalElmNums(int eshift);
  //! \brief Returns the global element numbers of this patch.
  const IntVec& getGlobalElementNums() const { return MLGE; }

  //! \brief Returns the nodal point correspondance array for an element.
  //! \param[in] iel 1-based element index local to current patch
  const IntVec& getElementNodes(int iel) const;
  //! \brief Returns number of bases of this patch.
  virtual size_t getNoBasis() const { return 1; }
  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;
  //! \brief Returns the total number of elements in this patch.
  //! \param[in] includeZeroVolElms If \e true, count all the regular elements
  //! in the patch, including the zero-volume elements due to multiple knots
  //! \param[in] includeXElms If \e true, include any extra-ordinary elements
  //! in the patch (contact and interface elements, but no zero-volume elements)
  size_t getNoElms(bool includeZeroVolElms = false,
                   bool includeXElms = false) const;
  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char, char) const { return 0; }
  //! \brief Returns the total number of MPC equations in this patch.
  size_t getNoMPCs() const { return mpcs.size(); }

  //! \brief Computes the total number of integration points in this patch.
  virtual void getNoIntPoints(size_t& nPt, size_t& nIPt);
  //! \brief Computes the number of boundary integration points in this patch.
  virtual void getNoBouPoints(size_t& nPt, char ldim, char lindx);

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

  //! \brief Returns \e true if this patch shares FE data with another patch.
  bool isShared() const { return shareFE == 'F'; }
  //! \brief Returns \e true if this patch has additional (extraordinary) nodes.
  bool hasXNodes() const { return MLGN.size() > nnod; }

  //! \brief Returns parameter values and node numbers of the domain corners.
  //! \param[out] u Parameter values of the domain corners
  //! \param[out] corners 1-based indices of the corner nodes (optional)
  virtual bool getParameterDomain(Real2DMat& u,
                                  IntVec* corners = nullptr) const = 0;

  //! \brief Obtain element neighbours.
  virtual void getElmConnectivities(IntMat& neighs) const = 0;

  // Various preprocessing methods
  // =============================

  //! \brief Copies the parameter domain from another patch.
  virtual void copyParameterDomain(const ASMbase*) {}

  //! \brief Makes two opposite boundaries periodic.
  //! \param[in] dir Parameter direction defining the periodic boundaries
  //! \param[in] basis Which basis to connect (mixed methods)
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeBoundaries(int dir = 1, int basis = 0, int master = 1) {}

  //! \brief Merges a given node in this patch with a given global node.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] globalNum Global number of the node to merge \a node with
  //! \param[in] verbose If \e true, print message on the merged nodes
  bool mergeNodes(size_t inod, int globalNum, bool verbose = true);

  //! \brief Renumbers all global node numbers in the entire model.
  //! \param[in] model All spline patches in the model
  //! \param[out] old2new Old-to-new node number mapping
  //! \return The number of unique nodes in the model
  //!
  //! \details After the renumbering, the global node numbers are in the range
  //! [1,\a nNod ], where \a nNod is the number of unique nodes in the model.
  //! The new node numbers computed by this method preserve the relative
  //! ordering of the nodes. That is not the case when the non-static version
  //! is used.
  static int renumberNodes(const ASMVec& model, std::map<int,int>& old2new);

  //! \brief Renumbers the global node numbers in this patch.
  //! \param old2new Old-to-new node number mapping
  //! \param nNod Number of unique nodes found so far
  //! \return The number of renumbered nodes in this patch
  //!
  //! \details After the renumbering, the global node numbers are in the range
  //! [1,\a nNod ], where \a nNod is the number of unique nodes in the model.
  int renumberNodes(std::map<int,int>& old2new, int& nNod);

  //! \brief Renumbers the global node numbers referred by this patch.
  //! \param[in] old2new Old-to-new node number mapping
  //! \param[in] renumGN Flag for renumbering the node number array \a MLGN
  bool renumberNodes(const std::map<int,int>& old2new, char renumNodes = 0);

  //! \brief Computes the set of all MPCs over the whole model.
  //! \param[in] model All spline patches in the model
  //! \param[out] allMPCs All multi-point constraint equations in the model
  //!
  //! \details The MPC equations are stored distributed over the patches,
  //! based on the slave DOF of the constraint. Therefore, constraints on the
  //! interface nodes between patches (to enforce higher order regularity) may
  //! be defined on both neighboring patches. This method will also merge such
  //! multiply defined equations into a single equation.
  static void mergeAndGetAllMPCs(const ASMVec& model, MPCSet& allMPCs);

  //! \brief Resolves (possibly multi-level) chaining in MPC equations.
  //! \param[in] allMPCs All multi-point constraint equations in the model
  //! \param[in] model All spline patches in the model
  //! \param[in] setPtrOnly If \e true, only set pointer to next MPC in chain
  static void resolveMPCchains(const MPCSet& allMPCs,
                               const ASMVec& model, bool setPtrOnly = false);

  //! \brief Initializes the multi-point constraint coefficients.
  virtual bool initConstraints() { return true; }

  //! \brief Sets the global node numbers for this patch.
  //! \param[in] nodes Vector of zero-based node numbers of this patch
  virtual void setNodeNumbers(const IntVec& nodes);

  //! \brief Checks for time-dependent in-homogeneous Dirichlet conditions.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  bool hasTimeDependentDirichlet(const std::map<int,RealFunc*>& func,
                                 const std::map<int,VecFunc*>& vfunc);

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  //! \param[in] time Current time
  //! \param[in] g2l Global-to-local mapping to apply to node numbers
  virtual bool updateDirichlet(const std::map<int,RealFunc*>& func,
                               const std::map<int,VecFunc*>& vfunc,
                               double time = 0.0,
                               const std::map<int,int>* g2l = nullptr);

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ) = 0;

  //! \brief Applies a transformation matrix from local to global system.
  virtual bool transform(const Matrix&) { return true; }

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int*) {}

  //! \brief Generates element groups for multi-threading of interior integrals.
  virtual void generateThreadGroups(const Integrand&, bool, bool) {}
  //! \brief Generates element groups for multi-threading of boundary integrals.
  virtual void generateThreadGroups(char, bool, bool) {}
  //! \brief Generate element-groups for multi-threading based on a partition.
  virtual void generateThreadGroupsFromElms(const IntVec&) {}


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

  //! \brief Evaluates an integral over element interfaces in the patch.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] iChk Object checking if an element interface has contributions
  virtual bool integrate(Integrand& integrand, GlobalIntegral& glbInt,
                         const TimeDomain& time,
                         const ASM::InterfaceChecker& iChk) { return false; }

  //! \brief Integrates a spatial dirac-delta function over a patch.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] u Parameters of the non-zero point of the dirac-delta function
  //! \param[in] p Function value at the specified point
  virtual bool diracPoint(Integrand& integrand, GlobalIntegral& glbInt,
                          const double* u, const Vec3& p) { return false; }


  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameters in range [0.0,1.0] of the point
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that matches the point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const = 0;
  //! \brief Returns the element that contains a specified spatial point.
  //! \param[in] param The parameters of the point in the knot-span domain
  //! \return Local element number within the patch that contains the point
  virtual int findElementContaining(const double* param) const = 0;

  //! \brief Creates a standard FE model of this patch for visualization.
  //! \param[out] grid The generated finite element grid
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int* npe) const = 0;
  //! \brief Returns an additional geometry to visualize (immersed boundaries).
  virtual ElementBlock* immersedGeometry() const { return nullptr; }
  //! \brief Filters out result point values that are outside physical domain.
  virtual void filterResults(Matrix&, const ElementBlock*) const {}

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
  //! \param[in] nf If nonzero, mixed evaluates nf fields on first basis
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int* npe, int nf = 0) const;

  //! \brief Evaluates the projected solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] nf If nonzero, mixed evaluates nf fields on first basis
  virtual bool evalProjSolution(Matrix& sField, const Vector& locSol,
                                const int* npe, int nf = 0) const
  { return this->evalSolution(sField, locSol, npe, nf); }

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //! \param[in] deriv Derivative order to return
  //! \param[in] nf If non-zero, mixed evaluates \a nf fields on first basis
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size() \a X \a gpar[2].size().
  //! Otherwise, we assume that it contains the \a u, \a v and \a w parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool regular = true,
                            int deriv = 0, int nf = 0) const;

  //! \brief Evaluates and interpolates a field over a given geometry.
  //! \param[in] basis The basis of the field to evaluate
  //! \param[in] locVec The coefficients of the field to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum The basis to evaluate for (mixed)
  virtual bool evaluate(const ASMbase* basis, const Vector& locVec,
                        RealArray& vec, int basisNum = 1) const;

  //! \brief Evaluates and interpolates a scalar field over a given geometry.
  //! \param[in] f The field to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum The basis to evaluate for (mixed)
  virtual bool evaluate(const Field* f, RealArray& vec, int basisNum = 1) const;

  //! \brief Evaluates and interpolates a scalar function over a given geometry.
  //! \param[in] f The function to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum The basis to evaluate for (mixed)
  //! \param[in] time Current time
  virtual bool evaluate(const FunctionBase* f, RealArray& vec,
                        int basisNum = 1, double time = 0.0) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] project Flag indicating result recovery method
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! If \a npe is null, the solution is recovered or evaluated at the Greville
  //! points and then projected onto the spline basis to obtain the control
  //! point values, which then are returned through \a sField.
  //! If \a npe is not null and \a project is defined, the solution is also
  //! projected onto the spline basis, and then evaluated at the \a npe points.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const int* npe = nullptr, char project = 0) const;

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
                                  const L2Integrand& integrand,
                                  bool continuous = false) const;

  //! \brief Projects the secondary solution using a continuous global L2-fit.
  //! \param[out] sField Secondary solution field control point values
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //!
  //! \details This method uses the integrate() interface to perform numerical
  //! integration of the projection matrices. It should use the same integration
  //! scheme as SIMbase::assembleSystem() and is therefore suitable for
  //! problems using internal integration point buffers (with history-dependent
  //! data, etc.) and that must be traversed in the same sequence each time.
  //! \note The implementation of this method is placed in GlbL2projector.C
  bool L2projection(Matrix& sField, IntegrandBase* integrand,
                    const TimeDomain& time);

  //! \brief Projects an explicit function using a continuous global L2-fit.
  //! \param[out] fVals Control point values of the function
  //! \param[in] function The function to project
  //! \param[in] t Current time
  //!
  //! \note The implementation of this method is placed in GlbL2projector.C
  bool L2projection(Matrix& fVals, FunctionBase* function, double t = 0.0);
  //! \brief Projects explicit functions using a continuous global L2-fit.
  //! \param[out] fVals Control point values of the functions
  //! \param[in] function The functions to project
  //! \param[in] t Current time
  //!
  //! \note The implementation of this method is placed in GlbL2projector.C
  bool L2projection(const std::vector<Matrix*>& fVals,
                    const std::vector<FunctionBase*>& function, double t = 0.0);

  //! \brief Returns the number of projection nodes for this patch.
  virtual size_t getNoProjectionNodes() const { return this->getNoNodes(1); }

  //! \brief Returns the number of nodes on refinement basis for this patch.
  virtual size_t getNoRefineNodes() const { return this->getNoNodes(1); }

  //! \brief Returns the number of elements on refinement basis for this patch.
  virtual size_t getNoRefineElms() const { return this->getNoElms(); }

  //! \brief Returns a field using the projection basis.
  virtual Field* getProjectedField(const Vector&) const
  { return nullptr; }

  //! \brief Returns a field using the projection basis.
  virtual Fields* getProjectedFields(const Vector&, size_t) const
  { return nullptr; }

  //! \brief Creates a separate projection basis for this patch.
  virtual bool createProjectionBasis(bool) { return false; }
  //! \brief Checks if a separate projection basis is used for this patch.
  virtual bool separateProjectionBasis() const { return false; }


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
  //! \param[in] ngnod Dimension of madof (the default -1 means unknown)
  bool extractNodalVec(const RealArray& globVec, RealArray& nodeVec,
                       const int* madof, int ngnod = -1) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] nndof Number of DOFs per node (the default is \a nf)
  //! \param[in] basis Which basis to extract nodal values for (mixed methods)
  virtual void extractNodeVec(const RealArray& globVec, RealArray& nodeVec,
                              unsigned char nndof = 0, int basis = 0) const;

  //! \brief Injects nodal results for this patch into the global vector.
  //! \param[in] nodeVec Nodal result vector for this patch
  //! \param globVec Global solution vector in DOF-order
  //! \param[in] nndof Number of DOFs per node (the default is \a nf)
  //! \param[in] basis Which basis to inject nodal values for (mixed methods)
  virtual bool injectNodeVec(const RealArray& nodeVec, RealArray& globVec,
                             unsigned char nndof = 0, int basis = 0) const;

  //! \brief Injects nodal results for this patch into the global vector.
  //! \param[in] nodeVec Nodal result vector for this patch
  //! \param globVec Global solution vector in DOF-order
  //! \param[in] madof Global Matrix of Accumulated DOFs
  //! \param[in] basis Which basis to inject nodal values for (mixed methods)
  bool injectNodalVec(const RealArray& nodeVec, RealArray& globVec,
                      const IntVec& madof, int basis = 0) const;

  //! \brief Creates and adds a two-point constraint to this patch.
  //! \param[in] slave Global node number of the node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  //! \param[in] master Global node number of the master node of the constraint
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  bool add2PC(int slave, int dir, int master, int code = 0);

  //! \brief Adds a general multi-point-constraint (MPC) equation to this patch.
  //! \param mpc Pointer to an MPC-object
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  //! \param[in] verbose If \e true, print out added constraint (debug build)
  bool addMPC(MPC*& mpc, int code = 0, bool verbose = false);

  //! \brief Adds MPCs representing a rigid coupling to this patch.
  //! \param[in] lindx Local index of the boundary item that should be rigid
  //! \param[in] ldim Dimension of the boundary item that should be rigid
  //! \param[in] basis Which basis to add rigid coupling for (mixed methods)
  //! \param gMaster Global node number of the master node
  //! \param[in] Xmaster Position of the master nodal point
  //! \param[in] extraPt If \e true, the master point is not a patch node
  //! \return \e true if a new global node was added, otherwise \e false
  virtual bool addRigidCpl(int lindx, int ldim, int basis,
                           int& gMaster, const Vec3& Xmaster,
                           bool extraPt = true);

protected:

  // Internal methods for preprocessing of boundary conditions
  // =========================================================

  //! \brief Creates constraint equations coupling global DOFs to local DOFs.
  //! \param[in] iSlave 0-based local index of node with global DOFs
  //! \param[in] master Global node number of node with local DOFs
  //! \param[in] Tlg Local-to-global transformation matrix
  void addLocal2GlobalCpl(int iSlave, int master, const Tensor& Tlg);
  //! \brief Creates and adds a three-point constraint to this patch.
  //! \param[in] slave Global node number of the node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  //! \param[in] master1 Global node number of 1st master node of the constraint
  //! \param[in] master2 Global node number of 2nd master node of the constraint
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  bool add3PC(int slave, int dir, int master1, int master2, int code = 0);
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
  //! \brief Adds a patch to the list of neighbors of this patch.
  //! \param[in] pch Pointer to the neighboring patch
  void addNeighbor(ASMbase* pch);
  //! \brief Creates an additional master node for a rigid coupling.
  //! \param gMaster Global node number of the master node
  //! \param[in] Xmaster Position of the master nodal point
  //! \return \e true if a new global node was added, otherwise \e false
  bool createRgdMasterNode(int& gMaster, const Vec3& Xmaster);
  //! \brief Adds MPC equations representing a rigid arm to a 6-DOF node
  //! \param[in] gSlave Global node number of the 3-DOF slave node
  //! \param[in] gMaster Global node number of the 6-DOF master node
  //! \param[in] dX Relative position of the slave w.r.t. the master node
  void addRigidMPC(int gSlave, int gMaster, const Vec3& dX);

  //! \brief Returns the number of Gauss points to use in one direction.
  //! \param[in] p Polynomial order of the basis functions
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  int getNoGaussPt(int p, bool neumann = false) const;


  // Miscellaneous methods for internal use
  // ======================================

  //! \brief Helper method used by evalPoint to search for a control point.
  //! \param[in] cit iterator of array of control point coordinates
  //! \param[in] end iterator of array of control point coordinates
  //! \param[in] X Coordinates of point to search for
  //! \param[in] dimension Number of spatial dimensions of the splines object
  //! \param[in] tol Zero tolerance
  int searchCtrlPt(RealArray::const_iterator cit, RealArray::const_iterator end,
                   const Vec3& X, int dimension, double tol = 0.001) const;

  //! \brief Assembles L2-projection matrices for the secondary solution.
  //! \param[out] A Left-hand-side matrix
  //! \param[out] B Right-hand-side vectors
  //! \param[in] obj Wrapper object for integrand/function
  //! \param[in] continuous If \e false, a discrete L2-projection is used
  virtual bool assembleL2matrices(SparseMatrix& A, StdVector& B,
                                  const L2Integrand& obj,
                                  bool continuous) const { return false; }

public:

  // More methods for preprocessing of Dirichlet boundary conditions
  // ===============================================================

  //! \brief Constrains all nodes in the patch.
  //! \param[in] dof Which DOFs to constrain at each node in the patch
  //! \param[in] code Inhomogeneous dirichlet condition code
  void constrainPatch(int dof, int code = 0);
  //! \brief Constrains a list of nodes in the patch.
  //! \param[in] nodes 1-based list of nodes to constrain
  //! \param[in] dof Which DOFs to constrain at each node
  //! \param[in] code Inhomogeneous dirichlet condition code
  void constrainNodes(const IntVec& nodes, int dof, int code = 0);
  //! \brief Constrains an extraordinary node in the patch.
  //! \param[in] node Global node number of the node to constrain
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \return \e true if \a node is a rigid master point in the patch,
  //! otherwise \e false
  bool constrainXnode(int node, int dof, int code = 0);
  //! \brief Constrains DOFs in the given node to the given value.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] dirs Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  //! \return Invalid local DOFs, or DOFs that already are constrained
  int prescribe(size_t inod, int dirs, int code);
  //! \brief Constrains DOFs in the given node to zero.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] dirs Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  //! \return Invalid local DOFs
  int fix(size_t inod, int dirs = 123);
  //! \brief Checks if the given DOFs are fixed.
  //! \param[in] node Global node number of the DOF to check
  //! \param[in] dof Local indices of the DOFs to check
  //! \param[in] all Returns \e true only if all DOFs are fixed
  bool isFixed(int node, int dof, bool all = false) const;

protected:
  //! \brief Returns \e true if \a dirs constains all local DOFs in the patch.
  bool allDofs(int dirs) const;

  //! \brief Calculated the deformed configuration for current element.
  //! \param[in] Xnod Array of nodal point coordinates for current element
  //! \param eVec Current element solution vectors
  //! \param[in] force2nd If \e true, put updated coordinates as the 2nd vector
  bool deformedConfig(const RealArray& Xnod, Vectors& eVec,
                      bool force2nd = false) const;

  //! \brief Collapses the given two nodes into one.
  //! \details The global node number of the node with the highest number
  //! is changed into the number of the other node.
  static bool collapseNodes(ASMbase& pch1, int node1, ASMbase& pch2, int node2);

  //! \brief Writes a Lagrangian basis to the given stream.
  bool writeLagBasis(std::ostream& os, const char* type) const;

public:
  static bool fixHomogeneousDirichlet; //!< If \e true, pre-eliminate fixed DOFs

  static int dbgElm; //!< One-based element index to print debugging info for

  static double modelSize; //!< Characteristic model size

  size_t idx; //!< Index of this patch in the multi-patch model

protected:
  // Standard finite element data structures
  unsigned char ndim;   //!< Number of parametric dimensions (1, 2 or 3)
  unsigned char nsd;    //!< Number of space dimensions (ndim <= nsd <= 3)
  unsigned char nf;     //!< Number of primary solution fields (1 or larger)
  unsigned char nLag;   //!< Number of Lagrange multipliers per node
  size_t        nel;    //!< Number of regular elements in this patch
  size_t        nnod;   //!< Number of regular nodes in this patch

  const IntVec& MLGE; //!< Matrix of Local to Global Element numbers
  const IntVec& MLGN; //!< Matrix of Local to Global Node numbers
  const IntMat& MNPC; //!< Matrix of Nodal Point Correspondance

  //! \brief Flag telling whether this patch shares its data with another patch.
  //! \details 'S' means this patch uses spline geometry of another patch.
  //! 'F' means this patch uses FE data and spline geometry of another patch.
  const char shareFE; //!< If \e true, this patch uses FE data of another patch

  BCVec  BCode; //!< Array of Boundary condition codes
  MPCMap dCode; //!< Inhomogeneous Dirichlet condition codes for the MPCs
  MPCSet mpcs;  //!< All multi-point constraints with the slave in this patch

  IntVec myMLGE; //!< The actual Matrix of Local to Global Element numbers
  IntVec myMLGN; //!< The actual Matrix of Local to Global Node numbers
  IntMat myMNPC; //!< The actual Matrix of Nodal Point Correspondance
  IntVec myElms; //!< Elements on patch - used with partitioning

  //! \brief Numerical integration scheme for this patch.
  //! \details A value in the range [1,10] means use that number of Gauss
  //! quadrature points in each parameter direction, regardless of polynomial
  //! order of the basis functions. If zero or negative, the number of Gauss
  //! quadrature points is set independently in each parameter direction to
  //! \a p+nGauss, where \a p is the polynomial order in that direction.
  //! If the value is set larger than 10, the number of quadrature points
  //! in each parameter direction is set to \a p+nGauss%10.
  int nGauss; //!< \sa getNoGaussPt

  size_t firstIp; //!< Global index to first interior integration point

  //! Global indices to first integration point for the Neumann boundaries
  std::map<char,size_t> firstBp;

  ASMVec neighbors; //!< Patches having nodes in common with this one

  //! Auxilliary node number map used when establishing Dirichlet constraints
  static std::map<int,int> xNode;

  static int gEl;  //!< Global element counter
  static int gNod; //!< Global node counter

private:
  std::vector<char> myLMTypes; //!< Type of Lagrange multiplier ('L' or 'G')
  std::set<size_t>  myLMs;     //!< Nodal indices of the Lagrange multipliers

protected:
  typedef std::array<double,3> XYZ; //!< Convenience type definition
  std::map<size_t,XYZ>   myRmaster; //!< Rigid master nodal points
};

#endif
