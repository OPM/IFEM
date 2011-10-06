// $Id$
//==============================================================================
//!
//! \file ASMbase.h
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for spline-based FE assembly drivers.
//!
//==============================================================================

#ifndef _ASM_BASE_H
#define _ASM_BASE_H

#include "MatVec.h"
#include "MPCLess.h"
#include "Function.h"
#include <vector>
#include <set>
#include <map>

typedef std::vector<int>       IntVec;  //!< General integer vector
typedef std::vector<IntVec>    IntMat;  //!< General 2D integer matrix
typedef std::set<MPC*,MPCLess> MPCSet;  //!< Sorted set of MPC equations
typedef std::map<MPC*,int>     MPCMap;  //!< MPC to function code mapping
typedef MPCSet::const_iterator MPCIter; //!< Iterator over an MPC equation set

struct TimeDomain;
struct ElementBlock;
class GlobalIntegral;
class LocalIntegral;
class Integrand;
class ASMbase;
class Vec3;

typedef std::vector<LocalIntegral*> LintegralVec; //!< Local integral container
typedef std::vector<ASMbase*>       ASMVec;       //!< Spline patch container


/*!
  \brief Base class for spline-based finite element assembly drivers.

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
    BC(int n, char x, char y, char z) : node(n), CX(x), CY(y), CZ(z) {}
  };

  typedef std::vector<BC> BCVec; //!< Nodal boundary condition container

protected:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMbase(unsigned char n_p, unsigned char n_s, unsigned char n_f);

public:
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~ASMbase();

  //! \brief Checks if this patch is empty.
  virtual bool empty() const = 0;

  //! \brief Generates the finite element topology data for this patch.
  virtual bool generateFEMTopology() = 0;

  //! \brief Clears the contents of this patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Writes the geometry/basis of the patch to the given stream.
  virtual bool write(std::ostream&, int = 0) const { return false; }


  // Service methods for query of various model data
  // ===============================================

  //! \brief Returns the number of spatial dimensions.
  unsigned char getNoSpaceDim() const { return nsd; }
  //! \brief Returns the number of parameter dimensions.
  unsigned char getNoParamDim() const { return ndim; }
  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int b = 0) const { return b > 1 ? 0 : nf; }

  //! \brief Returns local 1-based index of the node with given global number.
  //! \details If the given node number is not present, 0 is returned.
  //! \param[in] globalNum Global node number
  size_t getNodeIndex(int globalNum) const;
  //! \brief Returns the global node number for the given node.
  //! \param[in] inod 1-based node index local to current patch
  int getNodeID(size_t inod) const;
  //! \brief Returns the global element number for the given element
  //! \param[in] iel 1-based element index local to current patch
  int getElmID(size_t iel) const;
  //! \brief Returns the number of DOFs per node.
  virtual unsigned char getNodalDOFs(size_t) const { return nf; }
  //! \brief Returns which mixed field basis a node belongs to.
  virtual unsigned char getNodalBasis(size_t) const { return 0; }
  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const = 0;
  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X nsd\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const = 0;
  //! \brief Prints out the nodal coordinates of this patch to the given stream.
  void printNodes(std::ostream& os, const char* heading = 0) const;

  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int = 0) const { return MLGN.size(); }
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


  // Various preprocessing methods
  // =============================

  //! \brief Merges a given node in this patch with a given global node.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] globalNum Global number of the node to merge \a node with
  bool mergeNodes(size_t inod, int globalNum);

  //! \brief Renumbers the global node numbers referred by this patch.
  //! \details The node numbers referred by boundary condition and
  //! multi-point constraint objects in the patch are renumbered.
  //! The nodes themselves are assumed already to be up to date.
  //! \param[in] old2new Old-to-new node number mapping
  //! \param[in] silent If \e false, give error on missing node in \a old2new
  bool renumberNodes(const std::map<int,int>& old2new, bool silent);

  //! \brief Renumbers all global node numbers in the entire model.
  //! \param[in] model All spline patches in the model
  //! \param[out] l2gn Local-to-global node numbers (optional)
  //! \return The number of unique nodes in the model
  //!
  //! \details After the renumbering, the global node numbers are in the range
  //! [1,\a nnod ], where \a nnod is the number of unique nodes in the model.
  static int renumberNodes(const ASMVec& model, IntVec* l2gn = 0);

  //! \brief Resolves (possibly multi-level) chaining in MPC equations.
  //! \param[in] model All spline patches in the model
  //!
  //! \details If a master DOF in one MPC is specified as slave by another MPC,
  //! it is replaced by the master(s) of that other equation.
  //! Since an MPC-equation may couple nodes belonging to different patches,
  //! this method must have all patches in the model available.
  static void resolveMPCchains(const ASMVec& model);

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] time Current time
  bool updateDirichlet(const std::map<int,RealFunc*>& func, double time = 0.0);

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int*) {}


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec()) = 0;

  //! \brief Evaluates a boundary integral over a patch face/edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary face/edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec()) = 0;

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
  virtual bool tesselate(ElementBlock& grid, const int* npe) const;

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
  //! \param[in] project Flag indicating if the projected solution is wanted
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! If \a npe is NULL, the solution is evaluated at the Greville points and
  //! then projected onto the spline basis to obtain the control point values,
  //! which then are returned through \a sField.
  //! If \a npe is not NULL and \a project is \e true, the solution is also
  //! projected onto the spline basis, and then evaluated at the \a npe points.
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const int* npe = 0, bool project = false) const;

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
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const RealArray* gpar, bool regular = true) const;


  // Methods for result extraction
  // =============================

  //! \brief Extracts element results for this patch from a global vector.
  //! \param[in] globRes Global matrix of element results
  //! \param[out] elmRes Element results for this patch
  void extractElmRes(const Matrix& globRes, Matrix& elmRes) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] nndof Number of DOFs per node (the default is \a nf)
  //! \param[in] basis Which basis to extract nodal values for (mixed methods)
  virtual void extractNodeVec(const Vector& globVec, Vector& nodeVec,
			      unsigned char nndof = 0, int basis = 0) const;

  //! \brief Injects nodal results for this patch into the global vector.
  //! \param[in] nodeVec Nodal result vector for this patch
  //! \param globVec Global solution vector in DOF-order
  //! \param[in] nndof Number of DOFs per node (the default is \a nf)
  virtual bool injectNodeVec(const Vector& nodeVec, Vector& globVec,
			     unsigned char nndof = 0) const;

protected:

  // Internal methods for preprocessing of boundary conditions
  // =========================================================

  //! \brief Adds a general multi-point-constraint equation to this patch.
  //! \param[in] mpc Pointer to an MPC-object
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  bool addMPC(MPC* mpc, int code = 0);
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
  //! \brief Constrains DOFs in the given node to the given value.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] dirs Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  //! \param[in] code Identifier for inhomogeneous Dirichlet condition field
  void prescribe(size_t inod, int dirs, int code);
  //! \brief Constrains DOFs in the given node to zero.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] dirs Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  void fix(size_t inod, int dirs = 123);

  //! \brief Recursive method used by \a resolveMPCchains.
  //! \param[in] allMPCs All multi-point constraint equations in the model
  //! \param mpc Pointer to the multi-point constraint equation to resolve
  static bool resolveMPCchain(const MPCSet& allMPCs, MPC* mpc);

  //! \brief Collapses the given two nodes into one.
  //! \details The global node number of the node with the highest number
  //! is changed into the number of the other node.
  static void collapseNodes(int& node1, int& node2);

public:
  static bool fixHomogeneousDirichlet; //!< If \e true, pre-eliminate fixed DOFs

protected:
  // Standard finite element data structures
  unsigned char ndim; //!< Number of parametric dimensions (1, 2 or 3)
  unsigned char nsd;  //!< Number of space dimensions (ndim <= nsd <= 3)
  unsigned char nf;   //!< Number of primary solution fields (1 or larger)

  IntVec MLGE;  //!< Matrix of Local to Global Element numbers
  IntVec MLGN;  //!< Matrix of Local to Global Node numbers
  IntMat MNPC;  //!< Matrix of Nodal Point Correspondance
  BCVec  BCode; //!< Vector of Boundary condition codes
  MPCSet mpcs;  //!< All multi-point constraints with slave in this patch
  MPCMap dCode; //!< Inhomogeneous Dirichlet condition codes for the MPCs
};

#endif
