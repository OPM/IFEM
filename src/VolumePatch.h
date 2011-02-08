// $Id: VolumePatch.h,v 1.22 2010-06-22 11:38:05 kmo Exp $
//==============================================================================
//!
//! \file VolumePatch.h
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of a topological cube as a SplineVolume.
//! \details The class contains the necessary data and methods to establish
//! a finite element algebraic system of equations for the stiffness relation
//! of the volume, using NURBS as basis functions via the GoTools library.
//!
//==============================================================================

#ifndef _VOLUME_PATCH_H
#define _VOLUME_PATCH_H

#include "Function.h"
#include "MPCLess.h"
#include "MatVec.h"
#include <set>
#include <map>

namespace Go {
  class SplineVolume;
}

struct ElementBlock;
class SystemVector;
class LinEqSystem;
class LocalSystem;
class SAM;
class Vec3;
class Tensor;

typedef std::vector<int>        IntVec;  //!< General integer vector
typedef std::vector<IntVec>     IntMat;  //!< General 2D integer matrix
typedef std::set<MPC*,MPCLess>  MPCSet;  //!< Sorted set of MPC equations
typedef MPCSet::const_iterator  MPCIter; //!< Iterator over a MPC equation set


/*!
  \brief Assembly of finite element contributions from a SplineVolume object.
  \details This class stores data and provides methods for the calculation of
  element stiffness and mass matrices, as well as load vectors for a structured
  isogeometric grid, using NURBS or splines as basis functions.
*/

class VolumePatch
{
  //! \brief Struct for nodal point data.
  struct IJK
  {
    int I; //!< Index in first parameter direction
    int J; //!< Index in second parameter direction
    int K; //!< Index in third parameter direction
    int global; //!< Global node number [1,nnod]
  };

public:
  //! \brief Struct for boundary condition codes.
  struct BC
  {
    int node; //!< Global node number [1,nnod]
    char CX; //!< Boundary condition code for X-translation
    char CY; //!< Boundary condition code for Y-translation
    char CZ; //!< Boundary condition code for Z-translation
    //! \brief Constructor initializing a BC instance.
    BC(int n, char x, char y, char z) : node(n), CX(x), CY(y), CZ(z) {}
  };

  //! \brief Constructor creating an instance by reading the given file.
  VolumePatch(const char* fileName, bool checkRHS = false);
  //! \brief Constructor creating an instance by reading the given input stream.
  VolumePatch(std::istream& is, bool checkRHS = false);
  //! \brief Default constructor creating an empty patch.
  VolumePatch() { svol = 0; swapW = false; E = nu = rho = 0.0; }
  //! \brief The Destructor frees the dynamically allocated SplineVolume object.
  ~VolumePatch();

  //! \brief Creates an instance by reading the given input stream, \a is.
  bool read(std::istream& is);
  //! \brief Writes the geometry of the SplineVolume object to the given stream.
  bool write(std::ostream& os) const;

  //! \brief Check that the patch is modelled in a right-hand-side system.
  //! \details If it isn't, the w-parameter direction is swapped.
  bool checkRightHandSystem();

  //! \brief Refine the parametrization by inserting extra knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  bool uniformRefine(int dir, int nInsert);
  //! \brief Raise the order of the SplineVolume object for this patch.
  //! \param[in] ru Number of times to raise the order i u-direction
  //! \param[in] rv Number of times to raise the order i u-direction
  //! \param[in] rw Number of times to raise the order i u-direction
  bool raiseOrder(int ru, int rv, int rw);

  //! \brief Generates the finite element topology data for the patch.
  bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  void clear();

  //! \brief Checks if the patch is empty.
  bool empty() const { return svol == 0; }

  //! \brief Returns local 1-based index of the node with given global number.
  //! \details If the given node number is not present, 0 is returned.
  //! \param[in] globalNum Global node number
  int getNodeIndex(int globalNum) const;
  //! \brief Returns the global node number for the given node.
  //! \param[in] inod 1-based node index local to current patch
  int getNodeID(int inod) const;
  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  Vec3 getCoord(int inod) const;

  //! \brief Returns the total number of nodes in the patch.
  size_t getNoNodes() const { return nodeInd.size(); }
  //! \brief Returns the total number of elements in the patch.
  size_t getNoElms() const { return MNPC.size(); }
  //! \brief Returns the total number of MPC equations in the patch.
  size_t getNoMPCs() const { return mpcs.size(); }

  //! \brief Returns the beginning of the BC array.
  std::vector<BC>::const_iterator begin_BC() const { return BCode.begin(); }
  //! \brief Returns the end of the BC array.
  std::vector<BC>::const_iterator end_BC() const { return BCode.end(); }

  //! \brief Returns the beginning of the MNPC array.
  IntMat::const_iterator begin_elm() const { return MNPC.begin(); }
  //! \brief Returns the end of the MNPC array.
  IntMat::const_iterator end_elm() const { return MNPC.end(); }

  //! \brief Returns the beginning of the MPC set.
  MPCIter begin_MPC() const { return mpcs.begin(); }
  //! \brief Returns the end of the MPC set.
  MPCIter end_MPC() const { return mpcs.end(); }

  //! \brief Merges a given node in this patch with a given global node.
  //! \param[in] node 1-based node index local to current patch
  //! \param[in] globalNum Global number of the node to merge \a node with
  bool mergeNodes(int node, int globalNum);

  //! \brief Renumbers the global node numbers referred by this patch.
  //! \details The node numbers referred by boundary condition and
  //! multi-point constraint objects in the patch are renumbered.
  //! The noded themselves are assumed to already be up to date.
  //! \param[in] old2new Old-to-new node number mapping
  //! \param[in] silent If \e false, give error on missing node in \a old2new
  bool renumberNodes(const std::map<int,int>& old2new, bool silent);

  //! \brief Renumbers all global node numbers in the entire model.
  //! \param[in] model All volume patches in the model
  //! \return The number of unique nodes in the model
  //!
  //! \details After the renumbering, the global node numbers are in the range
  //! [1,\a nnod ], where \a nnod is the number of unique nodes in the model.
  static int renumberNodes(const std::vector<VolumePatch*>& model);

  //! \brief Resolves (possibly multi-level) chaining in MPC equations.
  //! \param[in] model All volume patches in the model
  //!
  //! \details If a master DOF in one MPC is specified as slave by another MPC,
  //! it is replaced by the master(s) of that other equation.
  //! Since an MPC-equation may couple nodes belonging to different patches,
  //! this method must have all patches in the model available.
  static void resolveMPCchains(const std::vector<VolumePatch*>& model);

  //! \brief Defines material properties for current volume patch.
  //! \param[in] Emod    Young's modulus
  //! \param[in] Poiss   Poisson's ratio
  //! \param[in] Density Mass density
  void setMaterial(double Emod, double Poiss, double Density)
  { E = Emod; nu = Poiss; rho = Density; }
  //! \brief Retrieves material properties for current volume patch.
  //! \param[out] Emod    Young's modulus
  //! \param[out] Poiss   Poisson's ratio
  //! \param[out] Density Mass density
  void getMaterial(double& Emod, double& Poiss, double& Density)
  { Emod = E; Poiss = nu; Density = rho; }

  //! \brief Makes two opposite boundary faces periodic.
  //! \param[in] dir Parameter direction defining the periodic faces
  void closeFaces(int dir);

  //! \brief Constrains all DOFs on a given boundary face.
  //! \param[in] dir Parameter direction defining the face to constrain
  //! \param[in] dof Which DOFs to constrain at each node on the face
  //! \param[in] value Prescribed value of the constrained DOFs
  void constrainFace(int dir, int dof = 123, double value = 0.0);

  //! \brief Constrains all DOFs along a line on a given boundary face.
  //! \param[in] fdir Parameter direction defining the face to constrain
  //! \param[in] ldir Parameter direction defining the line to constrain
  //! \param[in] xi Parameter value defining the line to constrain
  //! \param[in] dof Which DOFs to constrain at each node along the line
  //! \param[in] value Prescribed value of the constrained DOFs
  //!
  //! \details The parameter \a xi has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. The line to
  //! constrain goes along the parameter direction \a ldir in the face with
  //! with normal in parameter direction \a fdir, and positioned along the third
  //! parameter direction as indicated by \a xi. The actual value of \a xi
  //! is converted to the integer value closest to \a xi*n, where \a n is the
  //! number of nodes (control points) in that parameter direction.
  void constrainLine(int fdir, int ldir, double xi,
		     int dof = 123, double value = 0.0);

  //! \brief Constrains a corner node identified by the three parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] K Parameter index in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] value Prescribed value of the constrained DOFs
  //!
  //! \details The sign of the three indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  void constrainCorner(int I, int J, int K, int dof = 123, double value = 0.0);

  //! \brief Constrains a node identified by three relative parameter values.
  //! \param[in] xi Parameter in u-direction
  //! \param[in] eta Parameter in v-direction
  //! \param[in] zeta Parameter in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] value Prescribed value of the constrained DOFs
  //!
  //! \details The parameter values have to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  void constrainNode(double xi, double eta, double zeta,
		     int dof = 123, double value = 0.0);

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag (see below)
  //!
  //! \details The face orientation flag \a norient must be in range [0,7].
  //! When interpreted as a binary number, its 3 digits are decoded as follows:
  //! - left digit = 1: The u and v parameters of the neighbor face are swapped
  //! - middle digit = 1: Parameter \a u in neighbor patch face is reversed
  //! - right digit = 1: Parameter \a v in neighbor patch face is reversed
  bool connectPatch(int face, VolumePatch& neighbor,
		    int nface, int norient = 0);

  //! \brief Assembles coefficient matrices and right-hand-side vector.
  //! \param sys The linear system of equations
  //! \param[in] sam Data for finite element assembly management
  //! \param[in] gravity Gravitation vector
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  //! \param[in] displ Solution vector in DOF-order
  bool assembleSystem(LinEqSystem& sys, const SAM& sam, const Vec3& gravity,
		      int nGauss = 4, const Vector& displ = Vector());

  //! \brief Assembles right-hand-side vector due to surface traction on a face.
  //! \param S The right-hand-side vector
  //! \param[in] sam Data for finite element assembly management
  //! \param[in] t The surface traction function
  //! \param[in] dir Local index of the face subjected to the traction
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  //! \param[out] trac Evaluated tractions at integration points (optional)
  bool assembleForces(SystemVector& S, const SAM& sam, const TractionFunc& t,
		      int dir, int nGauss = 4, std::map<Vec3,Vec3>* trac = 0);

  //! \brief Evaluates some norms of the finite element solution.
  //! \details The energy norm of the solution is computed by numerical
  //! integration of \f$a({\bf u}^h,{\bf u}^h)\f$ over the patch.
  //! If an analytical solution is available, the norm of the exact error
  //! \f$a({\bf u}-{\bf u}^h,{\bf u}-{\bf u}^h)\f$ is computed as well.
  //! \param[out] gNorm Global norms
  //! \param[out] eNorm Element norms
  //! \param[in] nGauss Numerical integration scheme (number of points in 1D)
  //! \param[in] displ Solution vector in DOF-order
  //! \param[in] sol Pointer to analytical stress field (optional)
  bool solutionNorms(Vector& gNorm, Matrix& eNorm, int nGauss,
		     const Vector& displ, const TensorFunc* sol = 0);

  //! \brief Creates a hexahedron element model of this patch for visualization.
  //! \param[out] grid The generated hexahedron grid
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \note The number of element nodes must be set in \a grid on input.
  bool convertToElementBlock(ElementBlock& grid, const int* npe) const;
  //! \brief Extracts element results for this patch from a global vector.
  //! \param[in] globRes Global matrix of element results
  //! \param[out] elmRes Element results for this patch
  void extractElmRes(const Matrix& globRes, Matrix& elmRes) const;
  //! \brief Extracts nodal displacements for this patch from the global vector.
  //! \param[in] solution Global solution vector in DOF-order
  //! \param[out] displ Nodal displacement vector for this patch
  void extractSolution(const Vector& solution, Vector& displ) const;
  //! \brief Evaluates the displacement field at all visualization points.
  //! \param[out] dField Displacement field
  //! \param[in] displ Solution vector in DOF-order
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] cs Local coordinate system
  bool evalDisplField(Matrix& dField, const Vector& displ,
		      const int* npe, const LocalSystem* cs = 0) const;
  //! \brief Evaluates the stress field at all visualization points.
  //! \param[out] sField Stress field
  //! \param[in] displ Solution vector in DOF-order
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] cs Local coordinate system
  bool evalStressField(Matrix& sField, const Vector& displ,
		       const int* npe, const LocalSystem* cs = 0) const;

  //! \brief Calculates von Mises stress field from a stress tensor field.
  //! \param[in]  sigma Stress tensor field
  //! \param[out] vm von Mises stress field
  static void vonMises(const Matrix& sigma, Vector& vm);

  //! \brief If \e true, matching nodes in two adjacent patches are merged.
  static bool mergeDuplNodes; //!< If \e false, MPC-equations are used.
  static bool swapJac; //!< Should the sign of the Jacobian be swapped?
  static int  splineEvalMethod; //!< Which spline evaluation method to use

protected:
  //! \brief Adds a general multi-point-constraint equation to the patch.
  //! \param mpc Pointer to an MPC-object
  bool addMPC(MPC* mpc);
  //! \brief Creates and adds a single-point constraint to the patch.
  //! \param[in] node Global node number of the node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  //! \param[in] value The prescribed value
  bool addSPC(int node, int dir, double value = 0.0);
  //! \brief Creates and adds a periodicity constraint to the patch.
  //! \param[in] master Global node number of the master node
  //! \param[in] slave Global node number of the slave node to constrain
  //! \param[in] dir Which local DOF to constrain (1, 2, 3)
  bool addPeriodicity(int master, int slave, int dir);
  //! \brief Creates periodicity constraints between two nodes in the patch.
  //! \param[in] master Global node number of the master node
  //! \param[in] slave Global node number of the slave node to constrain
  //! \param[in] code Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  void makePeriodic(int master, int slave, int code = 123);
  //! \brief Constrains DOFs in the given node to the given value.
  //! \param[in] node 1-based node index local to current patch
  //! \param[in] code Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  //! \param[in] value The prescribed value
  void prescribe(int node, int code = 123, double value = 0.0);
  //! \brief Constrains DOFs in the given node to zero.
  //! \param[in] node 1-based node index local to current patch
  //! \param[in] code Which local DOFs to constrain (1, 2, 3, 12, 23, 123)
  void fix(int node, int code = 123);

  //! \brief Recursive method used by \a resolveMPCchains.
  static bool resolveMPCchain(const MPCSet& allMPCs, MPC* mpc);

  //! \brief Merges the given two nodes into one.
  //! \details The global node number of the node with the highest number
  //! is changed into the number of the other node.
  static void mergeNodes(IJK& node1, IJK& node2);

  //! \brief Calculates the parameter values for all visualization points in
  //! one direction.
  //! \param[out] par Parameter values for all visualization points.
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] nSegSpan Number of visualization segments over each non-zero
  //! knot-spans
  bool getGridParameters(RealArray& par, int dir, int nSegSpan) const;

  //! \brief Returns the volume in the parameter space for an element.
  //! \param[in] iel Element index
  double getParametricVolume(int iel) const;
  //! \brief Returns boundary face area in the parameter space for an element.
  //! \param[in] iel Element index
  //! \param[in] dir Local face index of the boundary face
  double getParametricArea(int iel, int dir) const;
  //! \brief Sets up the constitutive matrix for this patch.
  //! \param[out] C 6\f$\times\f$6-matrix, representing the material tensor
  //! \param[in] inverse If \e true, set up the inverse matrix instead
  bool getMaterialMatrix(Matrix& C, bool inverse = false) const;
  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel Element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  bool getElementCoordinates(Matrix& X, int iel) const;
  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  void getNodalCoordinates(Matrix& X) const;
  //! \brief Returns an index into the internal coefficient array for a node.
  //! \param[in] inod 1-based node index local to current patch
  int  coeffInd(int inod) const;

  //! \brief Evaluates the displacement field at all visualization points.
  //! \details Memory intensive version, used when \a splineEvalMethod = 1.
  //! \param[out] dField Displacement field
  //! \param[in] displ Solution vector in DOF-order
  //! \param[in] gpar Parameter values of the vizualization points
  //! \param[in] cs Local coordinate system
  bool evalDisplField1(Matrix& dField, const Vector& displ,
		       const RealArray* gpar, const LocalSystem* cs = 0) const;
  //! \brief Evaluates the displacement field at all visualization points.
  //! \details More efficient version, used when \a splineEvalMethod = 2.
  //! \param[out] dField Displacement field
  //! \param[in] displ Solution vector in DOF-order
  //! \param[in] gpar Parameter values of the vizualization points
  //! \param[in] cs Local coordinate system
  bool evalDisplField2(Matrix& dField, const Vector& displ,
		       const RealArray* gpar, const LocalSystem* cs = 0) const;

  //! \brief Evaluates the stress field at all visualization points.
  //! \details Memory intensive version, used when \a splineEvalMethod = 1.
  //! \param[out] sField Stress field
  //! \param[in] displ Solution vector in DOF-order
  //! \param[in] gpar Parameter values of the vizualization points
  //! \param[in] cs Local coordinate system
  bool evalStressField1(Matrix& sField, const Vector& displ,
		        const RealArray* gpar, const LocalSystem* cs = 0) const;
  //! \brief Evaluates the stress field at all visualization points.
  //! \details More efficient version, used when \a splineEvalMethod = 2.
  //! \param[out] sField Stress field
  //! \param[in] displ Solution vector in DOF-order
  //! \param[in] gpar Parameter values of the vizualization points
  //! \param[in] cs Local coordinate system
  bool evalStressField2(Matrix& sField, const Vector& displ,
		        const RealArray* gpar, const LocalSystem* cs = 0) const;

  //! \brief Calculates the Strain-displacement matrix \b B.
  static void formBmatrix(Matrix& B, const Matrix& dNdX);

private:
  // Geometry and physical properties
  Go::SplineVolume* svol;    //!< Pointer to the actual spline volume object
  bool              swapW;   //!< Has the w-parameter direction been swapped?
  double            E;       //!< Young's modulus
  double            nu;      //!< Poisson's ratio
  double            rho;     //!< Mass density
  // Finite element data structures
  IntMat            MNPC;    //!< Matrix of Nodal Point Correspondance
  std::vector<int>  MLGE;    //!< Matrix of Local to Global Element numbers
  std::vector<IJK>  nodeInd; //!< IJK-triplets for the control points (nodes)
  std::vector<BC>   BCode;   //!< Boundary condition codes
  MPCSet            mpcs;    //!< All multi-point constraints
  static int        gEl;     //!< Global element counter
  static int        gNod;    //!< Global node counter
};

#endif
