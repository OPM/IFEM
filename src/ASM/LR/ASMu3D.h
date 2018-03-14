// $Id$
//==============================================================================
//!
//! \file ASMu3D.h
//!
//! \date January 2013
//!
//! \author Kjetil A. Johannessen / NTNU
//!
//! \brief Driver for assembly of unstructured 3D spline FE models.
//!
//==============================================================================

#ifndef _ASM_U3D_H
#define _ASM_U3D_H

#include "ASMunstruct.h"
#include "ASM3D.h"
#include "LRSpline/LRSpline.h"
#include "ThreadGroups.h"
#include <memory>

class FiniteElement;

namespace Go {
  class SplineVolume;
}

namespace LR {
  class LRSplineVolume;
}


/*!
  \brief Driver for assembly of unstructured 3D spline FE models.
  \details This class contains methods common for 3D LR-spline patches.
*/

class ASMu3D : public ASMunstruct, public ASM3D
{
public:
  //! \brief Default constructor.
  explicit ASMu3D(unsigned char n_f = 3);
  //! \brief Copy constructor.
  ASMu3D(const ASMu3D& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMu3D() {}

  //! \brief Returns the spline volume representing the geometry of this patch.
  LR::LRSplineVolume* getVolume() const { return lrspline.get(); }

  //! \brief Returns the spline volume representing the basis of this patch.
  virtual const LR::LRSplineVolume* getBasis(int = 1) const { return lrspline.get(); }
  //! \brief Returns the spline volume representing the basis of this patch.
  virtual LR::LRSplineVolume* getBasis(int = 1) { return lrspline.get(); }


  // Methods for model generation and refinement
  // ===========================================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineVolume object to given stream.
  virtual bool write(std::ostream&, int = 0) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! and the arrays of global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel 1-based element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const;

  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const;

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ);

  //! \brief Returns the node indices for a given face.
  IntVec getFaceNodes(int face, int basis = 1, int orient = -1) const;

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (0 for all)
  //! \param[in] orient Orientation of boundary (used for sorting)
  //! \param[in] local If \e true, return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes, int basis,
                                int, int orient, bool local = false) const;

  //! \brief Returns the node index for a given corner.
  //! \param[in] I -1 or +1 for either umin or umax corner
  //! \param[in] J -1 or +1 for either vmin or vmax corner
  //! \param[in] K -1 or +1 for either wmin or wmax corner
  //! \param[in] basis which basis to consider (for mixed methods)
  virtual int getCorner(int I, int J, int K, int basis) const;

  //! \brief Returns the (1-indexed) node indices for a given edge.
  //! \param[in] lEdge index to local edge (1,2,...12)
  //! \param[in] open include end points or not
  //! \param[in] basis which basis to consider (for mixed methods)
  //! \param[in] orient local orientation of indices, see ASMunstruct::Sort
  virtual IntVec getEdge(int lEdge, bool open, int basis, int orient) const;

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction
  virtual bool getOrder(int& p1, int& p2, int& p3) const;

  //! \brief Refines the parametrization by inserting tensor knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int dir, int nInsert);
  using ASMunstruct::refine;
  //! \brief Refines the parametrization by inserting extra tensor knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  virtual bool refine(int dir, const RealArray& xi);
  //! \brief Raises the order of the tensor spline object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  //! \param[in] rw Number of times to raise the order in w-direction
  virtual bool raiseOrder(int ru, int rv, int rw);

  //! \brief Defines the minimum element volume for adaptive refinement.
  //! \param[in] nrefinements Maximum number of adaptive refinement levels
  virtual double getMinimumSize(int nrefinements) const;
  //! \brief Checks if the specified element is larger than the minimum size.
  //! \param[in] elmId Global/patch-local element index
  //! \param[in] globalNum If \e true, \a elmId is global otherwise patch-local
  virtual bool checkElementSize(int elmId, bool globalNum = true) const;


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains all DOFs on a given boundary face.
  //! \param[in] dir Parameter direction defining the face to constrain
  //! \param[in] open If \e true, exclude all points along the face boundary
  //! \param[in] dof Which DOFs to constrain at each node on the face
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain face for
  virtual void constrainFace(int dir, bool open, int dof,
                             int code = 0, char basis = 1);
  //! \brief Constrains all DOFs in local directions on a given boundary face.
  //! \param[in] dir Parameter direction defining the face to constrain
  //! \param[in] open If \e true, exclude all points along the face boundary
  //! \param[in] dof Which local DOFs to constrain at each node on the face
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] project If \e true, the local axis directions are projected
  //! \param[in] T1 Desired global direction of first local tangent direction
  //! \return Number of additional nodes added due to local axis constraints
  virtual size_t constrainFaceLocal(int dir, bool open, int dof, int code = 0,
                                    bool project = false, char T1 = '\0');

  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] lEdge Local index [1,12] of the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which DOFs to constrain at each node on the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  virtual void constrainEdge(int lEdge, bool open, int dof,
                             int code = 0, char = 1);

  //! \brief Constrains all DOFs along a line on a given boundary face.
  //! \param[in] fdir Parameter direction defining the face to constrain
  //! \param[in] ldir Parameter direction defining the line to constrain
  //! \param[in] xi Parameter value defining the line to constrain
  //! \param[in] dof Which DOFs to constrain at each node along the line
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The parameter \a xi has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. The line to
  //! constrain goes along the parameter direction \a ldir in the face with
  //! with normal in parameter direction \a fdir, and positioned along the third
  //! parameter direction as indicated by \a xi. The actual value of \a xi
  //! is converted to the integer value closest to \a xi*n, where \a n is the
  //! number of nodes (control points) in that parameter direction.
  virtual void constrainLine(int fdir, int ldir, double xi, int dof,
                             int code = 0, char = 1);

  //! \brief Constrains a corner node identified by the three parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] K Parameter index in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The sign of the three indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  virtual void constrainCorner(int I, int J, int K, int dof,
                               int code = 0, char = 1);
  //! \brief Constrains a node identified by three relative parameter values.
  //! \param[in] xi Parameter in u-direction
  //! \param[in] eta Parameter in v-direction
  //! \param[in] zeta Parameter in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The parameter values have to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual void constrainNode(double xi, double eta, double zeta, int dof,
                             int code = 0, char = 1);

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag (see below)
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //!
  //! \details The face orientation flag \a norient must be in range [0,7].
  //! When interpreted as a binary number, its 3 digits are decoded as follows:
  //! - left digit = 1: The u and v parameters of the neighbor face are swapped
  //! - middle digit = 1: Parameter \a u in neighbor patch face is reversed
  //! - right digit = 1: Parameter \a v in neighbor patch face is reversed
  virtual bool connectPatch(int face, ASM3D& neighbor, int nface, int norient,
                            int = 0, bool coordCheck = true, int = 1);

  /* More multipatch stuff, maybe later...
  //! \brief Makes two opposite boundary faces periodic.
  //! \param[in] dir Parameter direction defining the periodic faces
  //! \param[in] basis Which basis to connect (mixed methods), 0 means both
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeFaces(int dir, int basis = 0, int master = 1);
  */


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates a boundary integral over a patch face.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index [1,6] of the boundary face
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lEdge Local index [1,12] of the patch edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrateEdge(Integrand& integrand, int lEdge,
                             GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  //! \param[in] time Current time
  //! \param[in] g2l Pointer to global-to-local node number mapping
  virtual bool updateDirichlet(const std::map<int,RealFunc*>& func,
                               const std::map<int,VecFunc*>& vfunc, double time,
                               const std::map<int,int>* g2l = nullptr);

  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameters in range [0.0,1.0] of the point
  //! \param[out] param The (u,v,w) parameters of the point in knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that matches the point, if any
  //! \return 0 if no node (control point) matches this point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const;

  //! \brief Calculates parameter values for visualization nodal points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  virtual bool getGridParameters(RealArray& prm, int dir, int nSegSpan) const;

  //! \brief Creates a hexahedron element model of this patch for visualization.
  //! \param[out] grid The generated hexahedron grid
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int* npe) const;

  //! \brief Evaluates the primary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector in DOF-order
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] nf If nonzero, mixed evaluates nf fields on first basis
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int* npe, int nf) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] deriv Derivative order to return
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool = false,
                            int deriv = 0, int = 0) const;

  //! \brief Evaluates and interpolates a function over a given geometry.
  //! \param[in] func The function to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] time Current time
  virtual bool evaluate(const FunctionBase* func, RealArray& vec,
                        int, double time) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] project Flag indicating result recovery method
  //! (0=none, 'D'=direct evaluation, 'S'=superconvergent recovery)
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! If \a npe is null, the solution is recovered or evaluated at the Greville
  //! points and then projected onto the spline basis to obtain the control
  //! point values, which then are returned through \a sField.
  //! If \a npe is not null and \a project is defined, the solution is also
  //! projected onto the spline basis, and then evaluated at the \a npe points.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const int* npe, char project = '\0') const;

  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] old_basis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVars(const LR::LRSpline* old_basis,
                                   const RealArray& oldVar, RealArray& newVar,
                                   int nGauss) const;
  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] old_basis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVarsN(const LR::LRSpline* old_basis,
                                    const RealArray& oldVar, RealArray& newVar,
                                    int nGauss) const;
  using ASMunstruct::transferCntrlPtVars;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] old_basis The LR-spline basis to transfer from
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferCntrlPtVars(const LR::LRSpline* old_basis,
                                   RealArray& newVar, int nGauss) const;

  //! \brief Refines the parametrization based on a mesh density function.
  //! \param[in] refC Mesh refinement criteria function
  //! \param[in] refTol Mesh refinement threshold
  virtual bool refine(const RealFunc& refC, double refTol);

private:
  //! \brief Struct representing an inhomogeneous Dirichlet boundary condition.
  struct DirichletFace
  {
    LR::LRSplineVolume  *lr;       //!< Pointer to the right object (in case of multiple bases)
    LR::parameterEdge   edg;       //!< Which face is this
    IntVec              MLGE;      //!< Local-to-Global Element numbers
    IntVec              MLGN;      //!< Local-to-Global Nodal numbers
    IntMat              MNPC;      //!< Matrix of Nodal-Point Correpondanse
    int                 dof;       //!< Local DOF to constrain along the boundary
    int                 code;      //!< Inhomogeneous Dirichlet condition code
    int                 basis;     //!< Index to the basis used
    int                 corners[4];//!< Index of the four corners of this face

    //! \brief Default constructor.
    DirichletFace(int numbBasis, int numbElements, int d = 0, int c = 0, int b = 1)
    : MLGE(numbElements), MNPC(numbElements), dof(d), code(c), basis(b) {}
  };

  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  LR::LRSplineVolume* projectSolution(const IntegrandBase& integrand) const;
  //! \brief Projects the secondary solution using a superconvergent approach.
  //! \param[in] integrand Object with problem-specific data and methods
  LR::LRSplineVolume* scRecovery(const IntegrandBase& integrand) const;

  //! \brief Interpolates an LR spline at a given set of parametric coordinates.
  //! \param[in] upar Parametric interpolation points in the u-direction
  //! \param[in] vpar Parametric interpolation points in the v-direction
  //! \param[in] wpar Parametric interpolation points in the v-direction
  //! \param[in] points Interpolation values stored as one point per matrix row
  //! \return A LRSplineVolume representation of the interpolated points
  LR::LRSplineVolume* regularInterpolation(const RealArray& upar,
                                           const RealArray& vpar,
                                           const RealArray& wpar,
                                           const Matrix& points) const;

public:
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual LR::LRSpline* evalSolution(const IntegrandBase& integrand) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! We assume that the parameter value array \a gpar contains
  //! the \a u, \a v and \a w parameters directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const RealArray* gpar, bool = false) const;

  //! \brief Projects inhomogenuous dirichlet conditions by continuous L2-fit.
  //! \param[in] face low-level face information needed to do integration
  //! \param[in] values inhomogenuous function which is to be fitted
  //! \param[out] result fitted value in terms of control-point values
  //! \param[in] time time used in dynamic problems
  bool faceL2projection (const DirichletFace& face,
                         const FunctionBase& values,
                         Real2DMat& result,
                         double time) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Assembles L2-projection matrices for the secondary solution.
  //! \param[out] A Left-hand-side matrix
  //! \param[out] B Right-hand-side vectors
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e false, a discrete L2-projection is used
  virtual bool assembleL2matrices(SparseMatrix& A, StdVector& B,
                                  const IntegrandBase& integrand,
                                  bool continuous) const;

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag (see \a connectPatch)
  //! \param[in] basis Which basis to connect the nodes for (mixed methods)
  //! \param[in] slave 0-based index of the first slave node in this basis
  //! \param[in] master 0-based index of the first master node in this basis
  //! \param[in] coordCheck False to disable coordinate checks
  bool connectBasis(int face, ASMu3D& neighbor, int nface, int norient,
                    int basis = 1, int slave = 0, int master = 0,
                    bool coordCheck = true, int = 1);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] iel 1-based element index
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  void getGaussPointParameters(RealArray& uGP, int dir, int nGauss,
                               int iel, const double* xi) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  bool getGrevilleParameters(RealArray& prm, int dir) const;

  //! \brief Calculates parameter values for the Quasi-Interpolation points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  bool getQuasiInterplParameters(RealArray& prm, int dir) const;

  //! \brief Returns the volume in the parameter space for an element.
  //! \param[in] iel 1-based element index
  double getParametricVolume(int iel) const;

  //! \brief Returns boundary face area in the parameter space for an element.
  //! \param[in] iel 1-based element index
  //! \param[in] dir Local face index of the boundary face
  double getParametricArea(int iel, int dir) const;

  //! \brief Computes the element corner coordinates.
  //! \param[in] iel 1-based element index
  //! \param[out] XC Coordinates of the element corners
  //! \return Characteristic element size
  double getElementCorners(int iel, std::vector<Vec3>& XC) const;

  //! \brief Evaluate all basis functions and first derivatives on one element
  void evaluateBasis(int iel, double u, double v, double w,
                     Vector& N, Matrix& dNdu, int basis) const;

  //! \brief Evaluate all basis functions and \a derivs number of derivatives on one element
  void evaluateBasis(FiniteElement &el, int derivs, int basis = 1) const;

  //! \brief Evaluate all basis functions and first derivatives on one element
  void evaluateBasis(FiniteElement &el, Matrix &dNdu,
                     const Matrix &C, const Matrix &B, int basis = 1) const;

  //! \brief Evaluate all basis functions and first derivatives on one element
  void evaluateBasis(FiniteElement &el, Matrix &dNdu, int basis = 1) const;

  //! \brief Evaluate all basis functions and second order derivatives on one element
  void evaluateBasis(FiniteElement &el, Matrix &dNdu, Matrix3D& d2Ndu2, int basis = 1) const;

  using ASMunstruct::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM If \e true ignore global multipliers in sanity check
  void generateThreadGroups(const Integrand& integrand, bool silence,
                            bool ignoreGlobalLM);

  //! \brief Remap element wise errors to basis functions.
  //! \param     errors The remapped errors
  //! \param[in] origErr The element wise errors on the geometry mesh
  //! \param[in] elemErrors If true, map to elements instead of basis functions
  virtual void remapErrors(RealArray& errors,
                           const RealArray& origErr, bool elemErrors) const;

  //! \brief Extends the refinement domain with information for neighbors.
  //! \param refineIndices List of basis functions to refine
  //! \param neighborIndices Basis functions to refine from neighbor patches
  virtual void extendRefinementDomain(IntSet& refineIndices,
                                      const IntSet& neighborIndices) const;

public:
  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

protected:
  std::shared_ptr<LR::LRSplineVolume> lrspline; //!< Pointer to the LR-spline volume object

  Go::SplineVolume* tensorspline; //!< Pointer to original tensor spline object
  // The tensor spline object is kept for backward compatability with the REFINE
  // and RAISEORDER key-words, although we take note that there is a possibility
  // of optimization since all mapping values and Jacobians may be performed on
  // this object for increased efficiency.

  //! Inhomogeneous Dirichlet boundary condition data
  std::vector<DirichletFace> dirich;
  int myGeoBasis; //!< Used with mixed

  const Matrices& bezierExtract; //!< Bezier extraction matrices
  Matrices      myBezierExtract; //!< Bezier extraction matrices

  ThreadGroups threadGroups; //!< Element groups for multi-threaded assembly
  mutable double vMin; //!< Minimum element volume for adaptive refinement
};

#endif
