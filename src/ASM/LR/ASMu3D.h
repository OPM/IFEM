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

#include "ASMLRSpline.h"
#include "ASM3D.h"
#include "BasisFunctionCache.h"
#include "LRSpline/LRSpline.h"
#include "ThreadGroups.h"
#include <memory>

class FiniteElement;

namespace utl {
  class Point;
}

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

class ASMu3D : public ASMLRSpline, public ASM3D
{
protected:
  //! \brief Implementation of basis function cache.
  class BasisFunctionCache : public ::BasisFunctionCache<3>
  {
  public:
    //! \brief The constructor initializes the class.
    //! \param pch Patch the cache is for
    //! \param useBezier True to use bezier extraction
    BasisFunctionCache(const ASMu3D& pch, bool useBezier = true);

    //! \brief Constructor reusing quadrature info from another instance.
    //! \param cache Instance holding quadrature information
    //! \param b Basis to use
    BasisFunctionCache(const BasisFunctionCache& cache, int b);

    //! \brief Empty destructor.
    virtual ~BasisFunctionCache() = default;

  protected:
    //! \brief Implementation specific initialization.
    bool internalInit() override;

    //! \brief Calculates basis function info in a single integration point.
    //! \param el Element of integration point (0-indexed)
    //! \param gp Integration point on element (0-indexed)
    //! \param reduced If true, returns values for reduced integration scheme
    BasisFunctionVals calculatePt(size_t el, size_t gp, bool reduced) const override;

    //! \brief Calculates basis function info in a single integration point.
    //! \param fe Finite element information in integration point
    //! \param el Element of integration point (0-indexed)
    //! \param du Element size in parameter space
    //! \param gp Integration point on element (0-indexed)
    //! \param reduced If true, returns values for reduced integration scheme
    BasisFunctionVals calculatePrm(FiniteElement& fe,
                                   const std::array<double,3>& du,
                                   size_t el, size_t gp, bool reduced) const;

    //! \brief Calculates basis function info in all integration points.
    void calculateAll() override;

    //! \brief Struct holding bezier extraction matrices.
    struct BezierExtract {
      Matrix N;     //!< Basis function values
      Matrix dNdu; //!< Basis function u-derivatives
      Matrix dNdv; //!< Basis function v-derivatives
      Matrix dNdw; //!< Basis function w-derivatives
    };

    bool bezierEnabled; //!< True to enable Bezier extraction
    BezierExtract mainB; //!< Bezier extraction for main basis
    BezierExtract reducedB; //!< Bezier extraction for reduced basis

    const ASMu3D& patch; //!< Reference to patch cache is for

  private:
    //! \brief Configure quadratures.
    bool setupQuadrature();
  };

public:
  //! \brief Default constructor.
  explicit ASMu3D(unsigned char n_f = 3);
  //! \brief Copy constructor.
  ASMu3D(const ASMu3D& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMu3D() {}

  //! \brief Returns the spline volume representing a basis of this patch.
  virtual const LR::LRSplineVolume* getBasis(int basis = 1) const;
  //! \brief Returns the spline volume representing a basis of this patch.
  virtual LR::LRSplineVolume* getBasis(int basis = 1);


  // Methods for model generation and refinement
  // ===========================================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineVolume object to given stream.
  virtual bool write(std::ostream&, int) const;

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
  //! \param[in] forceItg If true return integration basis element coordinates
  virtual bool getElementCoordinates(Matrix& X, int iel, bool forceItg = false) const;

  //! \brief Obtain element neighbours.
  virtual void getElmConnectivities(IntMat& neighs, bool local = false) const;

  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X, bool = false) const;

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ);

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (for mixed methods)
  //! \param[in] orient Orientation of boundary (used for sorting)
  //! \param[in] local If \e true, return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes, int basis,
                                int, int orient, bool local) const;

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lEdge Local index of the boundary edge
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (for mixed methods)
  //! \param[in] orient Local orientation of indices, see ASMLRSpline::Sort()
  //! \param[in] local If \e true, return patch-local numbers
  //! \param[in] open If \e true, exclude edge end points
  virtual void getBoundary1Nodes(int lEdge, IntVec& nodes, int basis,
                                 int orient, bool local,
                                 bool open = false) const;

  //! \brief Returns the node index for a given corner.
  //! \param[in] I -1 or +1 for either umin or umax corner
  //! \param[in] J -1 or +1 for either vmin or vmax corner
  //! \param[in] K -1 or +1 for either wmin or wmax corner
  //! \param[in] basis which basis to consider (for mixed methods)
  virtual int getCorner(int I, int J, int K, int basis) const;

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction
  virtual bool getOrder(int& p1, int& p2, int& p3) const;

  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;

  //! \brief Returns the number of projection nodes for this patch.
  virtual size_t getNoProjectionNodes() const;

  //! \brief Refines the parametrization by inserting tensor knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int dir, int nInsert);
  using ASMLRSpline::refine;
  //! \brief Refines the parametrization by inserting extra tensor knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  virtual bool refine(int dir, const RealArray& xi);
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol);
  //! \brief Raises the order of the tensor spline object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  //! \param[in] rw Number of times to raise the order in w-direction
  //! \param[in] setOrder If \e true, raise order to \a ru, \a rv and \a rw
  virtual bool raiseOrder(int ru, int rv, int rw, bool setOrder);

  //! \brief Creates a separate projection basis for this patch.
  virtual bool createProjectionBasis(bool init);

  //! \brief Sets the minimum element volume for adaptive refinement.
  virtual void setMinimumSize(double size) { vMin = size; }
  //! \brief Defines the minimum element volume for adaptive refinement.
  //! \param[in] nrefinements Maximum number of adaptive refinement levels
  virtual double getMinimumSize(int nrefinements) const;
  //! \brief Checks if the specified element is larger than the minimum size.
  //! \param[in] elmId Global/patch-local element index
  //! \param[in] globalNum If \e true, \a elmId is global otherwise patch-local
  virtual bool checkElementSize(int elmId, bool globalNum = true) const;

  //! \brief Copies the refinement to another spline volume.
  //! \param basis Volume to copy refinement to
  //! \param[in] multiplicity Wanted multiplicity
  void copyRefinement(LR::LRSplineVolume* basis, int multiplicity = 1) const;


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
  //! \param[in] dof Which DOFs to constrain at each node along the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain edge for
  virtual void constrainEdge(int lEdge, bool open, int dof,
                             int code = 0, char basis = 1);

  //! \brief Constrains all DOFs along a line on a given boundary face.
  //! \param[in] fdir Parameter direction defining the face to constrain
  //! \param[in] ldir Parameter direction defining the line to constrain
  //! \param[in] xi Parameter value defining the line to constrain
  //! \param[in] dof Which DOFs to constrain at each node along the line
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain line for
  //!
  //! \details The parameter \a xi has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. The line to
  //! constrain goes along the parameter direction \a ldir in the face with
  //! with normal in parameter direction \a fdir, and positioned along the third
  //! parameter direction as indicated by \a xi. The actual value of \a xi
  //! is converted to the integer value closest to \a xi*n, where \a n is the
  //! number of nodes (control points) in that parameter direction.
  virtual void constrainLine(int fdir, int ldir, double xi, int dof,
                             int code = 0, char basis = 1);

  //! \brief Constrains a corner node identified by the three parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] K Parameter index in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain node for
  //!
  //! \details The sign of the three indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  virtual void constrainCorner(int I, int J, int K, int dof,
                               int code = 0, char basis = 1);
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
                             int code = 0);

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag (see below)
  //! \param[in] coordCheck False to disable coordinate checks
  //! \param[in] thick Thickness of connection
  //!
  //! \details The face orientation flag \a norient must be in range [0,7].
  //! When interpreted as a binary number, its 3 digits are decoded as follows:
  //! - left digit = 1: The u and v parameters of the neighbor face are swapped
  //! - middle digit = 1: Parameter \a u in neighbor patch face is reversed
  //! - right digit = 1: Parameter \a v in neighbor patch face is reversed
  virtual bool connectPatch(int face, ASM3D& neighbor, int nface, int norient,
                            int = 0, bool coordCheck = true, int thick = 1);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Computes the number of boundary integration points in this patch.
  virtual void getNoBouPoints(size_t& nPt, char ldim, char lindx);

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

  //! \brief Integrates a spatial dirac-delta function over a patch.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] param Parameters of the non-zero point of dirac-delta function
  //! \param[in] pval Function value at the specified point
  virtual bool diracPoint(Integrand& integrand, GlobalIntegral& glbInt,
                          const double* param, const Vec3& pval);

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

  //! \brief Returns the element that contains a specified spatial point.
  //! \param[in] param The parameters of the point in the knot-span domain
  //! \return Local element number within the patch that contains the point
  virtual int findElementContaining(const double* param) const;

  //! \brief Searches for the specified Cartesian point in the patch.
  //! \param X The Cartesian coordinates of the point, updated on output
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \return Distance from the point \a X to the found point
  virtual double findPoint(Vec3& X, double* param) const;

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
  //! \param[in] n_f If nonzero, mixed evaluates \a n_f fields on first basis
  //! \param[in] piola If \e true, use piola mapping
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int* npe, int n_f, bool piola) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] deriv Derivative order to return
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool = false,
                            int deriv = 0, int = 0) const;

  //! \brief Evaluates the projected solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] n_f If nonzero, mixed evaluates \a n_f fields on first basis
  virtual bool evalProjSolution(Matrix& sField, const Vector& locSol,
                                const int* npe, int n_f) const;

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

  //! \brief Checks if a separate projection basis is used for this patch.
  virtual bool separateProjectionBasis() const;

  //! \brief Returns a field using the projection basis.
  //! \param[in] coefs The coefficients for the field
  virtual Field* getProjectedField(const Vector& coefs) const;

  //! \brief Returns a field using the projection basis.
  //! \param[in] coefs The coefficients for the field
  virtual Fields* getProjectedFields(const Vector& coefs, size_t = 0) const;

  //! \brief Set master patch for VTF output.
  //! \param pch Master patch to use
  void setOutputMaster(const ASMu3D* pch)
  { outputMaster = pch; }

  //! \brief Extracts element results for this patch from a global vector.
  //! \param[in] globRes Global matrix of element results
  //! \param[out] elmRes Element results for this patch
  //! \param[in] internalFirst Global index of first element in the patch
  virtual void extractElmRes(const Matrix& globRes, Matrix& elmRes,
                             size_t internalFirst) const;

protected:
  //! \brief Struct representing an inhomogeneous Dirichlet boundary condition.
  struct DirichletFace
  {
    LR::LRSplineVolume* lr;        //!< Pointer to the right object (in case of multiple bases)
    LR::parameterEdge   edg;       //!< Which face is this
    IntVec              MLGE;      //!< Local-to-Global Element numbers
    IntVec              MLGN;      //!< Local-to-Global Nodal numbers
    IntMat              MNPC;      //!< Matrix of Nodal-Point Correpondanse
    int                 dof;       //!< Local DOF to constrain along the boundary
    int                 code;      //!< Inhomogeneous Dirichlet condition code
    int                 corners[4];//!< Index of the four corners of this face

    //! \brief The constructor detects the face corners.
    DirichletFace(LR::LRSplineVolume* sv, int dir,
                  int d = 0, int c = 0, int offset = 1);
    //! \brief Returns \e true if basis function \a b is at a corner point.
    bool isCorner(int b) const;
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
  //! \param[in] basis Which basis to interpolate onto
  //! \return A LRSplineVolume representation of the interpolated points
  LR::LRSplineVolume* regularInterpolation(const RealArray& upar,
                                           const RealArray& vpar,
                                           const RealArray& wpar,
                                           const Matrix& points,
                                           int basis) const;

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
  bool faceL2projection(const DirichletFace& face, const FunctionBase& values,
                        Real2DMat& result, double time) const;

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

  using ASMLRSpline::transferCntrlPtVars;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] old_basis The LR-spline basis to transfer from
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferCntrlPtVars(const LR::LRSpline* old_basis,
                                   RealArray& newVar, int nGauss) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Finds the patch-local element numbers on a patch boundary.
  //! \param[out] elms Array of element numbers
  //! \param[in] lIndex Local index of the boundary face
  //! \param[in] orient Orientation of boundary (used for sorting)
  virtual void findBoundaryElms(IntVec& elms, int lIndex, int orient) const;

  //! \brief Assembles L2-projection matrices for the secondary solution.
  //! \param[out] A Left-hand-side matrix
  //! \param[out] B Right-hand-side vectors
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e false, a discrete L2-projection is used
  virtual bool assembleL2matrices(SparseMatrix& A, StdVector& B,
                                  const L2Integrand& integrand,
                                  bool continuous) const;

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag, see connectPatch()
  //! \param[in] basis Which basis to connect the nodes for (mixed methods)
  //! \param[in] slave 0-based index of the first slave node in this basis
  //! \param[in] master 0-based index of the first master node in this basis
  //! \param[in] coordCheck False to disable coordinate checks
  //! \param[in] thick Thickness of connection
  bool connectBasis(int face, ASMu3D& neighbor, int nface, int norient,
                    int basis = 1, int slave = 0, int master = 0,
                    bool coordCheck = true, int thick = 1);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] iel 1-based element index
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  //! \param[in] spline If given get gauss points for this spline instead of integration basis
  void getGaussPointParameters(RealArray& uGP, int dir, int nGauss,
                               int iel, const double* xi,
                               const LR::LRSplineVolume* spline = nullptr) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] basisNum Which basis to get Greville point parameters for
  bool getGrevilleParameters(RealArray& prm, int dir, int basisNum = 1) const;

  //! \brief Calculates parameter values for the Quasi-Interpolation points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  bool getQuasiInterplParameters(RealArray& prm, int dir) const;

  //! \brief Returns boundary face area in the parameter space for an element.
  //! \param[in] iel 1-based element index
  //! \param[in] dir Local face index of the boundary face
  double getParametricArea(int iel, int dir) const;

  //! \brief Computes the element corner coordinates.
  //! \param[in] iel 1-based element index
  //! \param[out] XC Coordinates of the element corners
  //! \param[out] uC Spline parameters of the element corners (optional)
  //! \return Characteristic element size
  double getElementCorners(int iel, std::vector<Vec3>& XC,
                           RealArray* uC = nullptr) const;
  //! \brief Computes the element corner coordinates and parameters.
  //! \param[in] iel 1-based element index
  //! \param[out] XC Coordinates and parameters of the element corners
  void getCornerPoints(int iel, std::vector<utl::Point>& XC) const;

  //! \brief Evaluate all basis functions and first derivatives on one element
  void evaluateBasis(int iel, double u, double v, double w,
                     Vector& N, Matrix& dNdu, int basis = 1) const;

  //! \brief Evaluate all basis functions and first derivatives on one element
  void evaluateBasis(Vector& N, Matrix& dNdu,
                     const Matrix& C, const Matrix& B) const;

  //! \brief Evaluate all basis functions and first derivatives on one element
  void evaluateBasis(int iel, FiniteElement& fe, Matrix& dNdu,
                     int basis = 1) const;

  //! \brief Evaluate all basis functions and second order derivatives on one element
  void evaluateBasis(int iel, double u, double v, double w,
                     Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2, int basis = 1) const;

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] iel 0-based local element index
  //! \param[in] param The (u,v,w) parameters of the point in knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that matches the point, if any
  //! \return 0 if no node (control point) matches this point
  virtual int evalPoint(int iel, const double* param, Vec3& X) const;

  using ASMLRSpline::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM If \e true ignore global multipliers in sanity check
  void generateThreadGroups(const Integrand& integrand, bool silence,
                            bool ignoreGlobalLM);

  //! \brief Generate element groups from a partition.
  virtual void generateThreadGroupsFromElms(const std::vector<int>& elms);

  //! \brief Hook for changing number of threads.
  virtual void changeNumThreads();

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

  //! \brief Converts current tensor spline object to LR-spline.
  std::shared_ptr<LR::LRSplineVolume> createLRfromTensor();

public:
  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

protected:
  std::shared_ptr<LR::LRSplineVolume> lrspline; //!< The LR-spline volume object

  const ASMu3D* outputMaster = nullptr; //!< Master patch to use for VTF output

  Go::SplineVolume* tensorspline; //!< Pointer to original tensor spline object
  Go::SplineVolume* tensorPrjBas; //!< Pointer to tensor spline projection base
  // The tensor spline object is kept for backward compatability with the REFINE
  // and RAISEORDER key-words, although we take note that there is a possibility
  // of optimization since all mapping values and Jacobians may be performed on
  // this object for increased efficiency.

  //! Inhomogeneous Dirichlet boundary condition data
  std::vector<DirichletFace> dirich;

  const Matrices& bezierExtract; //!< Bezier extraction matrices
  Matrices      myBezierExtract; //!< Bezier extraction matrices

  //! Basis function cache
  std::vector<std::unique_ptr<BasisFunctionCache>> myCache;

private:
  mutable double vMin; //!< Minimum element volume for adaptive refinement
};

#endif
