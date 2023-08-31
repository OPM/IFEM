// $Id$
//==============================================================================
//!
//! \file ASMu2D.h
//!
//! \date September 2011
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D spline FE models.
//!
//==============================================================================

#ifndef _ASM_U2D_H
#define _ASM_U2D_H

#include "ASMLRSpline.h"
#include "ASM2D.h"
#include "BasisFunctionCache.h"
#include "Interface.h"
#include "LRSpline/LRSpline.h"
#include "ThreadGroups.h"
#include <memory>

class FiniteElement;

namespace utl {
  class Point;
}

namespace Go {
  class SplineCurve;
  class SplineSurface;
  struct BasisPtsSf;
  struct BasisDerivsSf;
  struct BasisDerivsSf2;
  struct BasisDerivsSf3;
}

namespace LR {
  class LRSplineSurface;
}


/*!
  \brief Driver for assembly of unstructured 2D spline FE models.
  \details This class contains methods common for 2D LR-spline patches.
*/

class ASMu2D : public ASMLRSpline, public ASM2D
{
protected:
  //! \brief Implementation of basis function cache.
  class BasisFunctionCache : public ::BasisFunctionCache<2>
  {
  public:
    //! \brief The constructor initializes the class.
    //! \param pch Patch the cache is for
    //! \param plcy Cache policy to use
    //! \param b Basis to use
    BasisFunctionCache(const ASMu2D& pch, ASM::CachePolicy plcy, int b);

    //! \brief Constructor reusing quadrature info from another instance.
    //! \param cache Instance holding quadrature information
    //! \param b Basis to use
    BasisFunctionCache(const BasisFunctionCache& cache, int b);

    //! \brief Empty destructor.
    virtual ~BasisFunctionCache() = default;

  protected:
    //! \brief Implementation specific initialization.
    bool internalInit() override;

    //! \brief Implementation specific cleanup.
    void internalCleanup() override;

    //! \brief Calculates basis function info in a single integration point.
    //! \param el Element of integration point (0-indexed)
    //! \param gp Integration point on element (0-indexed)
    //! \param reduced If true, returns values for reduced integration scheme
    BasisFunctionVals calculatePt(size_t el, size_t gp, bool reduced) const override;

    //! \brief Calculates basis function info in all integration points.
    void calculateAll() override;

  protected:
    const ASMu2D& patch; //!< Reference to patch cache is for
    int basis; //!< Basis to use

private:
    //! \brief Configure quadratures.
    bool setupQuadrature();
  };

public:
  //! \brief Base class that checks if an element has interface contributions.
  class InterfaceChecker : public ASM::InterfaceChecker
  {
  public:
    //! \brief The constructor initialises the reference to current patch.
    explicit InterfaceChecker(const ASMu2D& pch);
    //! \brief Empty destructor.
    virtual ~InterfaceChecker() {}
    //! \brief Returns non-zero if the specified element have contributions.
    //! \param[in] iel Element number
    virtual short int hasContribution(int iel, int = -1,
                                      int = -1, int = -1) const;
    //! \brief Get intersections for a given element edge.
    //! \param[in] iel Element index (1-based)
    //! \param[in] edge Edge to get intersections for (1..4)
    //! \param[out] cont If not null, the intersection continuity is given here
    const RealArray& getIntersections(int iel, int edge,
                                      int* cont = nullptr) const;

  protected:
    const ASMu2D& myPatch; //!< Reference to the patch being integrated

    //! \brief Struct describing an intersection of mesh lines.
    struct Intersection {
      int continuity = 0; //!< Continuity across intersection
      RealArray pts;      //!< Intersection points
    };

    //! Intersections for elements. Key: element << 4 + edge (1..4).
    std::map<int,Intersection> intersections;
  };

  //! \brief Default constructor.
  ASMu2D(unsigned char n_s = 2, unsigned char n_f = 2);
  //! \brief Copy constructor.
  ASMu2D(const ASMu2D& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMu2D() { geomB = nullptr; }

  //! \brief Returns the spline surface representing a basis of this patch.
  virtual const LR::LRSplineSurface* getBasis(int basis = 1) const;
  //! \brief Returns the spline surface representing a basis of this patch.
  virtual LR::LRSplineSurface* getBasis(int basis = 1);


  // Methods for model generation and refinement
  // ===========================================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineSurface object to given stream.
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
  virtual void getElmConnectivities(IntMat& neighs) const;

  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const
  { this->getCoordinates(X, nsd, *lrspline); }

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ);

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary edge
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (for mixed methods)
  //! \param[in] orient Orientation of boundary (used for sorting)
  //! \param[in] local If \e true, return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes, int basis, int = 1,
                                int orient = 0, bool local = false) const;

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face/edge
  //! \param[in] orient Orientation of boundary (used for sorting)
  //! \param[out] elms Array of element numbers
  virtual void getBoundaryElms(int lIndex, int orient, IntVec& elms) const;

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction
  virtual bool getOrder(int& p1, int& p2, int& p3) const;

  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;

  //! \brief Returns the number of projection nodes for this patch.
  virtual size_t getNoProjectionNodes() const;

  //! \brief Refines along the diagonal of the LR-spline patch.
  //! \details Progressively refine until the LR-spline object contains at least
  //! \a minBasisfunctions basis functions.
  //! \param[in] minBasisfunctions lower bound on number of basis functions
  bool diagonalRefine(int minBasisfunctions);
  //! \brief Refines the lower-left corner of the LR-spline patch.
  //! \brief Progressively refine until the LR-spline object contains at least
  //! \a minBasisfunctions basis functions.
  //! \param[in] minBasisfunctions lower bound on number of basis functions
  bool cornerRefine(int minBasisfunctions);
  //! \brief Refines the LR-spline patch uniformly.
  //! \details Inserts one (global) line at a time until the LR-spline object
  //! contains at least \a minBasisfunctions basis functions.
  //! \param[in] minBasisfunctions lower bound on number of basis functions
  bool uniformRefine(int minBasisfunctions);
  //! \brief Refines the parametrization by inserting tensor knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int dir, int nInsert);
  using ASMLRSpline::refine;
  //! \brief Refines the parametrization by inserting extra tensor knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  //! \param[in] scale Scaling factor for the added knot values
  virtual bool refine(int dir, const RealArray& xi, double scale);
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol);
  //! \brief Raises the order of the tensor spline object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  virtual bool raiseOrder(int ru, int rv);

  //! \brief Creates a separate projection basis for this patch.
  virtual bool createProjectionBasis(bool init);

  //! \brief Sets the minimum element area for adaptive refinement.
  virtual void setMinimumSize(double size) { aMin = size; }
  //! \brief Defines the minimum element area for adaptive refinement.
  //! \param[in] nrefinements Maximum number of adaptive refinement levels
  virtual double getMinimumSize(int nrefinements) const;
  //! \brief Checks if the specified element is larger than the minimum size.
  //! \param[in] elmId Global/patch-local element index
  //! \param[in] globalNum If \e true, \a elmId is global otherwise patch-local
  virtual bool checkElementSize(int elmId, bool globalNum = true) const;


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] dir Parameter direction defining the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which DOFs to constrain at each node along the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain edge for
  virtual void constrainEdge(int dir, bool open, int dof, int code, char basis);
  //! \brief Constrains all DOFs in local directions on a given boundary edge.
  //! \param[in] dir Parameter direction defining the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which local DOFs to constrain at each node on the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] project If \e true, the local axis directions are projected
  //! \return Number of additional nodes added due to local axis constraints
  virtual size_t constrainEdgeLocal(int dir, bool open, int dof, int code,
                                    bool project = false);

  //! \brief Constrains a corner node identified by the two parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain node for
  //!
  //! \details The sign of the two indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  virtual void constrainCorner(int I, int J, int dof,
                               int code = 0, char basis = 1);
  //! \brief Constrains a node identified by two relative parameter values.
  //! \param[in] xi Parameter in u-direction
  //! \param[in] eta Parameter in v-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The parameter values have to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual void constrainNode(double xi, double eta, int dof, int code = 0);

  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  //! \param[in] coordCheck False to disable coordinate checks
  //! \param[in] thick Thickness of connection
  virtual bool connectPatch(int edge, ASM2D& neighbor, int nedge, bool revers,
                            int = 0, bool coordCheck = true, int thick = 1);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  using ASMbase::integrate;
  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates an integral over element interfaces in the patch.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] iChk Object checking if an element interface has contributions
  virtual bool integrate(Integrand& integrand, GlobalIntegral& glbInt,
                         const TimeDomain& time,
                         const ASM::InterfaceChecker& iChk);

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

  //! \brief Returns the node index for a given corner.
  //! \param[in] I -1 or +1 for either umin or umax corner
  //! \param[in] J -1 or +1 for either vmin or vmax corner
  //! \param[in] basis which basis to consider (for mixed methods)
  virtual int getCorner(int I, int J, int basis) const;

protected:
  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] itgPts Parameters and weights of the integration points
  bool integrate(Integrand& integrand, GlobalIntegral& glbInt,
                 const TimeDomain& time, const Real3DMat& itgPts);

public:

  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameters in range [0.0,1.0] of the point
  //! \param[out] param The (u,v) parameters of the point in knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that matches the point, if any
  //! \return 0 if no node (control point) matches this point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const;

  //! \brief Returns the element that contains a specified spatial point.
  //! \param[in] param The parameters of the point in the knot-span domain
  //! \return Local element number within the patch that contains the point
  virtual int findElementContaining(const double* param) const;

  //! \brief Calculates parameter values for visualization nodal points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  virtual bool getGridParameters(RealArray& prm, int dir, int nSegSpan) const;

  //! \brief Creates a quad element model of this patch for visualization.
  //! \param[out] grid The generated quadrilateral grid
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

  //! \brief Evaluates the projected solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] nf If nonzero, mixed evaluates nf fields on first basis
  virtual bool evalProjSolution(Matrix& sField, const Vector& locSol,
                                const int* npe, int nf) const;

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

  //! \brief Stores the mesh basis to encapsulated postscript files.
  //! \param[in] fName Prefix for file names
  //! \param[in] fType Flag telling which file type(s) to write
  virtual void storeMesh(const std::string& fName, int fType) const;

protected:
  //! \brief Struct representing an inhomogeneous Dirichlet boundary condition.
  struct DirichletEdge
  {
    LR::LRSplineSurface* lr;       //!< Pointer to the right object (in case of multiple bases)
    LR::parameterEdge   edg;       //!< Which edge is this
    IntVec              MLGE;      //!< Local-to-Global Element numbers
    IntVec              MLGN;      //!< Local-to-Global Nodal numbers
    IntMat              MNPC;      //!< Matrix of Nodal-Point Correpondanse
    int                 dof;       //!< Local DOF to constrain along the boundary
    int                 code;      //!< Inhomogeneous Dirichlet condition code
    int                 corners[2];//!< Index of the two end-points of this line

    //! \brief The constructor detects the edge end points.
    DirichletEdge(LR::LRSplineSurface* sf, int dir,
                  int d = 0, int c = 0, int offset = 1);
    //! \brief Returns \e true if basis function \a b is at a corner point.
    bool isCorner(int b) const { return b == corners[0] || b == corners[1]; }
  };

  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  LR::LRSplineSurface* projectSolution(const IntegrandBase& integrand) const;
  //! \brief Projects the secondary solution using a superconvergent approach.
  //! \param[in] integrand Object with problem-specific data and methods
  LR::LRSplineSurface* scRecovery(const IntegrandBase& integrand) const;

  //! \brief Interpolates an LR spline at a given set of parametric coordinates.
  //! \param[in] upar Parametric interpolation points in the u-direction
  //! \param[in] vpar Parametric interpolation points in the v-direction
  //! \param[in] points Interpolation values stored as one point per matrix row
  //! \param[in] basis Which basis to interpolate onto
  //! \return A LRSplineSurface representation of the interpolated points
  LR::LRSplineSurface* regularInterpolation(const RealArray& upar,
                                            const RealArray& vpar,
                                            const Matrix& points, int basis) const;

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
  //! the \a u and \a v parameters directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const RealArray* gpar, bool = false) const;

  //! \brief Projects inhomogenuous dirichlet conditions by continuous L2-fit.
  //! \param[in] edge low-level edge information needed to do integration
  //! \param[in] values inhomogenuous function which is to be fitted
  //! \param[out] result fitted value in terms of control-point values
  //! \param[in] time time used in dynamic problems
  bool edgeL2projection(const DirichletEdge& edge, const FunctionBase& values,
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

  //! \brief Assembles L2-projection matrices for the secondary solution.
  //! \param[out] A Left-hand-side matrix
  //! \param[out] B Right-hand-side vectors
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e false, a discrete L2-projection is used
  virtual bool assembleL2matrices(SparseMatrix& A, StdVector& B,
                                  const L2Integrand& integrand,
                                  bool continuous) const;

  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  //! \param[in] basis Which basis to connect the nodes for (mixed methods)
  //! \param[in] slave 0-based index of the first slave node in this basis
  //! \param[in] master 0-based index of the first master node in this basis
  //! \param[in] coordCheck False to disable coordinate checks
  //! \param[in] thick Thickness of connection
  bool connectBasis(int edge, ASMu2D& neighbor, int nedge, bool revers,
                    int basis = 1, int slave = 0, int master = 0,
                    bool coordCheck = true, int thick = 1);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] iel 1-based element index
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  void getGaussPointParameters(RealArray& uGP, int dir, int nGauss,
                               int iel, const double* xi) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] basisNum Which basis to get Greville point parameters for
  bool getGrevilleParameters(RealArray& prm, int dir, int basisNum = 1) const;

  //! \brief Returns the area in the parameter space for an element.
  //! \param[in] iel 1-based element index
  double getParametricArea(int iel) const;
  //! \brief Returns boundary edge length in the parameter space for an element.
  //! \param[in] iel 1-based element index
  //! \param[in] dir Local index of the boundary edge
  double getParametricLength(int iel, int dir) const;

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

  //! \brief Returns the node indices for a given edge.
  //! \param nodes Array of node numbers
  //! \param[in] edge Local index of the boundary edge
  //! \param[in] basis Which basis to grab nodes for
  //! \param[in] orient Orientation of boundary (used for sorting)
  //! \param[in] local If \e true, return patch-local node numbers
  void getEdgeNodes(IntVec& nodes, int edge, int basis,
                    int orient, bool local) const;

  //! \brief Evaluates the basis functions and derivatives of an element.
  //! \param[in] iel 0-based element index
  //! \param fe Integration point data for current element
  //! \param[in] derivs Derivative order of the basis functions
  bool evaluateBasis(int iel, FiniteElement& fe, int derivs = 0) const;

  //! \brief Evaluate basis functions in a point.
  void computeBasis(double u, double v,
                    Go::BasisPtsSf& bas, int iel,
                    const LR::LRSplineSurface* spline = nullptr) const;
  //! \brief Evaluate basis functions and first derivatives in a point.
  void computeBasis(double u, double v,
                    Go::BasisDerivsSf& bas, int iel,
                    const LR::LRSplineSurface* spline = nullptr) const;
  //! \brief Evaluate basis functions and two derivatives in a point.
  void computeBasis(double u, double v,
                    Go::BasisDerivsSf2& bas, int iel,
                    const LR::LRSplineSurface* spline = nullptr) const;
  //! \brief Evaluate basis functions and three derivatives in a point.
  void computeBasis(double u, double v,
                    Go::BasisDerivsSf3& bas, int iel,
                    const LR::LRSplineSurface* spline = nullptr) const;

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] iel 0-based local element index
  //! \param[in] param The (u,v) parameters of the point in knot-span domain
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
  std::shared_ptr<LR::LRSplineSurface> createLRfromTensor();

  //! \brief Converts a rational spline surface to a LR NURBS surface.
  static std::shared_ptr<LR::LRSplineSurface> createLRNurbs(const Go::SplineSurface& srf);

  //! \brief Generate bezier basis.
  void generateBezierBasis();

  //! \brief Generate bezier extraction operators.
  void generateBezierExtraction();

  //! \brief Write NURBS elements as postscript file.
  void writePostscriptElementsNurbs (std::shared_ptr<LR::LRSplineSurface> mesh,
                                     std::ostream& out, bool close = true,
                                     int nu = 2, int nv = 2);
  //! \brief Write NURBS elements as postscript file.
  void writePostscriptMeshWithControlPointsNurbs (std::shared_ptr<LR::LRSplineSurface> mesh,
                                                  std::ostream& out, int nu = 2, int nv = 2);

public:
  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

  //! \brief Evaluate NURBS basis functions in a point.
  static void computeBasisNurbs(double u, double v,
                                Go::BasisPtsSf& bas, int iel,
                                const LR::LRSplineSurface& spline);
  //! \brief Evaluate NURBS basis functions and first derivatives in a point.
  static void computeBasisNurbs(double u, double v,
                                Go::BasisDerivsSf& bas, int iel,
                                const LR::LRSplineSurface& spline);
  //! \brief Evaluate NURBS basis functions and two derivatives in a point.
  static void computeBasisNurbs(double u, double v,
                                Go::BasisDerivsSf2& bas, int iel,
                                const LR::LRSplineSurface& spline);
  //! \brief Evaluate NURBS basis functions and three derivatives in a point.
  static void computeBasisNurbs(double u, double v,
                                Go::BasisDerivsSf3& bas, int iel,
                                const LR::LRSplineSurface& spline);

  //! \brief Query whether basis is rational or not.
  bool rational() const { return is_rational; }

  //! \brief Returns a matrix with control point coordinates.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of points
  //! \param[in] nsd Number of spatial dimensions
  //! \param[in] spline Spline to extract coefficients from
  //! \param[in] iel 1-based element index, if -1 return for all control points
  static bool getCoordinates(Matrix& X,
                             unsigned char nsd,
                             const LR::LRSplineSurface& spline,
                             int iel = -1);

protected:
  std::shared_ptr<LR::LRSplineSurface> lrspline; //!< Pointer to the LR-spline surface object

  bool is_rational = false; //!< True if basis is rational

  Go::SplineSurface* tensorspline; //!< Pointer to original tensor spline object
  Go::SplineSurface* tensorPrjBas; //!< Pointer to tensor spline projection base
  // The tensor spline object is kept for backward compatability with the REFINE
  // and RAISEORDER key-words, although we take note that there is a possibility
  // of optimization since all mapping values and Jacobians may be performed on
  // this object for increased efficiency.

  //! Inhomogeneous Dirichlet boundary condition data
  std::vector<DirichletEdge> dirich;

  ThreadGroups threadGroups; //!< Element groups for multi-threaded assembly
  ThreadGroups projThreadGroups; //!< Element groups for multi-threaded assembly - projection basis

  const Matrices& bezierExtract; //!< Bezier extraction matrices
  Matrices      myBezierExtract; //!< Bezier extraction matrices

  Go::BsplineBasis bezier_u; //!< Bezier basis in the u-direction
  Go::BsplineBasis bezier_v; //!< Bezier basis in the v-direction

  //! Basis function cache
  std::vector<std::unique_ptr<BasisFunctionCache>> myCache;

private:
  mutable double aMin; //!< Minimum element area for adaptive refinement

  //! \brief Evaluates the NURBS basis functions and derivatives of an element.
  //! \param[in] iel 0-based element index
  //! \param fe Integration point data for current element
  //! \param[in] derivs Derivative order of the basis functions
  bool evaluateBasisNurbs(int iel, FiniteElement& fe,
                          int derivs) const;
};

#endif
