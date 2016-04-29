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

#include "ASMunstruct.h"
#include "ASM2D.h"
#include <memory>

class FiniteElement;

namespace Go {
  class SplineCurve;
  class SplineSurface;
}

namespace LR {
  class LRSplineSurface;
}


/*!
  \brief Driver for assembly of unstructured 2D spline FE models.
  \details This class contains methods common for 2D LR-spline patches.
*/

class ASMu2D : public ASMunstruct, public ASM2D
{
public:
  //! \brief Default constructor.
  ASMu2D(unsigned char n_s = 2, unsigned char n_f = 2);
  //! \brief Copy constructor.
  ASMu2D(const ASMu2D& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMu2D() { geo = nullptr; }

  //! \brief Returns the spline surface representing this patch.
  LR::LRSplineSurface* getSurface() { return lrspline.get(); }

  //! \brief Returns the spline surface representing the basis of this patch.
  virtual const LR::LRSplineSurface* getBasis(int = 1) const { return lrspline.get(); }
  //! \brief Returns the spline surface representing the basis of this patch.
  virtual LR::LRSplineSurface* getBasis(int = 1) { return lrspline.get(); }

  // Methods for model generation and refinement
  // ===========================================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineSurface object to given stream.
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
  //! \param[in] iel Element index
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

  //! \brief Finds the global numbers of the nodes on a patch boundary.
  //! \param[in] lIndex Local index of the boundary edge
  //! \param glbNodes Array of global boundary node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& glbNodes, int basis = 0) const;

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction
  virtual bool getOrder(int& p1, int& p2, int& p3) const;

  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;

  //! \brief Checks that the patch is modelled in a right-hand-side system.
  virtual bool checkRightHandSystem();

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
  //! \details This method is mainly kept for backward compatability with the
  //! "REFINE" keyword in the input file.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int dir, int nInsert);
  using ASMunstruct::refine;
  //! \brief Refines the parametrization by inserting extra tensor knots.
  //! \details This method is mainly kept for backward compatability with the
  //! "REFINE" keyword in the input file.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  virtual bool refine(int dir, const RealArray& xi);
  //! \brief Refines the parametrization based on a mesh density function.
  //! \param[in] refC Mesh refinement criteria function
  //! \param[in] refTol Mesh refinement threshold
  virtual bool refine(const RealFunc& refC, double refTol);
  //! \brief Raises the order of the tensor spline object for this patch.
  //! \details This method is mainly kept for backward compatability with the
  //! "RAISEORDER" keyword in the input file.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  virtual bool raiseOrder(int ru, int rv);


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] dir Parameter direction defining the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which DOFs to constrain at each node on the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  virtual void constrainEdge(int dir, bool open, int dof, int code, char);
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
  //!
  //! \details The sign of the two indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  virtual void constrainCorner(int I, int J, int dof, int code = 0, char = 1);
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

  /* More multipatch stuff, maybe later...
  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  virtual bool connectPatch(int edge, ASMu2D& neighbor, int nedge,
                            bool revers = false);

  //! \brief Makes two opposite boundary edges periodic.
  //! \param[in] dir Parameter direction defining the periodic edges
  //! \param[in] basis Which basis to connect (mixed methods), 0 means both
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeEdges(int dir, int basis = 0, int master = 1);
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

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  //! \param[in] time Current time
  //! \param[in] g2l Pointer to global-to-local node number mapping
  virtual bool updateDirichlet(const std::map<int,RealFunc*>& func,
                               const std::map<int,VecFunc*>& vfunc, double time,
                               const std::map<int,int>* g2l = nullptr);

  //! \brief Returns the node index for a given corner.
  int getCorner(int I, int J, int basis = 1) const;
  //! \brief Rerturns the node indices for a given edge.
  std::vector<int> getEdgeNodes(int dir, int basis = 1) const;

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
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int* npe) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //! \param[in] deriv Derivative order to return
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size().
  //! Otherwise, we assume that it contains the \a u and \a v parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool regular = false,
                            int deriv = 0) const;

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

private:
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
  //! \return A LRSplineSurface representation of the interpolated points
  LR::LRSplineSurface* regularInterpolation(const RealArray& upar,
                                            const RealArray& vpar,
                                            const Matrix& points) const;

public:
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual LR::LRSpline* evalSolution(const IntegrandBase& integrand) const;

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
  //! \a gpar[0].size() \a X \a gpar[1].size().
  //! Otherwise, we assume that it contains the \a u and \a v parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const RealArray* gpar, bool regular = true) const;

  //! \brief Projects the secondary solution using a discrete global L2-norm.
  //! \param[out] sField Secondary solution field control point values
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e true, a continuous L2-projection is used
  virtual bool globalL2projection(Matrix& sField,
				  const IntegrandBase& integrand,
				  bool continuous = false) const;

  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVars Gauss point variables associated with \a oldBasis
  //! \param[out] newVars Gauss point variables associated with this patch.
  //! \param[in] nGauss Number of Gauss points along a knot-span
  bool transferGaussPtVars(const LR::LRSplineSurface* oldBasis,
                           const RealArray& oldVars, RealArray& newVars,
                           int nGauss) const;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVars Control point variables associated with \a oldBasis
  //! \param[out] newVars Gauss point variables associated with this patch.
  //! \param[in] nGauss Number of Gauss points along a knot-span
  bool transferCntrlPtVars(LR::LRSplineSurface* oldBasis,
                           const RealArray& oldVars, RealArray& newVars,
                           int nGauss) const;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[out] newVars Gauss point variables associated with this patch.
  //! \param[in] nGauss Number of Gauss points along a knot-span
  bool transferCntrlPtVars(const LR::LRSplineSurface* oldBasis,
                           RealArray& newVars, int nGauss) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  //! \param[in] basis Which basis to connect the nodes for (mixed methods)
  //! \param[in] slave 0-based index of the first slave node in this basis
  //! \param[in] master 0-based index of the first master node in this basis
  bool connectBasis(int edge, ASMu2D& neighbor, int nedge, bool revers,
                    int basis = 1, int slave = 0, int master = 0);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] iel Element index
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  void getGaussPointParameters(RealArray& uGP, int dir, int nGauss,
                               int iel, const double* xi) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  bool getGrevilleParameters(RealArray& prm, int dir) const;

  //! \brief Returns the area in the parameter space for an element.
  //! \param[in] iel Element index
  double getParametricArea(int iel) const;
  //! \brief Returns boundary edge length in the parameter space for an element.
  //! \param[in] iel Element index
  //! \param[in] dir Local index of the boundary edge
  double getParametricLength(int iel, int dir) const;

  //! \brief Computes the element corner coordinates.
  //! \param[in] iel Element index
  //! \param[out] XC Coordinates of the element corners
  void getElementCorners(int iel, std::vector<Vec3>& XC) const;

  //! \brief Evaluates the basis functions and derivatives of order \a derivs
  //! of an element.
  bool evaluateBasis(FiniteElement& el, int derivs = 0) const;

public:
  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

  typedef std::pair<int,int> Ipair; //!< Convenience type

protected:
  //! \brief Struct representing an inhomogeneous Dirichlet boundary condition.
  struct DirichletEdge
  {
    Go::SplineCurve*   curve; //!< Pointer to spline curve for the boundary
    int                dof;   //!< Local DOF to constrain along the boundary
    int                code;  //!< Inhomogeneous Dirichlet condition code
    std::vector<Ipair> nodes; //!< Nodes subjected to projection on the boundary
    //! \brief Default constructor.
    DirichletEdge(Go::SplineCurve* sc = nullptr, int d = 0, int c = 0)
    : curve(sc), dof(d), code(c) {}
  };

  std::shared_ptr<LR::LRSplineSurface> lrspline; //!< Pointer to the LR-spline surface object

  Go::SplineSurface* tensorspline; //!< Pointer to original tensor spline object
  // The tensor spline object is kept for backward compatability with the REFINE
  // and RAISEORDER key-words, although we take note that there is a possibility
  // of optimization since all mapping values and Jacobians may be performed on
  // this object for increased efficiency.

  //! Inhomogeneous Dirichlet boundary condition data
  std::vector<DirichletEdge> dirich;

  const std::vector<Matrix>& bezierExtract; //!< Bezier extraction matrices
  std::vector<Matrix>      myBezierExtract; //!< Bezier extraction matrices

  Go::BsplineBasis bezier_u;
  Go::BsplineBasis bezier_v;
};

#endif
