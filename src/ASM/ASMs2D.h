// $Id$
//==============================================================================
//!
//! \file ASMs2D.h
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D spline FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_H
#define _ASM_S2D_H

#include "ASMstruct.h"
#include "ASM2D.h"
#include "BasisFunctionCache.h"
#include "Interface.h"
#include "ThreadGroups.h"

#include <memory>

namespace utl {
  class Point;
}

namespace Go {
  class SplineCurve;
  class SplineSurface;
}


/*!
  \brief Driver for assembly of structured 2D spline FE models.
  \details This class contains methods common for structured 2D spline patches.
*/

class ASMs2D : public ASMstruct, public ASM2D
{
  //! \brief Struct for nodal point data.
  struct IJ
  {
    int I; //!< Index in first parameter direction
    int J; //!< Index in second parameter direction
  };

  typedef std::vector<IJ> IndexVec; //!< Node index container

  //! \brief Struct for edge node definitions.
  struct Edge
  {
    int icnod; //!< Global node number of first interior point along the edge
    int incr;  //!< Increment in the global numbering along the edge

    //! \brief Default constructor.
    Edge() { icnod = incr = 0; }
    //! \brief Returns \a icnod which then is incremented.
    int next();
  };

protected:
  //! \brief Implementation of basis function cache.
  class BasisFunctionCache : public ::BasisFunctionCache<2>
  {
  public:
    //! \brief The constructor initializes the class.
    //! \param pch Patch the cache is for
    BasisFunctionCache(const ASMs2D& pch);

    //! \brief Constructor reusing quadrature info from another instance.
    //! \param cache Instance holding quadrature information
    //! \param b Basis to use
    BasisFunctionCache(const BasisFunctionCache& cache, int b);

    //! \brief Empty destructor.
    virtual ~BasisFunctionCache() = default;

    //! \brief Returns number of elements in each direction.
    const std::array<size_t,2>& noElms() const { return nel; }

  protected:
    //! \brief Implementation specific initialization.
    bool internalInit() override;

    //! \brief Calculates basis function info in a single integration point.
    //! \param el Element of integration point (0-indexed)
    //! \param gp Integration point on element (0-indexed)
    //! \param reduced If true, returns values for reduced integration scheme
    BasisFunctionVals calculatePt(size_t el, size_t gp, bool reduced) const override;

    //! \brief Calculates basis function info in all integration points.
    void calculateAll() override;

    //! \brief Obtain global integration point index.
    //! \param el Element of integration point (0-indexed)
    //! \param gp Integration point on element (0-indexed)
    //! \param reduced True to return index for reduced quadrature
    size_t index(size_t el, size_t gp, bool reduced) const override;

    //! \brief Setup integration point parameters.
    virtual void setupParameters();

    //! \brief Configure quadratures.
    bool setupQuadrature();

    const ASMs2D& patch; //!< Reference to patch cache is for

    std::array<size_t,2> nel{}; //!< Number of elements in each direction

  private:
    //! \brief Obtain structured element indices.
    //! \param el Global element index
    std::array<size_t,2> elmIndex(size_t el) const;
  };

public:
  //! \brief Struct with data for definition of global node numbers of a patch.
  struct BlockNodes
  {
    int  ibnod[4]; //!< Vertex nodes
    Edge edges[4]; //!< Edge nodes
    int  iinod;    //!< Global node number of the first interior node
    int  inc[2];   //!< Increment in global node numbering in each direction
    int  nnodI;    //!< Number of nodes in parameter direction I
    int  indxI;    //!< Running node index in the local I-direction

    //! \brief Default constructor.
    BlockNodes()
    {
      memset(ibnod,0,4*sizeof(int));
      iinod = inc[0] = inc[1] = nnodI = 0;
      indxI = 1;
    }
    //! \brief Returns \a iinod which then is incremented.
    int next();
  };

private:
  typedef std::pair<int,int> Ipair; //!< Convenience type

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

public:
  //! \brief Base class that checks if an element has interface contributions.
  class InterfaceChecker : public ASM::InterfaceChecker
  {
  protected:
    const ASMs2D& myPatch; //!< Reference to the patch being integrated
  public:
    //! \brief The constructor initialises the reference to current patch.
    explicit InterfaceChecker(const ASMs2D& pch) : myPatch(pch) {}
    //! \brief Empty destructor.
    virtual ~InterfaceChecker() {}
    //! \brief Returns non-zero if the specified element have contributions.
    //! \param[in] I Index in first parameter direction of the element
    //! \param[in] J Index in second parameter direction of the element
    virtual short int hasContribution(int, int I, int J, int = -1) const;
  };

  //! \brief Default constructor.
  ASMs2D(unsigned char n_s = 2, unsigned char n_f = 2);
  //! \brief Special copy constructor for sharing of FE data.
  ASMs2D(const ASMs2D& patch, unsigned char n_f);
  //! \brief Default copy constructor copying everything.
  ASMs2D(const ASMs2D& patch);
  //! \brief The destructor frees the dynamically allocated boundary curves.
  virtual ~ASMs2D();

  //! \brief Returns the spline curve representing a boundary of this patch.
  //! \param[in] dir Parameter direction defining which boundary to return
  virtual Go::SplineCurve* getBoundary(int dir, int = 1);
  //! \brief Returns the spline surface representing a basis of this patch.
  virtual Go::SplineSurface* getBasis(int basis = 1);
  //! \brief Returns the spline surface representing a basis of this patch.
  virtual const Go::SplineSurface* getBasis(int basis = 1) const;
  //! \brief Copies the parameter domain from the \a other patch.
  virtual void copyParameterDomain(const ASMbase* other);

  // Methods for model generation
  // ============================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineSurface object to given stream.
  virtual bool write(std::ostream&, int) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the node-to-IJ-index array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Adds extraordinary elements associated with a patch boundary.
  //! \param[in] dim Dimension of the boundary (should be 1)
  //! \param[in] item Local index of the boundary edge
  //! \param[in] nXn Number of extraordinary nodes
  //! \param[out] nodes Global numbers assigned to the extraordinary nodes
  virtual bool addXElms(short int dim, short int item,
                        size_t nXn, IntVec& nodes);

  //! \brief Adds interface elements with coupling to all element DOFs.
  //! \param[in] iChk Object checking if an element interface has contributions
  bool addInterfaceElms(const ASM::InterfaceChecker& iChk);

  //! \brief Returns local 1-based index of the node with given global number.
  //! \details If the given node number is not present, 0 is returned.
  //! \param[in] globalNum Global node number
  //! \param[in] noAddedNodes If \e true, use \a xnMap to find the real node
  virtual size_t getNodeIndex(int globalNum, bool noAddedNodes = false) const;
  //! \brief Returns the global node number for the given node.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] noAddedNodes If \e true, use \a nxMap to find the real node
  virtual int getNodeID(size_t inod, bool noAddedNodes = false) const;

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  //! \param[in] iel 1-based element index local to current patch
  //! \param[in] forceItg If \e true, return the integration basis coordinates
  //! otherwise the geometry basis coordinates are returned
  virtual bool getElementCoordinates(Matrix& X, int iel, bool forceItg = false) const;

  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! \param[in] geo If \e true, coordinates for the geometry basis are returned
  //! otherwise the integration basis coordinates are returned
  virtual void getNodalCoordinates(Matrix& X, bool geo = false) const;

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Returns the coordinate of the element center.
  virtual Vec3 getElementCenter(int iel) const;

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ);

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary edge
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (for mixed methods)
  //! \param[in] thick Thickness of connection
  //! \param[in] local If \e true, return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int basis, int thick = 1,
                                int = 0, bool local = false) const;

  //! \brief Returns the node index for a given corner.
  //! \param[in] I -1 or +1 for either umin or umax corner
  //! \param[in] J -1 or +1 for either vmin or vmax corner
  //! \param[in] basis which basis to consider (for mixed methods)
  virtual int getCorner(int I, int J, int basis) const;

  //! \brief Assigns new global node numbers for all nodes of the patch.
  //! \param nodes Object with global nodes numbers to assign to this patch
  //! \param[in] basis Which basis to assign node numbers for in mixed methods
  //!
  //! \details The global node numbers generated by generateFEMTopology() are
  //! non-unique in the sense that a node that is shared by two (or more)
  //! patches along a common interface has a different number in each patch.
  //! This method therefore assigns a new global number to each node in the
  //! patch. The data provided through the \a nodes argument is sufficient
  //! to determine the unique global number under the assumption that they
  //! are ordered in the sequence determined by the local orientation of the
  //! patch and its edges.
  bool assignNodeNumbers(BlockNodes& nodes, int basis = 0);

  //! \brief Checks that the patch is modelled in a right-hand-side system.
  //! \details If it isn't, the v-parameter direction is swapped.
  virtual bool checkRightHandSystem();

  //! \brief Refines the parametrization by inserting extra knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  //! \param[in] scale Scaling factor for the added knot values
  virtual bool refine(int dir, const RealArray& xi, double scale = 1.0);
  //! \brief Refines the parametrization by inserting extra knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int dir, int nInsert);
  //! \brief Raises the order of the SplineSurface object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  virtual bool raiseOrder(int ru, int rv);

  //! \brief Creates a separate projection basis for this patch.
  virtual bool createProjectionBasis(bool init);

  //! \brief Checks if a separate projection basis is used for this patch.
  virtual bool separateProjectionBasis() const;


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
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //! \param[in] thick Thickness of connection
  virtual bool connectPatch(int edge, ASM2D& neighbor, int nedge, bool revers,
                            int = 0, bool coordCheck = true, int thick = 1);

  //! \brief Makes two opposite boundary edges periodic.
  //! \param[in] dir Parameter direction defining the periodic edges
  //! \param[in] basis Which basis to connect (mixed methods)
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeBoundaries(int dir, int basis, int master);

  //! \brief Collapses a degenereated edge into a single node.
  //! \param[in] dir Parameter direction defining the edge to collapse
  //! \param[in] basis Which basis to collapse edge for
  virtual bool collapseEdge(int dir, int basis = 1);

  //! \brief Adds MPCs representing a rigid coupling to this patch.
  //! \param[in] lindx Local index of the boundary item that should be rigid
  //! \param[in] ldim Dimension of the boundary item that should be rigid
  //! \param[in] basis Which basis to add rigid coupling for (mixed methods)
  //! \param gMaster Global node number of the master node
  //! \param[in] Xmaster Position of the master nodal point
  //! \param[in] extraPt If \e true, the master point is not a patch node
  //! \return \e true if a new global node was added, otherwise \e false
  virtual bool addRigidCpl(int lindx, int ldim, int basis,
                           int& gMaster, const Vec3& Xmaster, bool extraPt);

  //! \brief Sets the global node numbers for this patch.
  //! \param[in] nodes Vector of global node numbers (zero-based)
  virtual void setNodeNumbers(const IntVec& nodes);

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  //! \param[in] time Current time
  //! \param[in] g2l Pointer to global-to-local node number mapping
  virtual bool updateDirichlet(const std::map<int,RealFunc*>& func,
                               const std::map<int,VecFunc*>& vfunc, double time,
                               const std::map<int,int>* g2l = nullptr);


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
  //! \param[in] lIndex Local index [1,4] of the boundary edge
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
                         const TimeDomain& time, const ASM::InterfaceChecker& iChk);

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

  //! \brief Searches for the specified Cartesian point in the patch.
  //! \param X The Cartesian coordinates of the point, updated on output
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \return Distance from the point \a X to the found point
  virtual double findPoint(Vec3& X, double* param) const;

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
  //! \param[in] n_f If nonzero, mixed evaluates \a n_f fields on first basis
  //! \param[in] piola If \e true, use piola mapping
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int* npe, int n_f, bool piola) const;

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
                            const RealArray* gpar, bool regular = true,
                            int deriv = 0, int = 0) const;

  //! \brief Evaluates the projected solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] n_f If nonzero, mixed evaluates \a n_f fields on first basis
  virtual bool evalProjSolution(Matrix& sField, const Vector& locSol,
                                const int* npe, int n_f) const;

  //! \brief Evaluates and interpolates a field over a given geometry.
  //! \param[in] basis The basis of the field to evaluate
  //! \param[in] locVec The coefficients of the field to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum Basis number (mixed)
  virtual bool evaluate(const ASMbase* basis, const Vector& locVec,
                        RealArray& vec, int basisNum) const;

  //! \brief Evaluates and interpolates a field over a given geometry.
  //! \param[in] field The field to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum Basis number (mixed)
  virtual bool evaluate(const Field* field, RealArray& vec, int basisNum) const;

  //! \brief Evaluates and interpolates a function over a given geometry.
  //! \param[in] func The function to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum Basis number (mixed)
  //! \param[in] time Current time
  virtual bool evaluate(const FunctionBase* func, RealArray& vec,
                        int basisNum, double time) const;

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
                            const int* npe = nullptr, char project = 0) const;

private:
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  Go::SplineSurface* projectSolution(const IntegrandBase& integrand) const;

  //! \brief Projects the secondary solution using a superconvergent approach.
  Go::SplineSurface* scRecovery(const IntegrandBase&) const;

  //! \brief Projects the secondary solution using Quasi-Interpolation.
  Go::SplineSurface* projectSolutionLocal(const IntegrandBase&) const;

  //! \brief Projects the secondary solution using
  //! Variation Diminishing Spline Approximation.
  Go::SplineSurface* projectSolutionLocalApprox(const IntegrandBase&) const;

  //! \brief Projects the secondary solution using Least-Square Approximation.
  Go::SplineSurface* projectSolutionLeastSquare(const IntegrandBase&) const;

public:
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual Go::GeomObject* evalSolution(const IntegrandBase& integrand) const;

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

protected:

  // Internal utility methods
  // ========================

  //! \brief Finds the patch-local element numbers on a patch boundary.
  //! \param[out] elms Array of element numbers
  //! \param[in] lIndex Local index of the boundary edge
  virtual void findBoundaryElms(IntVec& elms, int lIndex, int = 0) const;

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
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //! \param[in] thick Thickness of connection
  bool connectBasis(int edge, ASMs2D& neighbor, int nedge, bool revers,
                    int basis = 1, int slave = 0, int master = 0,
                    bool coordCheck = true, int thick = 1);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  void getGaussPointParameters(RealArray& uGP, int dir,
                               int nGauss, const double* xi) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] basisNum Which basis to get Greville point parameters for
  bool getGrevilleParameters(RealArray& prm, int dir, int basisNum = 1) const;

  //! \brief Calculates parameter values for the Quasi-Interpolation points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  bool getQuasiInterplParameters(RealArray& prm, int dir) const;

  //! \brief Returns the area in the parameter space for an element.
  //! \param[in] iel Element index
  double getParametricArea(int iel) const;
  //! \brief Returns boundary edge length in the parameter space for an element.
  //! \param[in] iel Element index
  //! \param[in] dir Local index of the boundary edge
  double getParametricLength(int iel, int dir) const;

  //! \brief Returns a matrix with nodal coordinates for element containing given parameters.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! \param[in] u First parameter of point
  //! \param[in] v Second parameter of point
  bool getElementCoordinatesPrm(Matrix& X, double u, double v) const;

  //! \brief Computes the element border parameters.
  //! \param[in] i1 Parameter index in u-direction
  //! \param[in] i2 Parameter index in v-direction
  //! \param[out] u Parameter values of the west-east borders
  //! \param[out] v Parameter values of the south-north borders
  void getElementBorders(int i1, int i2, double* u, double* v) const;
  //! \brief Computes the element border parameters.
  //! \param[in] iel 1-based element index
  //! \param[out] u Parameter values of the element borders
  virtual void getElementBorders(int iel, double* u) const;

  //! \brief Computes the element corner coordinates.
  //! \param[in] i1 Parameter index in u-direction
  //! \param[in] i2 Parameter index in v-direction
  //! \param[out] XC Coordinates of the element corners
  //! \param[out] uC Spline parameters of the element corners (optional)
  //! \return Characteristic element size
  double getElementCorners(int i1, int i2, std::vector<Vec3>& XC,
                           RealArray* uC = nullptr) const;
  //! \brief Computes the element corner coordinates and parameters.
  //! \param[in] i1 Parameter index in u-direction
  //! \param[in] i2 Parameter index in v-direction
  //! \param[out] XC Coordinates and parameters of the element corners
  void getCornerPoints(int i1, int i2, std::vector<utl::Point>& XC) const;

  using ASMbase::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM Sanity check option
  virtual void generateThreadGroups(const Integrand& integrand, bool silence,
                                    bool ignoreGlobalLM);

  //! \brief Hook for changing number of threads.
  virtual void changeNumThreads();

  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] strip1 Strip width in first direction
  //! \param[in] strip2 Strip width in second direction
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM Sanity check option
  void generateThreadGroups(size_t strip1, size_t strip2,
                            bool silence, bool ignoreGlobalLM);

  //! \brief Generates element groups from a partition.
  virtual void generateThreadGroupsFromElms(const IntVec& elms);

  //! \brief Returns 0-based index of first node on integration basis.
  virtual int getFirstItgElmNode() const { return 0; }
  //! \brief Returns 0-based index of last node on integration basis.
  virtual int getLastItgElmNode() const;

public:
  //! \brief Auxilliary function for computation of basis function indices.
  static void scatterInd(int n1, int n2, int p1, int p2,
			 const int* start, IntVec& index);

private:
  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  bool getOrder(int& p1, int& p2) const;

public:
  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction (always zero)
  virtual bool getOrder(int& p1, int& p2, int& p3) const;

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[in] basis Which basis to return size parameters for (mixed methods)
  virtual bool getSize(int& n1, int& n2, int basis = 0) const;

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  //! \param[in] basis Which basis to return size parameters for (mixed methods)
  virtual bool getSize(int& n1, int& n2, int& n3, int basis) const;

  //! \brief Returns the number of elements in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  virtual bool getNoStructElms(int& n1, int& n2, int& n3) const;

  //! \brief Obtain element neighbours.
  virtual void getElmConnectivities(IntMat& neigh, bool local = false) const;

  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;
  //! \brief Returns the number of projection nodes for this patch.
  virtual size_t getNoProjectionNodes() const;

  //! \brief Returns parameter values and node numbers of the domain corners.
  //! \param[out] u Parameter values of the domain corners
  //! \param[out] corners 1-based indices of the corner nodes (optional)
  virtual bool getParameterDomain(Real2DMat& u, IntVec* corners) const;

  //! \brief Evaluates the basis functions at the specified point.
  //! \param[in] u First parameter value of evaluation point
  //! \param[in] v Second parameter value of evaluation point
  //! \param[out] N Basis function values
  virtual void evaluateBasis(double u, double v, double, Vector& N) const;

  //! \brief Establishes matrices with basis functions and 1st derivatives.
  //! \param[in] u First parameter value of current integration point
  //! \param[in] v Second parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  //! \param[in] fromRight If \e true, evaluate from right if at a knot
  void extractBasis(double u, double v, Vector& N, Matrix& dNdu,
                    bool fromRight = true) const;
  //! \brief Establishes matrices with basis functions, 1st and 2nd derivatives.
  //! \param[in] u First parameter value of current integration point
  //! \param[in] v Second parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  //! \param[out] d2Ndu2 Second derivatives of basis functions
  //! \param[in] fromRight If \e true, evaluate from right if at a knot
  void extractBasis(double u, double v, Vector& N,
                    Matrix& dNdu, Matrix3D& d2Ndu2,
                    bool fromRight = true) const;
  //! \brief Establishes a vector with basis function derivatives.
  //! \param[in] u First parameter value of current integration point
  //! \param[in] v Second parameter value of current integration point
  //! \param[in] dir Which parameter to establish derivatives with respect to
  //! \param[in] p The derivation order
  //! \param[out] dN Basis function derivatives
  //! \param[in] fromRight If \e true, evaluate from right if at a knot
  void extractBasis(double u, double v, int dir, int p, Vector& dN,
                    bool fromRight = true) const;

  //! \brief Returns a field using the projection basis.
  //! \param[in] coefs The coefficients for the field
  virtual Field* getProjectedField(const Vector& coefs) const;

  //! \brief Returns a field using the projection basis.
  //! \param[in] coefs The coefficients for the field
  virtual Fields* getProjectedFields(const Vector& coefs, size_t = 0) const;

private:
  //! \brief Returns an index into the internal coefficient array for a node.
  //! \param[in] inod 0-based node index local to current patch
  int coeffInd(size_t inod) const;

protected:
  std::shared_ptr<Go::SplineSurface> surf; //!< The actual spline surface object
  Go::SplineCurve* bou[4]; //!< Pointers to the four boundary curves
  bool              swapV; //!< Has the v-parameter direction been swapped?

  const IndexVec& nodeInd; //!< IJ-pairs for the control points (nodes)
  IndexVec      myNodeInd; //!< The actual IJ-pair container

  std::map<size_t,size_t> xnMap; //!< Node index map used by getCoord()
  std::map<size_t,size_t> nxMap; //!< Node index map used by getNodeID()

  //! Inhomogeneous Dirichlet boundary condition data
  std::vector<DirichletEdge> dirich;

  //! Element groups for multi-threaded assembly
  ThreadGroups threadGroups;

  //! Basis function cache
  std::vector<std::unique_ptr<BasisFunctionCache>> myCache;
};

#endif
