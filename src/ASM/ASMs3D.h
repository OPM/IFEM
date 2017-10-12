// $Id$
//==============================================================================
//!
//! \file ASMs3D.h
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D spline FE models.
//!
//==============================================================================

#ifndef _ASM_S3D_H
#define _ASM_S3D_H

#include "ASMstruct.h"
#include "ASM3D.h"
#include "ThreadGroups.h"

namespace Go {
  class SplineSurface;
  class SplineVolume;
}


/*!
  \brief Driver for assembly of structured 3D spline FE models.
  \details This class contains methods common for structured 3D spline patches.
*/

class ASMs3D : public ASMstruct, public ASM3D
{
  //! \brief Struct for nodal point data.
  struct IJK
  {
    int I; //!< Index in first parameter direction
    int J; //!< Index in second parameter direction
    int K; //!< Index in third parameter direction
  };

  typedef std::vector<IJK> IndexVec; //!< Node index container

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

  //! \brief Struct for face node definitions.
  struct Face
  {
    int isnod; //!< Global node number of the first interior point on the face
    int incrI; //!< Increment in global numbering in the I-direction on the face
    int incrJ; //!< Increment in global numbering in the J-direction on the face
    int nnodI; //!< Number of nodes in the local I-direction on the face
    int indxI; //!< Running node index in the local I-direction

    //! \brief Default constructor.
    Face() { isnod = incrI = incrJ = nnodI = 0; indxI = 1; }
    //! \brief Returns \a isnod which then is incremented.
    int next();
  };

public:
  //! \brief Struct with data for definition of global node numbers of a patch.
  struct BlockNodes
  {
    int  ibnod[8];  //!< Vertex nodes
    Edge edges[12]; //!< Edge nodes
    Face faces[6];  //!< Face nodes
    int  iinod;     //!< Global node number of the first interior node
    int  inc[3];    //!< Increment in global node numbering in each direction
    int  nnodI;     //!< Number of nodes in parameter direction I
    int  nnodJ;     //!< Number of nodes in parameter direction J
    int  indxI;     //!< Running node index in the local I-direction
    int  indxJ;     //!< Running node index in the local J-direction

    //! \brief Default constructor.
    BlockNodes() { iinod = inc[0] = inc[1] = inc[2] = 0; indxI = indxJ = 1; }
    //! \brief Returns \a iinod which then is incremented.
    int next();
  };

private:
  typedef std::pair<int,int> Ipair; //!< Convenience type

  //! \brief Struct representing an inhomogeneous Dirichlet boundary condition.
  struct DirichletFace
  {
    Go::SplineSurface* surf;  //!< Pointer to spline surface for the boundary
    int                dof;   //!< Local DOF to constrain along the boundary
    int                code;  //!< Inhomogeneous Dirichlet condition code
    std::vector<Ipair> nodes; //!< Nodes subjected to projection on the boundary
    //! \brief Default constructor.
    DirichletFace(Go::SplineSurface* ss = nullptr, int d = 0, int c = 0)
    : surf(ss), dof(d), code(c) {}
  };

protected:
  //! \brief Base class that checks if an element has interface contributions.
  class InterfaceChecker
  {
    const ASMs3D& myPatch; //!< Reference to the patch being integrated
  public:
    //! \brief The constructor initialises the reference to current patch.
    explicit InterfaceChecker(const ASMs3D& pch) : myPatch(pch) {}
    //! \brief Empty destructor.
    virtual ~InterfaceChecker() {}
    //! \brief Returns non-zero if the specified element have contributions.
    //! \param[in] I Index in first parameter direction of the element
    //! \param[in] J Index in second parameter direction of the element
    //! \param[in] K Index in third parameter direction of the element
    virtual short int hasContribution(int I, int J, int K) const;
  };

public:
  //! \brief Default constructor.
  explicit ASMs3D(unsigned char n_f = 3);
  //! \brief Special copy constructor for sharing of FE data.
  ASMs3D(const ASMs3D& patch, unsigned char n_f);
  //! \brief Default copy constructor copying everything.
  ASMs3D(const ASMs3D& patch);
  //! \brief Empty destructor.
  virtual ~ASMs3D() {}

  //! \brief Returns the spline volume representing the geometry of this patch.
  Go::SplineVolume* getVolume() const { return svol; }
  //! \brief Returns the spline surface representing a boundary of this patch.
  //! \param[in] dir Parameter direction defining which boundary to return
  virtual Go::SplineSurface* getBoundary(int dir, int = 1);
  //! \brief Returns the spline volume representing the basis of this patch.
  virtual Go::SplineVolume* getBasis(int = 1) const { return svol; }
  //! \brief Copies the parameter domain from the \a other patch.
  virtual void copyParameterDomain(const ASMbase* other);

  // Methods for model generation
  // ============================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineVolume object to given stream.
  virtual bool write(std::ostream&, int = 0) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the node-to-IJK-index array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Adds extraordinary elements associated with a patch boundary.
  //! \param[in] dim Dimension of the boundary (should be 2)
  //! \param[in] item Local index of the boundary face
  //! \param[in] nXn Number of extraordinary nodes
  //! \param[out] nodes Global numbers assigned to the extraordinary nodes
  virtual bool addXElms(short int dim, short int item,
                        size_t nXn, IntVec& nodes);

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

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face
  //! \param nodes Array of node numbers
  //! \param basis Which basis to grab nodes for (0 for all)
  //! \param thick Thickness of connection
  //! \param local If true return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int basis, int thick = 1,
                                int = 0, bool local = false) const;

  //! \brief Returns the node index for a given corner.
  //! \param I     -1 or +1 for either umin or umax corner
  //! \param J     -1 or +1 for either vmin or vmax corner
  //! \param K     -1 or +1 for either wmin or wmax corner
  //! \param basis which basis to consider (only applicable for mixed methods)
  virtual int getCorner(int I, int J, int K, int basis) const;

  //! \brief Returns the (1-indexed) node indices for a given edge.
  //! \param lEdge  index to local edge (1,2,...12)
  //! \param open   include end points or not
  //! \param basis  which basis to consider (only applicable for mixed methods)
  //! \param orient local orientation of indices, see ASMunstruct::Sort
  virtual std::vector<int> getEdge(int lEdge, bool open, int basis, int = 0) const;

  //! \brief Assigns new global node numbers for all nodes of the patch.
  //! \param nodes Object with global nodes numbers to assign to this patch
  //! \param[in] basis Which basis to assign node numbers for in mixed methods
  //!
  //! \details The global node numbers generated by \a generateFEMTopology are
  //! non-unique in the sense that a node that is shared by two (or more)
  //! patches along a common interface has a different number in each patch.
  //! This method therefore assigns a new global number to each node in the
  //! patch. The data provided through the \a nodes argument is sufficient
  //! to determine the unique global number under the assumption that they
  //! are ordered in the sequence determined by the local orientation of the
  //! patch, its faces and edges.
  bool assignNodeNumbers(BlockNodes& nodes, int basis = 0);

  //! \brief Checks that the patch is modelled in a right-hand-side system.
  //! \details If it isn't, the w-parameter direction is swapped.
  virtual bool checkRightHandSystem();

  //! \brief Refines the parametrization by inserting extra knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  bool refine(int dir, const RealArray& xi);
  //! \brief Refines the parametrization by inserting extra knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  bool uniformRefine(int dir, int nInsert);
  //! \brief Raises the order of the SplineVolume object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  //! \param[in] rw Number of times to raise the order in w-direction
  bool raiseOrder(int ru, int rv, int rw);


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains all DOFs on a given boundary face.
  //! \param[in] dir Parameter direction defining the face to constrain
  //! \param[in] open If \e true, exclude all points along the face boundary
  //! \param[in] dof Which DOFs to constrain at each node on the face
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Basis to constrain edge for
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
  //! \param[in] basis Basis to constrain edge for
  virtual void constrainEdge(int lEdge, bool open, int dof,
                             int code = 0, char basis = 1);

  //! \brief Constrains all DOFs along a line on a given boundary face.
  //! \param[in] fdir Parameter direction defining the face to constrain
  //! \param[in] ldir Parameter direction defining the line to constrain
  //! \param[in] xi Parameter value defining the line to constrain
  //! \param[in] dof Which DOFs to constrain at each node along the line
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Basis to constrain line for
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
  //! \param[in] basis Basis to constrain corner for
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
  //! \param[in] basis Basis to constrain node for
  //!
  //! \details The parameter values have to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual void constrainNode(double xi, double eta, double zeta, int dof,
                             int code = 0, char basis = 1);

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag (see below)
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //! \param[in] thick Thickness of connection
  //!
  //! \details The face orientation flag \a norient must be in range [0,7].
  //! When interpreted as a binary number, its 3 digits are decoded as follows:
  //! - left digit = 1: The u and v parameters of the neighbor face are swapped
  //! - middle digit = 1: Parameter \a u in neighbor patch face is reversed
  //! - right digit = 1: Parameter \a v in neighbor patch face is reversed
  virtual bool connectPatch(int face, ASM3D& neighbor, int nface, int norient,
                            int = 0, bool coordCheck = true, int thick = 1);

  //! \brief Makes two opposite boundary faces periodic.
  //! \param[in] dir Parameter direction defining the periodic faces
  //! \param[in] basis Which basis to connect (mixed methods), 0 means both
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeFaces(int dir, int basis = 0, int master = 1);

  //! \brief Sets the global node numbers for this patch.
  //! \param[in] nodes Vector of global node numbers (zero-based)
  virtual void setNodeNumbers(const std::vector<int>& nodes);

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

protected:
  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] itgPts Parameters and weights of the integration points
  bool integrate(Integrand& integrand, GlobalIntegral& glbInt,
                 const TimeDomain& time, const Real3DMat& itgPts);

  //! \brief Evaluates an integral over element interfaces in the patch.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] iChk Object checking if an element interface has contributions
  bool integrate(Integrand& integrand, GlobalIntegral& glbInt,
                 const TimeDomain& time, const InterfaceChecker& iChk);

public:

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
  //! \param[in] regular Flag indicating how the sampling points are defined
  //! \param[in] deriv Derivative order to return
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size() \a X \a gpar[2].size().
  //! Otherwise, we assume that it contains the \a u, \a v and \a w parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool regular = true,
                            int deriv = 0, int = 0) const;

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
  Go::SplineVolume* projectSolution(const IntegrandBase& integrand) const;

  //! \brief Projects the secondary solution using Quasi-Interpolation.
  Go::SplineVolume* projectSolutionLocal(const IntegrandBase&) const;

  //! \brief Projects the secondary solution using
  //! Variation Diminishing Spline Approximation.
  Go::SplineVolume* projectSolutionLocalApprox(const IntegrandBase&) const;

  //! \brief Projects the secondary solution using Least-Square Approximation.
  Go::SplineVolume* projectSolutionLeastSquare(const IntegrandBase&) const;

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
  //! \a gpar[0].size() \a X \a gpar[1].size() \a X \a gpar[2].size().
  //! Otherwise, we assume that it contains the \a u, \a v and \a w parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
			    const RealArray* gpar, bool regular = true) const;

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
  //! \param[in] coordCheck False to turn off coordinate checks
  //! \param[in] thick Thickness of connection
  bool connectBasis(int face, ASMs3D& neighbor, int nface, int norient,
                    int basis = 1, int slave = 0, int master = 0,
                    bool coordCheck = true, int thick = 1);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  //! \return The parameter value matrix casted into a one-dimensional vector
  const Vector& getGaussPointParameters(Matrix& uGP, int dir, int nGauss,
					const double* xi) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] basis Basis number (mixed)
  bool getGrevilleParameters(RealArray& prm, int dir, int basis=1) const;

  //! \brief Calculates parameter values for the Quasi-Interpolation points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  bool getQuasiInterplParameters(RealArray& prm, int dir) const;

  //! \brief Returns the volume in the parameter space for an element.
  //! \param[in] iel Element index
  double getParametricVolume(int iel) const;
  //! \brief Returns boundary face area in the parameter space for an element.
  //! \param[in] iel Element index
  //! \param[in] dir Local face index of the boundary face
  double getParametricArea(int iel, int dir) const;

  //! \brief Computes the element border parameters.
  //! \param[in] i1 Parameter index in u-direction
  //! \param[in] i2 Parameter index in v-direction
  //! \param[in] i3 Parameter index in w-direction
  //! \param[out] u Parameter values of the west-east borders
  //! \param[out] v Parameter values of the south-north borders
  //! \param[out] w Parameter values of the back-front borders
  void getElementBorders(int i1, int i2, int i3,
                         double* u, double* v, double* w) const;

  //! \brief Computes the element corner coordinates.
  //! \param[in] i1 Parameter index in u-direction
  //! \param[in] i2 Parameter index in v-direction
  //! \param[in] i3 Parameter index in w-direction
  //! \param[out] XC Coordinates of the element corners
  //! \return Characteristic element size
  double getElementCorners(int i1, int i2, int i3, std::vector<Vec3>& XC) const;

  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM If \e true, ignore global multipliers in sanity check
  virtual void generateThreadGroups(const Integrand& integrand, bool silence,
                                    bool ignoreGlobalLM);
  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] strip1 Strip width in first direction
  //! \param[in] strip2 Strip width in second direction
  //! \param[in] strip3 Strip width in third direction
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM If \e true, ignore global multipliers in sanity check
  void generateThreadGroups(size_t strip1, size_t strip2, size_t strip3,
                            bool silence, bool ignoreGlobalLM);
  //! \brief Generates element groups for multi-threading of boundary integrals.
  //! \param[in] lIndex Local index [1,6] of the boundary face
  //! \param[in] silence If \e true, suppress threading group outprint
  virtual void generateThreadGroups(char lIndex, bool silence, bool);

public:
  //! \brief Auxilliary function for computation of basis function indices.
  static void scatterInd(int n1, int n2, int n3, int p1, int p2, int p3,
			 const int* start, IntVec& index);

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction
  virtual bool getOrder(int& p1, int& p2, int& p3) const;

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  //! \param[in] basis Which basis to return size parameters for (mixed methods)
  virtual bool getSize(int& n1, int& n2, int& n3, int basis = 0) const;

  //! \brief Returns the number of elements in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  virtual bool getNoStructElms(int& n1, int& n2, int& n3) const;

  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

  //! \brief Establishes matrices with basis functions and 1st derivatives.
  //! \param[in] u First parameter value of current integration point
  //! \param[in] v Second parameter value of current integration point
  //! \param[in] w Third parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  //! \param[in] fromRight If \e true, evaluate from right if at a knot
  void extractBasis(double u, double v, double w, Vector& N, Matrix& dNdu,
                    bool fromRight = true) const;
  //! \brief Establishes matrices with basis functions, 1st and 2nd derivatives.
  //! \param[in] u First parameter value of current integration point
  //! \param[in] v Second parameter value of current integration point
  //! \param[in] w Third parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  //! \param[out] d2Ndu2 Second derivatives of basis functions
  //! \param[in] fromRight If \e true, evaluate from right if at a knot
  void extractBasis(double u, double v, double w, Vector& N,
                    Matrix& dNdu, Matrix3D& d2Ndu2,
                    bool fromRight = true) const;
  //! \brief Establishes a vector with basis function derivatives.
  //! \param[in] u First parameter value of current integration point
  //! \param[in] v Second parameter value of current integration point
  //! \param[in] w Third parameter value of current integration point
  //! \param[in] dir Which parameter to establish derivatives with respect to
  //! \param[in] p The derivation order
  //! \param[out] dN Basis function derivatives
  //! \param[in] fromRight If \e true, evaluate from right if at a knot
  void extractBasis(double u, double v, double w, int dir, int p, Vector& dN,
                    bool fromRight = true) const;

private:
  //! \brief Returns an index into the internal coefficient array for a node.
  //! \param[in] inod 0-based node index local to current patch
  int coeffInd(size_t inod) const;

  //! \brief Find the start node and size for a basis.
  //! \param[out] n1 Size of basis in first parameter direction
  //! \param[out] n2 Size of basis in second parameter direction
  //! \param[out] n3 Size of basis in third parameter direction
  //! \param[in] basis Basis to find data for
  //! \return 1-based index of start node for basis
  int findStartNode(int& n1, int& n2, int& n3, char basis) const;

  //! \brief Find local sizes for a given face.
  //! \param[out] n1 Number of nodes in first local parameter direction on face
  //! \param[out] n2 Number of nodes in second local parameter direction face
  //! \param[in] basis Basis to obtain sizes for
  //! \param[in] face Face to obtain sizes for
  bool getFaceSize(int& n1, int& n2, int basis, int face) const;

protected:
  Go::SplineVolume* svol;  //!< Pointer to the actual spline volume object
  bool              swapW; //!< Has the w-parameter direction been swapped?

  const IndexVec& nodeInd; //!< IJK-triplets for the control points (nodes)
  IndexVec      myNodeInd; //!< The actual IJK-triplet container

  std::map<size_t,size_t> xnMap; //!< Node index map used by \a getCoord
  std::map<size_t,size_t> nxMap; //!< Node index map used by \a getNodeID

  //! Inhomogeneous Dirichlet boundary condition data
  std::vector<DirichletFace> dirich;

  //! Element groups for multi-threaded volume assembly
  ThreadGroups                threadGroupsVol;
  //! Element groups for multi-threaded face assembly
  std::map<char,ThreadGroups> threadGroupsFace;
};

#endif
