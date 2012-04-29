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
#include "Utilities.h"

namespace Go {
  class SplineCurve;
  class SplineSurface;
  class BasisDerivsSf;
  class BasisDerivsSf2;
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
    BlockNodes() { iinod = inc[0] = inc[1] = 0; indxI = 1; }
    //! \brief Returns \a iinod which then is incremented.
    int next();
  };

private:
  typedef std::pair<int,int> Ipair; //!< Convenience declaration

  //! \brief Struct representing an inhomogeneous Dirichlet boundary condition.
  struct DirichletEdge
  {
    Go::SplineCurve*   curve; //!< Pointer to spline curve for the boundary
    int                dof;   //!< Local DOF to constrain along the boundary
    int                code;  //!< Inhomogeneous Dirichlet condition code
    std::vector<Ipair> nodes; //!< Nodes subjected to projection on the boundary
    //! \brief Default constructor.
    DirichletEdge(Go::SplineCurve* sc = 0, int d = 0, int c = 0)
    : curve(sc), dof(d), code(c) {}
  };

public:
  //! \brief Default constructor.
  ASMs2D(unsigned char n_s = 2, unsigned char n_f = 2);
  //! \brief Copy constructor.
  ASMs2D(const ASMs2D& patch, unsigned char n_f = 0);
  //! \brief The destructor frees the dynamically allocated boundary curves.
  virtual ~ASMs2D();

  //! \brief Returns the spline surface representing the geometry of this patch.
  Go::SplineSurface* getSurface() const { return surf; }
  //! \brief Returns the spline curve representing a boundary of this patch.
  //! \param[in] dir Parameter direction defining which boundary to return
  Go::SplineCurve* getBoundary(int dir);
  //! \brief Returns the spline surface representing the basis of this patch.
  virtual Go::SplineSurface* getBasis(int = 1) const { return surf; }


  // Methods for model generation
  // ============================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineSurface object to given stream.
  virtual bool write(std::ostream&, int = 0) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the node-to-IJ-index array, as well as global node and element numbers.
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
  //! patch and its edges.
  bool assignNodeNumbers(BlockNodes& nodes, int basis = 0);

  //! \brief Refines the parametrization by inserting extra knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  bool refine(int dir, const RealArray& xi);
  //! \brief Refines the parametrization by inserting extra knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  bool uniformRefine(int dir, int nInsert);
  //! \brief Raises the order of the SplineSurface object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  bool raiseOrder(int ru, int rv);


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] dir Parameter direction defining the edge to constrain
  //! \param[in] dof Which DOFs to constrain at each node on the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  virtual void constrainEdge(int dir, int dof, int code = 0);
  //! \brief Constrains all DOFs in local directions on a given boundary edge.
  //! \param[in] dir Parameter direction defining the edge to constrain
  //! \param[in] dof Which local DOFs to constrain at each node on the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] project If \e true, the local axis directions are projected
  //! \return Number of additional nodes added due to local axis constraints
  virtual size_t constrainEdgeLocal(int dir, int dof, int code, bool project);

  //! \brief Constrains a corner node identified by the two parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The sign of the two indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  virtual void constrainCorner(int I, int J, int dof, int code = 0);
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
  virtual bool connectPatch(int edge, ASMs2D& neighbor, int nedge,
			    bool revers = false);

  //! \brief Makes two opposite boundary edges periodic.
  //! \param[in] dir Parameter direction defining the periodic edges
  //! \param[in] basis Which basis to connect (mixed methods), 0 means both
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeEdges(int dir, int basis = 0, int master = 1);

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  //! \param[in] time Current time
  virtual bool updateDirichlet(const std::map<int,RealFunc*>& func,
			       const std::map<int,VecFunc*>& vfunc,
			       double time = 0.0);


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
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size().
  //! Otherwise, we assume that it contains the \a u and \a v parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
			    const RealArray* gpar, bool regular = true) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] project Flag indicating result recovery method
  //! (0=none, 'D'=direct evaluation, 'S'=superconvergent recovery)
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

private:
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  Go::SplineSurface* projectSolution(const IntegrandBase& integrand) const;

  //! \brief Projects the secondary solution using a superconvergent approach.
  Go::SplineSurface* scRecovery(const IntegrandBase&) const;

  //! \brief Projects the secondary solution field using Quasi-Interpolation.
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

  //! \brief Projects the secondary solution using a discrete global L2-norm.
  //! \param[out] sField Secondary solution field control point values
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e true, a continuous L2-projection is used
  virtual bool globalL2projection(Matrix& sField,
				  const IntegrandBase& integrand,
				  bool continuous = false) const;

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
  bool connectBasis(int edge, ASMs2D& neighbor, int nedge, bool revers,
		    int basis = 1, int slave = 0, int master = 0);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  //! \return The parameter value matrix casted into a one-dimensional vector
  const Vector& getGaussPointParameters(Matrix& uGP, int dir, int nGauss,
					const double* xi) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  bool getGrevilleParameters(RealArray& prm, int dir) const;

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

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[in] basis Which basis to return size parameters for (mixed methods)
  virtual bool getSize(int& n1, int& n2, int basis = 0) const;

  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] silence If \e true, suppress threading group outprint
  virtual void generateThreadGroups(bool silence);

  //! \brief Establishes matrices with basis functions and 1st derivatives.
  static void extractBasis(const Go::BasisDerivsSf& spline,
			   Vector& N, Matrix& dNdu);
  //! \brief Establishes matrices with basis functions, 1st and 2nd derivatives.
  static void extractBasis(const Go::BasisDerivsSf2& spline,
			   Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2);

public:
  //! \brief Auxilliary function for computation of basis function indices.
  static void scatterInd(int n1, int n2, int p1, int p2,
			 const int* start, IntVec& index);

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  bool getOrder(int& p1, int& p2) const;
  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction
  //! \param[out] p3 Order in third (w) direction (always zero)
  virtual bool getOrder(int& p1, int& p2, int& p3) const
  { p3 = 0; return this->getOrder(p1,p2); }

  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

private:
  //! \brief Returns an index into the internal coefficient array for a node.
  //! \param[in] inod 0-based node index local to current patch
  int coeffInd(size_t inod) const;

protected:
  Go::SplineSurface* surf; //!< Pointer to the actual spline surface object
  Go::SplineCurve* bou[4]; //!< Pointers to the four boundary curves

  const IndexVec& nodeInd; //!< IJ-pairs for the control points (nodes)
  IndexVec      myNodeInd; //!< The actual IJ-pair container

  //! Inhomogeneous Dirichlet boundary condition data
  std::vector<DirichletEdge> dirich;

  //! Element groups for multi-threaded assembly
  utl::ThreadGroups threadGroups;
};

#endif
