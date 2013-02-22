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

#include "ASM3D.h"
#include "ASMunstruct.h"

class FiniteElement;

namespace Go {
	class SplineSurface;
	class SplineVolume;
}

namespace LR {
	class LRSplineSurface;
	class LRSplineVolume;
}


/*!
	\brief Driver for assembly of unstructured 3D spline FE models.
	\details This class contains methods common for unstructured 3D spline patches.
*/

class ASMu3D : public ASMunstruct, public ASM3D
{

public:
	//! \brief Default constructor.
	ASMu3D(unsigned char n_f = 3);
	//! \brief Copy constructor.
	ASMu3D(const ASMu3D& patch, unsigned char n_f = 0);
	//! \brief Empty destructor.
	virtual ~ASMu3D() {}

	//! \brief Returns the spline volume representing the geometry of this patch.
	LR::LRSplineVolume* getVolume() const { return lrspline; }
	//! \brief Returns the spline surface representing a boundary of this patch.
	//! \param[in] dir Parameter direction defining which boundary to return
	virtual LR::LRSplineSurface* getBoundary(int dir);
	//! \brief Returns the spline volume representing the basis of this patch.
	virtual LR::LRSplineVolume* getBasis(int = 1) const { return lrspline; }


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

	//! \brief Returns the total number of nodes in this patch.
	virtual size_t getNoNodes(int basis = 0) const;

	//! \brief Returns a matrix with nodal coordinates for an element.
	//! \param[in] iel Element index
	//! \param[out] X 3\f$\times\f$n-matrix,
	//! where \a n is the number of nodes in one element
	virtual bool getElementCoordinates(Matrix& X, int iel) const;

	//! \brief Returns a matrix with all nodal coordinates within the patch.
	//! \param[out] X 3\f$\times\f$n-matrix,
	//! where \a n is the number of nodes in the patch
	virtual void getNodalCoordinates(Matrix& X) const;

	//! \brief Returns the global coordinates for the given node.
	//! \param[in] inod 1-based node index local to current patch
	virtual Vec3 getCoord(size_t inod) const;


	//! \brief Updates the nodal coordinates for this patch.
	//! \param[in] displ Incremental displacements to update the coordinates with
	virtual bool updateCoords(const Vector& displ);

	//! \brief Checks that the patch is modelled in a right-hand-side system.
	//! \details If it isn't, the w-parameter direction is swapped.
	bool checkRightHandSystem();

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

	bool refine(const std::vector<int>& elements,
	            const std::vector<int>& options, const char* fName = 0) {
		std::cout << "ASMu3D::refine()" << std::endl;
		return ASMunstruct::refine(elements, options, fName);
	}


	// Various methods for preprocessing of boundary conditions and patch topology
	// ===========================================================================

	//! \brief Constrains all DOFs on a given boundary face.
	//! \param[in] dir Parameter direction defining the face to constrain
	//! \param[in] open If \e true, exclude all points along the face boundary
	//! \param[in] dof Which DOFs to constrain at each node on the face
	//! \param[in] code Inhomogeneous dirichlet condition code
	void constrainFace(int dir, bool open, int dof = 123, int code = 0);
	//! \brief Constrains all DOFs in local directions on a given boundary face.
	//! \param[in] dir Parameter direction defining the face to constrain
	//! \param[in] open If \e true, exclude all points along the face boundary
	//! \param[in] dof Which local DOFs to constrain at each node on the face
	//! \param[in] code Inhomogeneous dirichlet condition code
	//! \param[in] project If \e true, the local axis directions are projected
	//! \param[in] T1 Desired global direction of first local tangent direction
	//! \return Number of additional nodes added due to local axis constraints
	size_t constrainFaceLocal(int dir, bool open, int dof = 3, int code = 0,
	                          bool project = false, char T1 = '\0');

	//! \brief Constrains all DOFs on a given boundary edge.
	//! \param[in] lEdge Local index [1,12] of the edge to constrain
	//! \param[in] open If \e true, exclude the end points of the edge
	//! \param[in] dof Which DOFs to constrain at each node on the edge
	//! \param[in] code Inhomogeneous dirichlet condition code
	void constrainEdge(int lEdge, bool open, int dof = 123, int code = 0);

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
	void constrainLine(int fdir, int ldir, double xi,
	                   int dof = 123, int code = 0);

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
	void constrainCorner(int I, int J, int K,
	                     int dof = 123, int code = 0);
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
	void constrainNode(double xi, double eta, double zeta,
	                   int dof = 123, int code = 0);

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
	virtual bool connectPatch(int face, ASMu3D& neighbor, int nface,
	                          int norient = 0);

	//! \brief Makes two opposite boundary faces periodic.
	//! \param[in] dir Parameter direction defining the periodic faces
	//! \param[in] basis Which basis to connect (mixed methods), 0 means both
	//! \param[in] master 1-based index of the first master node in this basis
	virtual void closeFaces(int dir, int basis = 0, int master = 1);

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

public:
	//! \brief Projects the secondary solution field onto the primary basis.
	//! \param[in] integrand Object with problem-specific data and methods
	virtual Go::GeomObject* evalSolution(const IntegrandBase& integrand) const { return NULL; }

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

	//! \brief Projects the secondary solution using a discrete global L2-norm.
	//! \param[out] sField Secondary solution field control point values
	//! \param[in] integrand Object with problem-specific data and methods
	//! \param[in] continuous If \e true, a continuous L2-projection is used
	virtual bool globalL2projection(Matrix& sField,
	                                const IntegrandBase& integrand,
	                                bool continuous = false) const { return false; }

	//! \brief Evaluate and interpolate a field over a given geometry
	//! \param[in] input The basis of the field to evaluate
	//! \param[in] locVec The coefficients of the field
	//! \param[out] vec The obtained coefficients after interpolation
	virtual bool evaluate(const ASMbase* input,
				                const Vector& locVec, Vector& vec);

protected:

	// Internal utility methods
	// ========================

	//! \brief Connects all matching nodes on two adjacent boundary faces.
	//! \param[in] face Local face index of this patch, in range [1,6]
	//! \param neighbor The neighbor patch
	//! \param[in] nface Local face index of neighbor patch, in range [1,6]
	//! \param[in] norient Relative face orientation flag (see \a connectPatch)
	//! \param[in] basis Which basis to connect the nodes for (mixed methods)
	//! \param[in] slave 0-based index of the first slave node in this basis
	//! \param[in] master 0-based index of the first master node in this basis
	bool connectBasis(int face, ASMu3D& neighbor, int nface, int norient,
	                  int basis = 1, int slave = 0, int master = 0);

	//! \brief Extracts parameter values of the Gauss points in one direction.
	//! \param[out] uGP Parameter values in given direction for all points
	//! \param[in] dir Parameter direction (0,1,2)
	//! \param[in] nGauss Number of Gauss points along a knot-span
	//! \param[in] iel Element index
	//! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
	//! \return The parameter value matrix casted into a one-dimensional vector
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
	//! \param[in] iel Element index
	double getParametricVolume(int iel) const;

	//! \brief Returns boundary face area in the parameter space for an element.
	//! \param[in] iel Element index
	//! \param[in] dir Local face index of the boundary face
	double getParametricArea(int iel, int dir) const;

	//! \brief Computes the element corner coordinates.
	//! \param[in] iel index (0 based) of the element on which corners are requested
	//! \param[out] XC Coordinates of the element corners
	void getElementCorners(int iel, std::vector<Vec3>& XC) const;

	//! \brief Evaluate all basis functions and \a derivs number of derivatives on one element
	virtual void evaluateBasis(FiniteElement &el, int derivs) const;

	//! \brief Evaluate all basis functions and first derivatives on one element
	virtual void evaluateBasis(FiniteElement &el, Matrix &dNdu, Matrix &C, Matrix &B) const ;

	//! \brief Evaluate all basis functions and first derivatives on one element
	virtual void evaluateBasis(FiniteElement &el, Matrix &dNdu) const;

	//! \brief Evaluate all basis functions and second order derivatives on one element
	virtual void evaluateBasis(FiniteElement &el, Matrix &dNdu, Matrix3D& d2Ndu2) const;

public:
	//! \brief Returns the polynomial order in each parameter direction.
	//! \param[out] p1 Order in first (u) direction
	//! \param[out] p2 Order in second (v) direction
	//! \param[out] p3 Order in third (w) direction
	virtual bool getOrder(int& p1, int& p2, int& p3) const;

#if 0
	//! \brief Returns the number of elements on a boundary.
	virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

	//! \brief Returns the number of nodal points in each parameter direction.
	//! \param[out] n1 Number of nodes in first (u) direction
	//! \param[out] n2 Number of nodes in second (v) direction
	//! \param[out] n3 Number of nodes in third (w) direction
	//! \param[in] basis Which basis to return size parameters for (mixed methods)
	virtual bool getSize(int& n1, int& n2, int& n3, int basis = 0) const;
#endif

	// int getNodeID(size_t inod, bool = false)          const { return inod; };
	// unsigned char getNodalDOFs(size_t inod)           const { return 0; };
	// void getNoIntPoints(size_t& nPt)                        {  };
	// void getNoBouPoints(size_t& nPt, char ldim, char lindx) {  };
	// bool getSolution(Matrix& sField, const Vector& locSol,
	//                  const IntVec& nodes)             const { return false; };
	// void extractNodeVec(const Vector& globVec, Vector& nodeVec,
	//                     unsigned char nndof = 0, int basis = 1) const { };
	// bool injectNodeVec(const Vector& nodeVec, Vector& globVec,
	//                    unsigned char nndof = 0) const { return false;};

private:
	//! \brief Returns an index into the internal coefficient array for a node.
	//! \param[in] inod 0-based node index local to current patch
	int coeffInd(size_t inod) const;

protected:
	LR::LRSplineVolume* lrspline;  //!< Pointer to the actual spline volume object
	bool                 swapW; //!< Has the w-parameter direction been swapped?

	std::map<size_t,size_t> xnMap; //!< Node index map used by \a getCoord
	std::map<size_t,size_t> nxMap; //!< Node index map used by \a getNodeID
	std::vector<Matrix> bezierExtract; //!< The Bezier exctraction matrices for all elements


private:
	Go::SplineVolume* tensorspline; //!< Pointer to original tensor spline object
	// The tensor spline object is kept for backward compatability with the REFINE
	// and RAISEORDER key-words, although we take note that there is a possibility
	// of optimization since all mapping values and Jacobians may be performed on
	// this object for increased efficiency.

	mutable int workingEl; //!< here to keep track of element across function calls (to avoid topological searches all the time)
};

#endif
