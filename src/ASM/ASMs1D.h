// $Id$
//==============================================================================
//!
//! \file ASMs1D.h
//!
//! \date Apr 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 1D spline FE models.
//!
//==============================================================================

#ifndef _ASM_S1D_H
#define _ASM_S1D_H

#include "ASMstruct.h"
#include "ASM1D.h"
#include "Tensor.h"

typedef std::vector<Tensor> TensorVec; //!< An array of non-symmetric tensors

namespace Go {
  class SplineCurve;
}


/*!
  \brief Driver for assembly of structured 1D spline FE models.
  \details This class contains methods common for structured 1D spline patches.
*/

class ASMs1D : public ASMstruct, public ASM1D
{
public:
  //! \brief Default constructor.
  ASMs1D(unsigned char n_s = 1, unsigned char n_f = 1);
  //! \brief Copy constructor.
  ASMs1D(const ASMs1D& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMs1D() {}

  //! \brief Returns the spline curve representing the geometry of this patch.
  Go::SplineCurve* getCurve() const { return curv; }


  // Methods for model generation
  // ============================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);
  //! \brief Writes the geometry of the SplineCurve object to given stream.
  virtual bool write(std::ostream&, int = 0) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! and the global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Generates a beam finite element model for the patch.
  //! \param[in] Zaxis Vector defining a point in the local XZ-plane
  virtual bool generateOrientedFEModel(const Vec3& Zaxis);

  //! \brief Generates a twisted beam finite element model for the patch.
  //! \param[in] twist Function describing the twist angle along the beam
  //! \param[in] Zaxis Vector defining a point in the local XZ-plane
  bool generateTwistedFEModel(const RealFunc& twist, const Vec3& Zaxis);

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Returns the number of nodal points in the patch.
  //! \param[in] basis Which basis to return size parameters for (mixed methods)
  virtual int getSize(int basis = 0) const;

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Returns the current rotation tensor the given node.
  //! \param[in] inod 1-based node index local to current patch
  Tensor getRotation(size_t inod) const;

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel Element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const;

  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const;

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ);

  //! \brief Updates the nodal rotations for this patch.
  //! \param[in] displ Incremental displacements to update the rotations with
  //! \param[in] reInit If \e true, reinitialize rotations from unity
  bool updateRotations(const Vector& displ, bool reInit = false);
  //! \brief Updates the previous nodal rotations for this patch at convergence.
  void updateRotations() { prevT = myT; }

  //! \brief Finds the global (or patch-local) node number on a patch end.
  //! \param[in] lIndex Local index of the end point
  //! \param nodes Array of global boundary node numbers
  //! \param[in] thick Thickness of connection
  //! \param[in] local If \e true, return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int, int thick, int, bool local) const;

  //! \brief Finds the node that is closest to the given point.
  //! \param[in] X Global coordinates of point to search for
  //! \return 1-based nodal index and distance to to the found node
  virtual std::pair<size_t,double> findClosestNode(const Vec3& X) const;

  //! \brief Refines the parametrization by inserting extra knots.
  //! \param[in] xi Relative positions of added knots in each existing knot span
  bool refine(const RealArray& xi);
  //! \brief Refines the parametrization by inserting extra knots uniformly.
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  bool uniformRefine(int nInsert);
  //! \brief Raises the order of the SplineCurve object for this patch.
  //! \param[in] ru Number of times to raise the order
  bool raiseOrder(int ru);


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains a node identified by a relative parameter value.
  //! \param[in] xi Parameter value along the curve
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis The basis to constrain node for
  //! \return 1-based index of the constrained node
  //!
  //! \details The parameter value has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual int constrainNode(double xi, int dof, int code = 0, char basis = 1);

  //! \brief Connects matching nodes on two adjacent vertices.
  //! \param[in] vertex Local vertex index of this patch, in range [1,2]
  //! \param neighbor The neighbor patch
  //! \param[in] nvertex Local vertex index of neighbor patch, in range [1,2]
  //! \param[in] thick Thickness of connection
  virtual bool connectPatch(int vertex, ASMs1D& neighbor, int nvertex,
                            int thick = 1);

  //! \brief Makes the two end vertices of the curve periodic.
  //! \param[in] basis Which basis to connect (mixed methods), 0 means both
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeEnds(int basis = 0, int master = 1);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
			 GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates a boundary integral over a patch end.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the end point
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time);


  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameter in range [0.0,1.0] of the point
  //! \param[out] param The parameter of the point in knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that matches the point, if any
  //! \return 0 if no node (control point) matches this point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const;

  //! \brief Returns the element that contains a specified spatial point.
  //! \param[in] param The parameter of the point in the knot-span domain
  //! \return Local element number within the patch that contains the point
  virtual int findElementContaining(const double* param) const;

  //! \brief Creates a line element model of this patch for visualization.
  //! \param[out] grid The generated line grid
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
  //! \param[in] locSol Solution vector in DOF-order
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int* npe, int = 0) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] deriv Derivative order to return
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool = true,
                            int deriv = 0, int = 0) const;

  using ASMbase::evaluate;
  //! \brief Evaluates and interpolates a function over a given geometry.
  //! \param[in] func The function to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basis Basis number (mixed)
  //! \param[in] time Current time
  virtual bool evaluate(const FunctionBase* func, RealArray& vec,
                        int basis, double time) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] project Flag indicating the projection method
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! If \a npe is null, the solution is evaluated at the Greville points and
  //! then projected onto the spline basis to obtain the control point values,
  //! which then are returned through \a sField.
  //! If \a npe is not null and \a project is defined, the solution is also
  //! projected onto the spline basis, and then evaluated at the \a npe points.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const int* npe = nullptr, char project = 0) const;

  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  Go::SplineCurve* projectSolution(const IntegrandBase& integrand) const;
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual Go::GeomObject* evalSolution(const IntegrandBase& integrand) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
			    const RealArray* gpar, bool = true) const;

  //! \brief Calculates parameter values for visualization nodal points.
  //! \param[out] prm Parameter values for all points
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  virtual bool getGridParameters(RealArray& prm, int nSegSpan) const;

  //! \brief Evaluates the basis functions at the specified point.
  //! \param[in] u Parameter value of evaluation point
  //! \param[out] N Basis function values
  virtual void evaluateBasis(double u, double, double, Vector& N) const;
  //! \brief Establishes the vector with basis function values.
  //! \param[in] u Parameter value of current integration point
  //! \param[out] N Basis function values
  void extractBasis(double u, Vector& N) const;
  //! \brief Establishes matrices with basis functions and 1st derivatives.
  //! \param[in] u Parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  void extractBasis(double u, Vector& N, Matrix& dNdu) const;
  //! \brief Establishes matrices with basis functions, 1st and 2nd derivatives.
  //! \param[in] u Parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  //! \param[out] d2Ndu2 Second derivatives of basis functions
  void extractBasis(double u, Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2) const;
  //! \brief Establishes matrices with basis functions, 1st, 2nd and 3rd derivatives.
  //! \param[in] u Parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  //! \param[out] d2Ndu2 Second derivatives of basis functions
  //! \param[out] d3Ndu3 Third derivatives of basis functions
  void extractBasis(double u, Vector& N, Matrix& dNdu,
                    Matrix3D& d2Ndu2, Matrix4D& d3Ndu3) const;

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

  //! \brief Initializes the local element axes for a patch of beam elements.
  //! \param[in] Zaxis Vector defining a point in the local XZ-plane
  bool initLocalElementAxes(const Vec3& Zaxis);

  //! \brief Connects matching nodes on two adjacent vertices.
  //! \param[in] vertex Local vertex index of this patch, in range [1,2]
  //! \param neighbor The neighbor patch
  //! \param[in] nvertex Local vertex index of neighbor patch, in range [1,2]
  //! \param[in] basis Which basis to connect the nodes for (mixed methods)
  //! \param[in] slave 0-based index of the first slave node in this basis
  //! \param[in] master 0-based index of the first master node in this basis
  //! \param[in] thick Thickness of connection
  bool connectBasis(int vertex, ASMs1D& neighbor, int nvertex,
                    int basis = 1, int slave = 0, int master = 0,
                    int thick = 1);

  //! \brief Extracts parameter values of the Gauss points.
  //! \param[out] uGP Parameter values for all points
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  //! \return The parameter value matrix casted into a one-dimensional vector
  const Vector& getGaussPointParameters(Matrix& uGP, int nGauss,
					const double* xi) const;

  //! \brief Calculates parameter values for the Greville points.
  //! \param[out] prm Parameter values for all points
  bool getGrevilleParameters(RealArray& prm) const;

  //! \brief Returns the length in the parameter space for an element.
  //! \param[in] iel Element index
  double getParametricLength(int iel) const;

  //! \brief Returns the parametric length on the \a i'th knot-span.
  double getKnotSpan(int i) const;

  //! \brief Computes the element border parameters.
  //! \param[in] iel 1-based element index
  //! \param[out] u Parameter values of the element borders
  virtual void getElementBorders(int iel, double* u) const;

  //! \brief Computes the element end coordinates.
  //! \param[in] i Parameter index for the knot-span element
  //! \param[out] XC Coordinates of the element corners
  //! \return Element length
  double getElementEnds(int i, std::vector<Vec3>& XC) const;

  //! \brief Returns nodal rotation matrices for an element, if any.
  //! \param[out] T Array of nodal rotation matrices
  //! \param[in] iel 0-based element index
  bool getElementNodalRotations(TensorVec& T, size_t iel) const;

public:
  //! \brief Auxilliary function for computation of basis function indices.
  static void scatterInd(int p1, int start, IntVec& index);

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] p1 Order in first (u) direction
  //! \param[out] p2 Order in second (v) direction (always zero)
  //! \param[out] p3 Order in third (w) direction (always zero)
  virtual bool getOrder(int& p1, int& p2, int& p3) const;

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction (always zero)
  //! \param[out] n3 Number of nodes in third (w) direction (always zero)
  //! \param[in] basis Which basis to return size parameters for (mixed methods)
  virtual bool getSize(int& n1, int& n2, int& n3, int basis) const;

  //! \brief Returns the number of elements in each parameter direction.
  //! \param[out] n1 Number of elements in first (u) direction
  //! \param[out] n2 Number of elements in second (v) direction (always zero)
  //! \param[out] n3 Number of elements in third (w) direction (always zero)
  virtual bool getNoStructElms(int& n1, int& n2, int& n3) const;

  //! \brief Returns parameter values and node numbers of the domain corners.
  //! \param[out] u Parameter values of the domain corners
  //! \param[out] corners 1-based indices of the corner nodes (optional)
  virtual bool getParameterDomain(Real2DMat& u, IntVec* corners) const;

protected:
  Go::SplineCurve* curv; //!< Pointer to the actual spline curve object

  const TensorVec& elmCS;  //!< Element coordinate systems (for 3D beams)
  const TensorVec& nodalT; //!< Nodal rotation tensors (for 3D beams)

  TensorVec myCS;  //!< The actual element coordinate systems
  TensorVec myT;   //!< The actual nodal rotation tensors
  TensorVec prevT; //!< Nodal rotation tensors of last converged configuration
};

#endif
