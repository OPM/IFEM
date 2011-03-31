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

namespace Go {
  class SplineCurve;
}


/*!
  \brief Driver for assembly of structured 1D spline FE models.
  \details This class contains methods common for structured 1D spline patches.
*/

class ASMs1D : public ASMstruct
{
public:
  //! \brief Constructor creating an instance by reading the given file.
  ASMs1D(const char* fileName, unsigned char n_s = 1, unsigned char n_f = 1);
  //! \brief Constructor creating an instance by reading the given input stream.
  ASMs1D(std::istream& is, unsigned char n_s = 1, unsigned char n_f = 1);
  //! \brief Default constructor creating an empty patch.
  ASMs1D(unsigned char n_s = 1, unsigned char n_f = 1);
  //! \brief Empty destructor.
  virtual ~ASMs1D() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! and the global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  virtual void clear();

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Creates an instance by reading the given input stream, \a is.
  bool read(std::istream& is);
  //! \brief Writes the geometry of the SplineCurve object to given stream.
  bool write(std::ostream& os) const;

  //! \brief Refine the parametrization by inserting extra knots.
  //! \param[in] xi Relative positions of added knots in each existing knot span
  bool refine(const RealArray& xi);
  //! \brief Refine the parametrization by inserting extra knots uniformly.
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  bool uniformRefine(int nInsert);
  //! \brief Raise the order of the SplineCurve object for this patch.
  //! \param[in] ru Number of times to raise the order
  bool raiseOrder(int ru);


  // Various methods for preprocessing of boundary conditions
  // ========================================================

  //! \brief Makes the two end vertices of the curve periodic.
  void closeEnds();

  //! \brief Constrains a node identified by a relative parameter value.
  //! \param[in] xi Parameter value along the curve
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The parameter value has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  void constrainNode(double xi, int dof = 123, int code = 0);

  //! \brief Connects matching nodes on two adjacent vertices.
  //! \param[in] vertex Local vertex index of this patch, in range [1,2]
  //! \param neighbor The neighbor patch
  //! \param[in] nvertex Local vertex index of neighbor patch, in range [1,2]
  bool connectPatch(int vertex, ASMs1D& neighbor, int nvertex);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());

  //! \brief Evaluates a boundary integral over a patch end.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the end point
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());


  // Post-processing methods
  // =======================

  //! \brief Creates a line element model of this patch for visualization.
  //! \param[out] grid The generated line grid
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
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
			    const RealArray* gpar, bool = true) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  //!
  //! \details If \a npe is NULL, the solution is evaluated at the Greville
  //! points and then projected onto the spline basis to obtain the control
  //! point values, which then are returned through \a sField.
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const int* npe = 0) const;

  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  Go::SplineCurve* projectSolution(const Integrand& integrand) const;
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual Go::GeomObject* evalSolution(const Integrand& integrand) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const RealArray* gpar, bool = true) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Calculates parameter values for the visualization nodal points.
  //! \param[out] prm Parameter values for all points
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  virtual bool getGridParameters(RealArray& prm, int nSegSpan) const;

  //! \brief Returns the length in the parameter space for an element.
  //! \param[in] iel Element index
  double getParametricLength(int iel) const;

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel Element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const;
  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const;

  //! \brief Returns the number of nodal points in the patch.
  virtual int getSize() const;

private:
  //! \brief Establishes vectors with basis functions and 1st derivatives.
  //! \param[in] u Parameter value of current integration point
  //! \param[out] N Basis function values
  //! \param[out] dNdu First derivatives of basis functions
  void extractBasis(double u, Vector& N, Matrix& dNdu) const;

  //! \brief Returns the parametric length on the \a i'th knot-span
  double getKnotSpan(int i) const;

protected:
  Go::SplineCurve* curv; //!< Pointer to the actual spline curve object
};

#endif
