// $Id$
//==============================================================================
//!
//! \file ASMs2DLag.h
//!
//! \date Mar 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 2D %Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_LAG_H
#define _ASM_S2D_LAG_H

#include "ASMs2D.h"
#include "ASMLagBase.h"


/*!
  \brief Driver for assembly of structured 2D %Lagrange FE models.
  \details This class contains methods for structured 2D %Lagrange patches.
*/

class ASMs2DLag : public ASMs2D, protected ASMLagBase
{
protected:
  //! \brief Implementation of basis function cache.
  class BasisFunctionCache : public ASMs2D::BasisFunctionCache
  {
  public:
    //! \brief The constructor initializes the class.
    //! \param pch Patch the cache is for
    BasisFunctionCache(const ASMs2DLag& pch);

    //! \brief Constructor reusing quadrature info from another instance.
    //! \param cache Instance holding quadrature information
    //! \param b Basis to use
    BasisFunctionCache(const ASMs2D::BasisFunctionCache& cache, int b);

    //! \brief Empty destructor.
    virtual ~BasisFunctionCache() = default;

    //! \brief Obtain a single integration point parameter.
    //! \param dir Direction of for integration point
    //! \param el Element number in given direction
    //! \param gp Integration point in given direction
    //! \param reduced True to return parameter for reduced quadrature
    double getParam(int dir, size_t el, size_t gp, bool reduced) const override;

  protected:
    //! \brief Implementation specific initialization.
    bool internalInit() override;

    //! \brief Obtain global integration point index.
    //! \param[in] gp Integration point on element (0-indexed)
    size_t index(size_t, size_t gp, bool) const override { return gp; }

    //! \brief Calculates basis function info in a single integration point.
    //! \param[in] gp Integration point on element (0-indexed)
    //! \param[in] red If \e true, returns for reduced integration scheme
    BasisFunctionVals calculatePt(size_t, size_t gp, bool red) const override;

    //! \brief Calculates basis function info in all integration points.
    void calculateAll() override;

    //! \brief Configure quadrature points.
    void setupParameters() override;
  };

public:
  //! \brief Default constructor.
  explicit ASMs2DLag(unsigned char n_s = 2, unsigned char n_f = 2);
  //! \brief Special copy constructor for sharing of FE data.
  ASMs2DLag(const ASMs2DLag& patch, unsigned char n_f);
  //! \brief Default copy constructor copying everything.
  ASMs2DLag(const ASMs2DLag& patch);
  //! \brief Empty destructor.
  virtual ~ASMs2DLag() {}


  // Methods for model generation
  // ============================

  //! \brief Writes the FEM basis to given stream.
  virtual bool write(std::ostream& os, int = 0) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the nodal coordinate array, as well as global node and element numbers.
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

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  //! \param[in] iel Element index
  virtual bool getElementCoordinates(Matrix& X, int iel, bool = true) const;

  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X, bool = false) const;

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Returns the geometric center of an element.
  virtual Vec3 getElementCenter(int iel) const;

  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] dir Parameter direction defining the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which DOFs to constrain at each node along the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain edge for
  virtual void constrainEdge(int dir, bool open, int dof, int code, char basis);

  //! \brief Returns the number of projection nodes for this patch.
  virtual size_t getNoProjectionNodes() const;

  //! \brief Assembles L2-projection matrices for the secondary solution.
  //! \param[out] A Left-hand-side matrix
  //! \param[out] B Right-hand-side vectors
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e false, a discrete L2-projection is used
  virtual bool assembleL2matrices(SystemMatrix& A, SystemVector& B,
                                  const L2Integrand& integrand,
                                  bool continuous) const;

protected:
  //! \brief Assigned global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  //! \param[in] Xnod Coordinates of the node
  void setCoord(size_t inod, const Vec3& Xnod);

  //! \brief Finds the patch-local element numbers on a patch boundary.
  //! \param[out] elms Array of element numbers
  //! \param[in] lIndex Local index of the boundary edge
  virtual void findBoundaryElms(IntVec& elms, int lIndex, int = 0) const;

  //! \brief Finds the element containing specified parametric point.
  //! \details Optionally, the local coordinates of the point are calculated.
  int findElement(double u, double v,
                  double* xi = nullptr, double* eta = nullptr) const;

  //! \brief Evaluates an integral over an interior element domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] iel Element index
  //! \param cache Basis function cache associated with the element
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  bool integrateElm(Integrand& integrand, GlobalIntegral& glbInt, int iel,
                    ASMs2D::BasisFunctionCache& cache, const TimeDomain& time);

public:
  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ);


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


  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameters in range [0.0,1.0] of the point
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that is closest to the point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const;

  //! \brief Creates a quad element model of this patch for visualization.
  //! \param[out] grid The generated quadrilateral grid
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int*) const;

  using ASMs2D::evalSolution;
  //! \brief Evaluates the primary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector in DOF-order
  //! \param[in] n_f If nonzero, mixed evaluates \a n_f fields on first basis
  //!
  //! \details The number of visualization points is the same as the order of
  //! the %Lagrange elements by default.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int*, int n_f, bool) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //!
  //! \details If \a gpar is null, the nodal point values for the solution are
  //! returned if \a regular is \e false and the element center values are
  //! returned if \a regular is \e true. If \a gpar is not null, we assume that
  //! it contains the \a u and \a v parameters for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool regular,
                            int, int) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] elms List of elements to evaluate at
  //! \param[in] lpar Local parameter values of the result sampling points
  //!
  //! \details If \a lpar is null, the solution is evaluated at the center of
  //! each element. Otherwise, it is assumed to contain the local &xi; and &eta;
  //! parameters of the evaluation points, w.r.t. the provided elements.
  //! \note The resuls are appended as extra columns to output matrix \a sField.
  //! Therefore, any existing results is preserved.
  bool evalSolution(Matrix& sField, const Vector& locSol,
                    const IntVec& elms, const RealArray* lpar = nullptr) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //!
  //! \details The number of visualization points is the same as the order of
  //! the %Lagrange elements by default.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const int*, char) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //!
  //! \details We assume that the parameter value array \a gpar contains
  //! the \a u and \a v parameters directly for each sampling point.
  //! If \a gpar is null or empty and \a regular is \e true,
  //! the solution is instead evaluated at all element centers.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const RealArray* gpar, bool regular) const;

  //! \brief Evaluates and interpolates a field over a given geometry.
  //! \param[in] basis The basis of the field to evaluate
  //! \param[in] locVec The coefficients of the field to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum The basis to evaluate for (mixed)
  virtual bool evaluate(const ASMbase* basis, const Vector& locVec,
                        RealArray& vec, int basisNum) const;

  //! \brief Returns the polynomial order in each parameter direction.
  //! \param[out] pu Order in first (u) direction
  //! \param[out] pv Order in second (v) direction
  //! \param[out] pw Order in third (w) direction (always zero)
  virtual bool getOrder(int& pu, int& pv, int& pw) const;

  using ASMs2D::getSize;
  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  virtual bool getSize(int& n1, int& n2, int = 0) const;

  using ASMs2D::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  virtual void generateThreadGroups(const Integrand&, bool, bool);

  //! \brief Returns the matrix of nodal point correspondance for given basis.
  virtual IntMat getElmNodes(int basis) const;

  //! \brief Returns the number of elements on a boundary.
  virtual size_t getNoBoundaryElms(char lIndex, char ldim) const;

private:
  //! \brief Returns matrix of nodal point correspondance for a structured grid.
  static void createMNPC(size_t nx, size_t ny, int p1, int p2, IntMat& MNPC);

  //! \brief Evaluates a nodal solution field at specified point in an element.
  bool evalSolPt(int iel, double xi, double eta, size_t nCmp,
                 const Vector& pchSol, RealArray& ptSol, RealArray& N) const;

protected:
  size_t nx; //!< Number of nodes in first parameter direction
  size_t ny; //!< Number of nodes in second parameter direction
  int    p1; //!< Polynomial order in first parameter direction
  int    p2; //!< Polynomial order in second parameter direction
};

#endif
