// $Id$
//==============================================================================
//!
//! \file ASMu3Dmx.h
//!
//! \date Mar 8 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 3D spline mixed FE models.
//!
//==============================================================================

#ifndef _ASM_U3D_MX_H
#define _ASM_U3D_MX_H

#include "ASMu3D.h"
#include "ASMmxBase.h"


/*!
  \brief Driver for assembly of unstructured 3D spline mixed FE models.
  \details This class implements a three-field mixed formulation with splines as
  basis functions. The first field is of one order higher than the second field,
  and its basis is obtained by order-elevating the input spline object once.
  By default, the geometry is represented by the second (lower order) basis,
  however, by negating the \a n_f1 argument to the constructor, the first
  basis is used instead.
*/

class ASMu3Dmx : public ASMu3D, private ASMmxBase
{
public:
  //! \brief The constructor initializes the dimension of each basis.
  ASMu3Dmx(const CharVec& n_f);
  //! \brief Copy constructor.
  ASMu3Dmx(const ASMu3Dmx& patch, const CharVec& n_f = CharVec(3,0));
  //! \brief Empty destructor.
  virtual ~ASMu3Dmx() { lrspline = nullptr; geo = nullptr; }

  //! \brief Returns the spline surface representing the basis of this patch.
  virtual LR::LRSplineVolume* getBasis(int basis = 1);

  //! \brief Returns the spline surface representing the basis of this patch.
  virtual const LR::LRSplineVolume* getBasis(int basis = 1) const;

  // Methods for model generation
  // ============================

  //! \brief Writes the geometry/basis of the patch to given stream.
  virtual bool write(std::ostream& os, int basis = 0) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the node-to-IJ-index array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Returns the number of bases.
  virtual size_t getNoBasis() const { return m_basis.size(); }
  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;
  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int basis = 0) const;
  //! \brief Returns the number of DOFs per node.
  //! \param[in] inod 1-based node index local to current patch
  virtual unsigned char getNodalDOFs(size_t inod) const;
  //! \brief Returns the classification of a node.
  //! \param[in] inod 1-based node index local to current patch
  virtual char getNodeType(size_t inod) const;
  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int* sysMadof);

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

  //! \brief Extract the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  virtual bool getSolution(Matrix& sField, const Vector& locSol,
                           const IntVec& nodes) const;

  using ASMu3D::evalSolution;
  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] deriv Derivative order to return
  //! \param[in] nf If nonzero, evaluates nf fields on first basis
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size().
  //! Otherwise, we assume that it contains the \a u and \a v parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool = false,
                            int deriv = 0, int nf = 0) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! When \a regular is \e true, it is assumed that the parameter value array
  //! \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size().
  //! Otherwise, we assume that it contains the \a u and \a v parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const RealArray* gpar, bool = false) const;

  //! \brief Evaluates the projected solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \param[in] nf If nonzero, mixed evaluates nf fields on first basis
  virtual bool evalProjSolution(Matrix& sField, const Vector& locSol,
                                const int* npe, int nf = 0) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual void extractNodeVec(const Vector& globVec, Vector& nodeVec,
                              unsigned char = 0, int basis = 0) const;

  //! \brief Inject nodal results for this patch into a global vector.
  //! \param[in] nodeVec Nodal result vector for this patch
  //! \param[out] globVec Global solution vector in DOF-order
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual bool injectNodeVec(const Vector& nodeVec, Vector& globVec,
                             unsigned char = 0, int basis = 0) const;

  //! \brief Returns the number of projection nodes for this patch.
  virtual size_t getNoProjectionNodes() const;

  //! \brief Returns the number of refinement nodes for this patch.
  virtual size_t getNoRefineNodes() const;

  //! \brief Returns the number of refinement elements for this patch.
  virtual size_t getNoRefineElms() const;

  //! \brief Returns a field using the projection basis.
  //! \param[in] coefs The coefficients for the field
  //! \param[in] nf Number of components
  virtual Fields* getProjectedFields(const Vector& coefs, size_t nf) const;

  using ASMu3D::refine;
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the refinement
  //! \param sol Control point results values that are transferred to new mesh
  //! \param[in] fName Optional file name for an image of the resulting mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol,
                      const char* fName = nullptr);

  //! \brief Remap (geometry) element wise errors to refinement basis functions.
  //! \param     errors The remapped errors
  //! \param[in] origErr The element wise errors on the geometry mesh
  //! \param[in] elemErrors If true, map to elements instead of basis functions
  virtual void remapErrors(RealArray& errors,
                           const RealArray& origErr, bool elemErrors) const;

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag (see below)
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //!
  //! \details The face orientation flag \a norient must be in range [0,7].
  //! When interpreted as a binary number, its 3 digits are decoded as follows:
  //! - left digit = 1: The u and v parameters of the neighbor face are swapped
  //! - middle digit = 1: Parameter \a u in neighbor patch face is reversed
  //! - right digit = 1: Parameter \a v in neighbor patch face is reversed
  virtual bool connectPatch(int face, ASM3D& neighbor, int nface, int norient,
                            int = 0, bool coordCheck = true, int = 1);

  //! \brief Obtain the refinement basis.
  virtual const LR::LRSpline* getRefinementBasis() const;

protected:
  //! \brief Assembles L2-projection matrices for the secondary solution.
  //! \param[out] A Left-hand-side matrix
  //! \param[out] B Right-hand-side vectors
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e false, a discrete L2-projection is used
  virtual bool assembleL2matrices(SparseMatrix& A, StdVector& B,
                                  const IntegrandBase& integrand,
                                  bool continuous) const;

private:
  std::vector<std::shared_ptr<LR::LRSplineVolume>> m_basis; //!< Spline bases
  std::shared_ptr<LR::LRSplineVolume> projBasis; //!< Basis to project onto
  std::shared_ptr<LR::LRSplineVolume> refBasis; //!< Basis to refine based on
  const std::vector<Matrices>& bezierExtract;   //!< Bezier extraction matrices
  std::vector<Matrices>        myBezierExtract; //!< Bezier extraction matrices
};

#endif
