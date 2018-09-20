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
  \details This class implements a two-field mixed formulation with splines as
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
  explicit ASMu3Dmx(const CharVec& n_f);
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

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&, int basis);
  //! \brief Writes the geometry/basis of the patch to given stream.
  virtual bool write(std::ostream& os, int basis) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the node-to-IJ-index array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry);

  //! \brief Returns the number of bases.
  virtual size_t getNoBasis() const { return m_basis.size(); }
  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis) const;
  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int basis) const;
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

  //! \brief Constrains all DOFs on a given boundary face.
  //! \param[in] dir Parameter direction defining the face to constrain
  //! \param[in] open If \e true, exclude all points along the face boundary
  //! \param[in] dof Which DOFs to constrain at each node on the face
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain face for (0 means check all)
  virtual void constrainFace(int dir, bool open, int dof,
                             int code, char basis);
  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] lEdge Local index [1,12] of the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which DOFs to constrain at each node along the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain edge for (0 means check all)
  virtual void constrainEdge(int lEdge, bool open, int dof,
                             int code, char basis);
  //! \brief Constrains a corner node identified by the three parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] K Parameter index in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain node for (0 means check all)
  virtual void constrainCorner(int I, int J, int K, int dof,
                               int code, char basis);


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

  //! \brief Extracts the primary solution field at the specified nodes.
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
  //! \param[in] nf If nonzero, evaluates \a nf fields on first basis
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size().
  //! Otherwise, we assume that it contains the \a u and \a v parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool,
                            int deriv, int nf) const;

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
                            const RealArray* gpar, bool) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual void extractNodeVec(const RealArray& globVec, RealArray& nodeVec,
                              unsigned char, int basis) const;

  //! \brief Injects nodal results for this patch into a global vector.
  //! \param[in] nodeVec Nodal result vector for this patch
  //! \param[out] globVec Global solution vector in DOF-order
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual bool injectNodeVec(const RealArray& nodeVec, RealArray& globVec,
                             unsigned char, int basis) const;

  //! \brief Returns the number of refinement nodes for this patch.
  virtual size_t getNoRefineNodes() const;

  //! \brief Returns the number of refinement elements for this patch.
  virtual size_t getNoRefineElms() const;

  using ASMu3D::refine;
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the refinement
  //! \param sol Control point results values that are transferred to new mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol);

  //! \brief Remaps element-wise errors from geometry mesh to refinement mesh.
  //! \param[out] errors The remapped errors
  //! \param[in] origErr The element-wise errors on the geometry mesh
  //! \param[in] elemErrors If \e true, map to elements and not basis functions
  virtual void remapErrors(RealArray& errors,
                           const RealArray& origErr, bool elemErrors) const;

  //! \brief Copies the refinement to another spline volume.
  //! \param basis Volume to copy refinement to
  //! \param[in] multiplicity Wanted multiplicity
  void copyRefinement(LR::LRSplineVolume* basis, int multiplicity) const;

  //! \brief Swaps between the main and alternative projection basis.
  virtual void swapProjectionBasis();

private:
  //! \brief Finds the elements and associted sizes at given parametric point.
  //! \param[in] param Parametric point to find elements at
  //! \param[out] elms List of elements on each basis containign given point
  //! \param[out] sizes List of element sizes (numer of element nodes)
  void getElementsAt(const RealArray& param,
                     std::vector<int>& elms,
                     std::vector<size_t>& sizes) const;

protected:
  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM If \e true, ignore global multipliers in sanity check
  void generateThreadGroups(const Integrand& integrand, bool silence,
                            bool ignoreGlobalLM);

private:
  typedef std::shared_ptr<LR::LRSplineVolume> SplinePtr; //!< Pointer to spline

  std::vector<SplinePtr> m_basis;      //!< All bases
  LR::LRSplineVolume*    threadBasis;  //!< Basis for thread groups
  SplinePtr              refBasis;     //!< Basis to refine based on
  SplinePtr              altProjBasis; //!< Alternative projection basis

  const std::vector<Matrices>& bezierExtractmx;  //!< Bezier extraction matrices
  std::vector<Matrices>        myBezierExtractmx; //!< Bezier extraction matrices

  ThreadGroups altProjThreadGroups; //!< Element groups for multi-threaded assembly - alternative projection basis
};

#endif
