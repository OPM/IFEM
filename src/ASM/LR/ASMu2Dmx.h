// $Id$
//==============================================================================
//!
//! \file ASMu2Dmx.h
//!
//! \date May 7 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D spline mixed FE models.
//!
//==============================================================================

#ifndef _ASM_U2D_MX_H
#define _ASM_U2D_MX_H

#include "ASMu2D.h"
#include "ASMmxBase.h"


/*!
  \brief Driver for assembly of unstructured 2D spline mixed FE models.
  \details This class implements a two-field mixed formulation with splines as
  basis functions. The first field is of one order higher than the second field,
  and its basis is obtained by order-elevating the input spline object once.
  By default, the geometry is represented by the second (lower order) basis,
  however, by negating the \a n_f1 argument to the constructor, the first
  basis is used instead.
*/

class ASMu2Dmx : public ASMu2D, private ASMmxBase
{
public:
  //! \brief The constructor initializes the dimension of each basis.
  ASMu2Dmx(unsigned char n_s, const CharVec& n_f);
  //! \brief Copy constructor.
  ASMu2Dmx(const ASMu2Dmx& patch, const CharVec& n_f = CharVec(2,0));
  //! \brief Empty destructor.
  virtual ~ASMu2Dmx() { lrspline = nullptr; geo = nullptr; }

  //! \brief Returns the spline surface representing the basis of this patch.
  virtual LR::LRSplineSurface* getBasis(int basis = 1);

  //! \brief Returns the spline surface representing the basis of this patch.
  virtual const LR::LRSplineSurface* getBasis(int basis = 1) const;

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

  //! \brief Evaluates an integral over element interfaces in the patch.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] iChk Object checking if an element interface has contributions
  virtual bool integrate(Integrand& integrand, GlobalIntegral& glbInt,
                         const TimeDomain& time, const ASM::InterfaceChecker& iChk);


  // Post-processing methods
  // =======================

  //! \brief Extract the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  virtual bool getSolution(Matrix& sField, const Vector& locSol,
                           const IntVec& nodes) const;

  using ASMu2D::evalSolution;
  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] deriv Derivative order to return
  //! \param[in] nf If nonzero, evaluate nf fields on first basis
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
                            const RealArray* gpar, bool = true) const;

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

  //! \brief Projects the secondary solution using a discrete global L2-norm.
  //! \param[out] sField Secondary solution field control point values
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e true, a continuous L2-projection is used
  virtual bool globalL2projection(Matrix& sField,
                                  const IntegrandBase& integrand,
                                  bool continuous = false) const;

  using ASMu2D::refine;
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the refinement
  //! \param sol Control point results values that are transferred to new mesh
  //! \param[in] fName Optional file name for an image of the resulting mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol,
                      const char* fName = nullptr);

  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  //! \param[in] basis Which basis to connect
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //! \param[in] thick Thickness of connection
  virtual bool connectPatch(int edge, ASM2D& neighbor, int nedge, bool revers,
                            int basis, bool coordCheck = true, int thick = 1);

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary edge
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (0 for all)
  //! \param[in] thick Thickness of connection
  //! \param[in] orient Orientation for returned boundary nodes
  //! \param[in] local If \e true return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int basis, int thick = 1,
                                int orient = 0, bool local = false) const;

  //! \brief Remap (geometry) element wise errors to refinement basis functions.
  //! \param     errors The remapped errors
  //! \param[in] origErr The element wise errors on the geometry mesh
  virtual void remapErrors(RealArray& errors, const RealArray& origErr) const;

protected:
  //! \brief Assembles L2-projection matrices for the secondary solution.
  //! \param[out] A Left-hand-side matrix
  //! \param[out] B Right-hand-side vectors
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] continuous If \e false, a discrete L2-projection is used
  virtual bool assembleL2matrices(SparseMatrix& A, StdVector& B,
                                  const IntegrandBase& integrand,
                                  bool continuous) const;

  using ASMu2D::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM If \e true, ignore global multipliers in sanity check
  void generateThreadGroups(const Integrand& integrand, bool silence,
                            bool ignoreGlobalLM);

private:
  std::vector<std::shared_ptr<LR::LRSplineSurface>> m_basis; //!< All bases
  LR::LRSplineSurface* threadBasis; //!< Basis for thread groups
};

#endif
