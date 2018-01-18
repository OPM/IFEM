// $Id$
//==============================================================================
//!
//! \file ASMunstruct.h
//!
//! \date September 2011
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#ifndef _ASM_UNSTRUCT_H
#define _ASM_UNSTRUCT_H

#include "ASMbase.h"
#include "GoTools/geometry/BsplineBasis.h"

typedef std::set<int> IntSet; //!< General integer set

class ThreadGroups;


namespace LR //! Utilities for LR-splines.
{
  class Basisfunction;
  class LRSpline;

  /*!
    \brief A struct of data to control the mesh refinement.
    \details The \a options parameters have the following interpretation:
    options[0] is the beta percentage of elements to refine,
    options[1] is the knotline multiplicity (default 1),
    options[2] is the refinement scheme (default 0),
    (FULLSPAN=0, MINSPAN=1, ISOTROPIC_ELEMENTS=2, ISOTROPIC_FUNCTIONS=3),
    options[3] is nonzero if testing for linear independence at all iterations,
    options[4] is the maximum number of T-joints allowed in the model,
    options[5] is the maximum allowed parametric aspect ratio of an element,
    options[6] is one if all "gaps" are to be closed,
    options[7] is one if using "true beta".
  */
  struct RefineData
  {
    bool      refShare; //!< If \e true, force refinement of shared FE grids
    IntVec    options;  //!< Parameters used to control the refinement
    IntVec    elements; //!< 0-based indices of the elements to refine
    RealArray errors;   //!< List of error indicators for the elements

    //! \brief Default constructor.
    explicit RefineData(bool rs = false) : refShare(rs) {}
  };

  //! \brief Expands the basis coefficients of an LR-spline object.
  //! \param basis The spline object to extend
  //! \param[in] v The vector to append to the basis coefficients
  //! \param[in] nf Number of fields in the given vector
  //! \param[in] ofs Offset in vector
  int extendControlPoints(LRSpline* basis, const Vector& v,
                          int nf, int ofs = 0);

  //! \brief Contracts the basis coefficients of an LR-spline object.
  //! \param basis The spline object to contact
  //! \param[out] v Vector containing the extracted basis coefficients
  //! \param[in] nf Number of fields in the given vector
  //! \param[in] ofs Offset in vector
  void contractControlPoints(LRSpline* basis, Vector& v,
                             int nf, int ofs = 0);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[in] spline The LR-spline object to get parameter values for
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] d Parameter direction (0,1,2)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] iel Element index
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  void getGaussPointParameters(const LRSpline* spline, RealArray& uGP,
                               int d, int nGauss, int iel, const double* xi);

  //! \brief Generates thread groups for a LR-spline mesh.
  //! \param[out] threadGroups The generated thread groups
  //! \param[in] lr The LR-spline to generate thread groups for
  void generateThreadGroups(ThreadGroups& threadGroups, const LRSpline* lr);
}


/*!
  \brief Base class for unstructured spline-based FE assembly drivers.
  \details This class contains methods common for unstructured spline patches.
*/

class ASMunstruct : public ASMbase
{
protected:
  //! \brief The constructor sets the space dimensions.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMunstruct(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Copy constructor.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  ASMunstruct(const ASMunstruct& patch, unsigned char n_f = 0);

public:
  //! \brief Empty destructor.
  virtual ~ASMunstruct() {}

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geo == nullptr; }

  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  //! \param[in] fName Optional file name for an image of the resulting mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol,
                      const char* fName = nullptr);

  //! \brief Resets global element and node counters.
  static void resetNumbering() { gEl = gNod = 0; }

  using ASMbase::evalSolution;
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual LR::LRSpline* evalSolution(const IntegrandBase& integrand) const = 0;

  //! \brief Returns a Bezier basis of order \a p.
  static Go::BsplineBasis getBezierBasis(int p,
                                         double start = -1.0, double end = 1.0);

  //! \brief Returns a list of basis functions having support on given elements.
  IntVec getFunctionsForElements(const IntVec& elements,
                                 bool globalId = false) const;
  //! \brief Returns a list of basis functions having support on given elements.
  void getFunctionsForElements(IntSet& functions, const IntVec& elements,
                               bool globalId = true) const;

  //! \brief Sort basis functions based on local knot vectors.
  static void Sort(int u, int v, int orient,
                   std::vector<LR::Basisfunction*>& functions);

  //! \brief Returns all boundary functions that are covered by the given nodes.
  //! \param[in] nodes Set of (0-based) patch local node IDs
  //! \return 0-based node IDs for boundary functions whose support is
  //! completely covered by the union of the support of the input nodes
  IntVec getBoundaryNodesCovered(const IntSet& nodes) const;

  //! \brief Returns all functions whose support overlap with the input nodes.
  //! \param[in] nodes List of (0-based) patch local node IDs
  //! (typically requested by adaptive refinement)
  //! \param[in] dir 3-bit binary mask on which parameter directions are allowed
  //! to grow; i.e. bin(011)=dec(3) allows u-direction and v-direction to grow,
  //! default is bin(111)=dec(7) all directions
  //! \return 0-based node IDs for functions with overlapping support with
  //! the ones in boundary
  IntVec getOverlappingNodes(const IntSet& nodes, int dir = 7) const;

  //! \brief Returns all functions whose support overlap with the input node.
  //! \param[in] node 0-based patch local node ID
  //! \param[in] dir 3-bit binary mask for which parameter directions can grow
  //! \return 0-based node IDs for functions with overlapping support
  IntVec getOverlappingNodes(int node, int dir = 7) const
  {
    return this->getOverlappingNodes(IntSet(&node,(&node)+1),dir);
  }

  //! \brief Remaps element-wise errors from geometry mesh to refinement mesh.
  //! \param[out] errors The remapped errors
  //! \param[in] orig The element-wise errors on the geometry mesh
  //! \param[in] elemErrors If true, map to elements instead of basis functions
  virtual void remapErrors(RealArray& errors, const RealArray& orig,
                           bool elemErrors = false) const = 0;

  //! \brief Extends the refinement domain with information for neighbors.
  //! \param refineIndices List of basis functions to refine
  //! \param neighborIndices Basis functions to refine from neighbor patches
  virtual void extendRefinementDomain(IntSet& refineIndices,
                                      const IntSet& neighborIndices) const = 0;

  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVars(const LR::LRSpline* oldBasis,
                                   const RealArray& oldVar, RealArray& newVar,
                                   int nGauss) const = 0;
  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVarsN(const LR::LRSpline* oldBasis,
                                    const RealArray& oldVar, RealArray& newVar,
                                    int nGauss) const = 0;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferCntrlPtVars(const LR::LRSpline* oldBasis,
                                   RealArray& newVar, int nGauss) const = 0;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Control point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] nf Number of field components
  bool transferCntrlPtVars(LR::LRSpline* oldBasis,
                           const RealArray& oldVar, RealArray& newVar,
                           int nGauss, int nf = 1) const;

  //! \brief Refines the parametrization based on a mesh density function.
  //! \param[in] refC Mesh refinement criteria function
  //! \param[in] refTol Mesh refinement threshold
  virtual bool refine(const RealFunc& refC, double refTol) = 0;

protected:
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param lrspline The spline to perform adaptation for
  bool doRefine(const LR::RefineData& prm, LR::LRSpline* lrspline);

  LR::LRSpline* geo; //!< Pointer to the actual spline geometry object

  static int gEl;  //!< Global element counter
  static int gNod; //!< Global node counter
};

#endif
