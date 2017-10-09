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
    (FULLSPAN=0, MINSPAN=1, ISOTROPIC ELEMENTS=2, ISOTROPIC FUNCTIONS=3),
    options[3] is the symmetry, i.e., always refine a multiple of this value,
    options[4] is nonzero if testing for linear independence at all iterations,
    options[5] is the maximum number of T-joints allowed in the model,
    options[6] is the maximum allowed parametric aspect ratio of an element,
    options[7] is one if all "gaps" are to be closed,
    options[8] is one if using true beta.
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
  //! \brief The destructor frees the dynamically allocated spline object.
  virtual ~ASMunstruct();

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
  IntVec getFunctionsForElements(const IntVec& elements);

  //! \brief Sort basis functions based on local knot vectors.
  static void Sort(int u, int v, int orient,
                   std::vector<LR::Basisfunction*>& functions);

  //! \brief Return all boundary functions that are covered by the given set of input nodes
  //! \param nodes List of (0-indexed) patch local node IDs (typically requested by adaptive refinement)
  //! \returns     Node IDs (0-indexed) for boundary functions whos support is
  //               completely covered by the union of the support in the input array
  virtual IntVec getBoundaryNodesCovered(const IntVec& nodes) const ;

  //! \brief Return all functions whos support overlap with the input functions
  //! \param nodes List of (0-indexed) patch local node IDs (typically requested by adaptive refinement)
  //! \returns     Node IDs (0-indexed) for functions with overlapping support with the ones in boundary
  virtual IntVec getOverlappingNodes(const IntVec& nodes) const ;

  //! \brief Remap element wise errors from geometry mesh to refinement mesh.
  //! \param     errors The remapped errors
  //! \param[in] origErr The element wise errors on the geometry mesh
  virtual void remapErrors(RealArray& errors, const RealArray& origErr) const = 0;

  //! \brief Match neighbours after refinement in multipatch models.
  //! \param neigh Neigbouring patch
  //! \param[in] midx Index of face/edge on this patch
  //! \param[in] sidx  Index of face/edge on neighbour
  //! \param[in] orient Orientation flag for connection
  virtual bool matchNeighbour(ASMunstruct* neigh,
                              int midx, int sidx, int orient) = 0;

protected:
  LR::LRSpline* geo; //!< Pointer to the actual spline geometry object

  static int gEl;  //!< Global element counter
  static int gNod; //!< Global node counter
};

#endif
