// $Id$
//==============================================================================
//!
//! \file ASMLRSpline.h
//!
//! \date December 2010
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for FE assembly drivers using LR B-splines.
//!
//==============================================================================

#ifndef _ASM_LR_SPLINE_H
#define _ASM_LR_SPLINE_H

#include "ASMbase.h"
#include "ASMunstruct.h"
#include "GoTools/geometry/BsplineBasis.h"

class ThreadGroups;


namespace LR //! Utilities for LR-splines.
{
  class Basisfunction;
  class LRSpline;

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
  //! \param[in] addConstraints If given, additional constraint bases
  void generateThreadGroups(ThreadGroups& threadGroups,
                            const LRSpline* lr,
                            const std::vector<LRSpline*>& addConstraints = {});
}


/*!
  \brief Base class for LR B-spline FE assembly drivers.
  \details This class contains methods common for unstructured spline patches.
*/

class ASMLRSpline : public ASMbase, public ASMunstruct
{
protected:
  //! \brief The constructor sets the space dimensions.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMLRSpline(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Special copy constructor for sharing of FE data.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  ASMLRSpline(const ASMLRSpline& patch, unsigned char n_f);

public:
  //! \brief Empty destructor.
  virtual ~ASMLRSpline() {}

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geomB == nullptr; }

  //! \brief Returns parameter values and node numbers of the domain corners.
  virtual bool getParameterDomain(Real2DMat&, IntVec*) const;

  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol);

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
  virtual IntVec getBoundaryCovered(const IntSet& nodes) const;

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

  //! \brief Finds the node that is closest to the given point \b X.
  virtual std::pair<size_t,double> findClosestNode(const Vec3& X) const;

  //! \brief Returns the coordinates of the element center.
  virtual Vec3 getElementCenter(int iel) const;

  //! \brief Computes the total number of integration points in this patch.
  virtual void getNoIntPoints(size_t& nPt, size_t& nIPt);

protected:
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param lrspline The spline to perform adaptation for
  bool doRefine(const LR::RefineData& prm, LR::LRSpline* lrspline);

  using ASMbase::evalPoint;
  //! \brief Evaluates the geometry at a specified point.
  virtual int evalPoint(int iel, const double* param, Vec3& X) const = 0;

  //! \brief Santity check thread groups.
  //! \param[in] groups The generated thread groups
  //! \param[in] bases The bases to check for
  //! \param[in] threadBasis The LRSpline the element groups are derived from
  static bool checkThreadGroups(const IntMat& groups,
                                const std::vector<const LR::LRSpline*>& bases,
                                const LR::LRSpline* threadBasis);

  //! \brief Analyze and print thread group statistics.
  //! \param[in] groups The generated thread groups
  static void analyzeThreadGroups(const IntMat& groups);

  std::shared_ptr<LR::LRSpline> geomB;    //!< Pointer to the spline geometry object
  std::shared_ptr<LR::LRSpline> projB;    //!< Pointer to the projection spline object
  std::shared_ptr<LR::LRSpline> altProjB; //!< Alternative projection basis
  std::shared_ptr<LR::LRSpline> refB;     //!< Basis to refine based on
};

#endif
