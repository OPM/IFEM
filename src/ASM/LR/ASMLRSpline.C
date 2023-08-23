// $Id$
//==============================================================================
//!
//! \file ASMLRSpline.C
//!
//! \date December 2010
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Basisfunction.h"

#include "ASMLRSpline.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "ThreadGroups.h"
#include "Utilities.h"
#include "Profiler.h"
#include "IFEM.h"
#include <fstream>

#ifdef USE_OPENMP
#include <omp.h>
#endif


int LR::extendControlPoints (LRSpline* basis, const Vector& v, int nf, int ofs)
{
  if (v.empty())
    return 0;
  else if (v.size() == (size_t)basis->nBasisFunctions())
    nf = 1; // This is a scalar field
  else if (v.size() != (size_t)(basis->nBasisFunctions()*nf)) {
    std::cerr <<" *** LR::extendControlPoints: Invalid vector size "
              << v.size() <<", nBasis = "<< basis->nBasisFunctions()
              << std::endl;
    return -1;
  }

  RealArray cpts, cpp;
  for (LR::Basisfunction* b : basis->getAllBasisfunctions()) {
    int istart = ofs + b->getId()*nf;
    b->getControlPoint(cpp);
    cpts.insert(cpts.end(), cpp.begin(), cpp.end());
    cpts.insert(cpts.end(), v.begin()+istart, v.begin()+istart+nf);
  }
  basis->rebuildDimension(basis->dimension()+nf);
  basis->setControlPoints(cpts);
  return nf;
}


void LR::contractControlPoints (LRSpline* basis, Vector& v, int nf, int ofs)
{
  RealArray cpts, cpp;
  for (LR::Basisfunction* b : basis->getAllBasisfunctions()) {
    int istart = ofs + b->getId()*nf;
    b->getControlPoint(cpp);
    for (int i = 0; i < nf; i++)
      v[istart+i] = cpp[basis->dimension()-nf+i];
    cpts.insert(cpts.end(), cpp.begin(), cpp.end()-nf);
  }
  basis->rebuildDimension(basis->dimension()-nf);
  basis->setControlPoints(cpts);
}


void LR::getGaussPointParameters (const LRSpline* lrspline, RealArray& uGP,
                                  int d, int nGauss, int iel, const double* xi)
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** LR::getGaussPointParameters: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return;
  }
#endif

  const LR::Element* el = lrspline->getElement(iel-1);
  double ustart = el->getParmin(d);
  double ustop  = el->getParmax(d);

  uGP.resize(nGauss);
  for (int i = 0; i < nGauss; i++)
    uGP[i] = 0.5*((ustop-ustart)*xi[i] + ustop+ustart);
}


void LR::generateThreadGroups (ThreadGroups& threadGroups,
                               const LRSpline* lr,
                               const std::vector<LRSpline*>& addConstraints)
{
  int nElement = lr->nElements();
#ifdef USE_OPENMP
  if (omp_get_max_threads() > 1)
  {
    threadGroups[0].clear();
    threadGroups[1].clear();

    IntVec status(nElement,0); // status vector for elements:

    std::vector<IntSet> additionals;
    if (!addConstraints.empty()) {
      additionals.resize(nElement);
      for (LR::Element* e : lr->getAllElements()) {
        for (LRSpline* lr2 : addConstraints) {
           int elB = lr2->getElementContaining(e->midpoint());
           for (LR::Basisfunction* b2 : lr2->getElement(elB)->support())
             for (LR::Element* el3 : b2->support()) {
               RealArray midpoint = el3->midpoint();
               std::vector<RealArray> points;
               if (lr2->nElements() != lr->nElements()) {
                 RealArray diff(midpoint.size());
                 for (size_t j = 0; j < midpoint.size(); ++j)
                   diff[j] = (el3->getParmax(j) - el3->getParmin(j)) / 4.0;
                 if (lr->dimension() == 2)
                   points = {{midpoint[0] + diff[0], midpoint[1] + diff[1]},
                             {midpoint[0] - diff[0], midpoint[1] - diff[1]},
                             {midpoint[0] - diff[0], midpoint[1] + diff[1]},
                             {midpoint[0] + diff[0], midpoint[1] - diff[1]}};
                 else
                   points = {{midpoint[0] + diff[0], midpoint[1] - diff[1], midpoint[2] - diff[2]},
                             {midpoint[0] + diff[0], midpoint[1] + diff[1], midpoint[2] - diff[2]},
                             {midpoint[0] + diff[0], midpoint[1] + diff[1], midpoint[2] + diff[2]},
                             {midpoint[0] - diff[0], midpoint[1] - diff[1], midpoint[2] - diff[2]},
                             {midpoint[0] - diff[0], midpoint[1] + diff[1], midpoint[2] - diff[2]},
                             {midpoint[0] - diff[0], midpoint[1] + diff[1], midpoint[2] + diff[2]},
                             {midpoint[0] + diff[0], midpoint[1] - diff[1], midpoint[2] + diff[2]},
                             {midpoint[0] - diff[0], midpoint[1] - diff[1], midpoint[2] + diff[2]}};
               } else
                 points = {midpoint};

               for (const RealArray& vec : points)
                 additionals[e->getId()].insert(lr->getElementContaining(vec));
             }
        }
      }
    }

    // -1 is unusable for current color, 0 is available,
    // any other value is the assigned color

    int fixedElements = 0;
    for (int nColors = 0; fixedElements < nElement; nColors++)
    {
      // reset un-assigned element tags
      for (int i=0; i<nElement; i++)
        if (status[i]<0)
          status[i] = 0;

      // look for available elements
      IntVec thisColor;
      for (LR::Element* e : lr->getAllElements()) {
        int i = e->getId();
        if (status[i] == 0) {
          status[i] = nColors+1;
          thisColor.push_back(i);
          fixedElements++;
          for (LR::Basisfunction* b : e->support())
            for (LR::Element* el2 : b->support()) {
              int j = el2->getId();
              if (status[j] == 0)  // if not assigned a color yet
                status[j] = -1; // set as unavailable (with current color)
              if (static_cast<size_t>(j) < additionals.size())
                for (int extra : additionals[j])
                  if (status[extra] == 0)
                    status[extra] = -1;
            }
        }
      }
      threadGroups[0].push_back(thisColor);
    }
    return;
  }
#endif

  threadGroups.oneGroup(nElement); // No threading, all elements in one group
}


ASMLRSpline::ASMLRSpline (unsigned char n_p, unsigned char n_s,
                          unsigned char n_f)
  : ASMbase(n_p,n_s,n_f)
{
  geomB = nullptr;
}


ASMLRSpline::ASMLRSpline (const ASMLRSpline& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  geomB = patch.geomB;
}


bool ASMLRSpline::refine (const LR::RefineData& prm, Vectors& sol)
{
  if (!geomB)
    return false;

  if (shareFE && !prm.refShare)
  {
    // This patch shares spline object with another patch
    // (in another simulator on the same mesh),
    // and is assumed to have been refined already
    nnod = geomB->nBasisFunctions();
    nel  = geomB->nElements();
    return true;
  }

  if (prm.errors.empty() && prm.elements.empty())
    return true;

  PROFILE2("ASMLRSpline::refine");

  IntVec nf(sol.size());
  for (size_t j = 0; j < sol.size(); j++)
    if ((nf[j] = LR::extendControlPoints(geomB.get(),sol[j],this->getNoFields(1))) < 0)
      return false;

  if (!this->doRefine(prm,geomB.get()))
    return false;

  nnod = geomB->nBasisFunctions();
  nel  = geomB->nElements();
  IFEM::cout <<"Refined mesh: "<< nel <<" elements "<< nnod <<" nodes."<< std::endl;

  for (int i = sol.size()-1; i >= 0; i--)
    if (nf[i] > 0) {
      sol[i].resize(nf[i]*nnod);
      LR::contractControlPoints(geomB.get(),sol[i],nf[i]);
    }

  bool linIndepTest = prm.options.size() > 3 ? prm.options[3] != 0 : false;
  if (linIndepTest)
  {
    std::cout <<"Testing for linear independence by overloading"<< std::endl;
    bool isLinIndep = geomB->isLinearIndepByOverloading(false);
    if (!isLinIndep) {
      std::cout <<"Inconclusive..."<< std::endl;
#ifdef HAS_BOOST
      std::cout <<"Testing for linear independence by full tensor expansion"<< std::endl;
      isLinIndep = geomB->isLinearIndepByMappingMatrix(false);
#endif
    }
    if (isLinIndep)
      std::cout <<"...Passed."<< std::endl;
    else
      std::cout <<"FAILED!!!"<< std::endl;
    return isLinIndep;
  }

  return true;
}


bool ASMLRSpline::doRefine (const LR::RefineData& prm, LR::LRSpline* lrspline)
{
  // to pick up if LR splines get stuck while doing refinement,
  // print entry and exit point of this function

  char doRefine = 0;
  if (!prm.errors.empty())
    doRefine = 'E'; // Refine based on error indicators
  else if (!prm.elements.empty())
    doRefine = 'I'; // Refine the specified elements
  else
    return doRefine;

  double beta         = prm.options.size() > 0 ? prm.options[0]/100.0 : 0.10;
  int    multiplicity = prm.options.size() > 1 ? prm.options[1]       : 1;
  if (multiplicity > 1)
    for (int d = 0; d < lrspline->nVariate(); d++) {
      int p = lrspline->order(d) - 1;
      if (multiplicity > p) multiplicity = p;
    }

  enum refinementStrategy strat = LR_FULLSPAN;
  if (prm.options.size() > 2)
    switch (prm.options[2]) {
    case 1: strat = LR_MINSPAN; break;
    case 2: strat = LR_STRUCTURED_MESH; break;
    }

  int  maxTjoints     = prm.options.size() > 4 ? prm.options[4]      : -1;
  int  maxAspectRatio = prm.options.size() > 5 ? prm.options[5]      : -1;
  bool closeGaps      = prm.options.size() > 6 ? prm.options[6] != 0 : false;

  // set refinement parameters
  if (maxTjoints > 0)
    lrspline->setMaxTjoints(maxTjoints);
  if (maxAspectRatio > 0)
    lrspline->setMaxAspectRatio((double)maxAspectRatio);
  lrspline->setCloseGaps(closeGaps);
  lrspline->setRefMultiplicity(multiplicity);
  lrspline->setRefStrat(strat);

  // do actual refinement
  if (doRefine == 'E')
    lrspline->refineByDimensionIncrease(prm.errors,beta);
  else if (strat == LR_STRUCTURED_MESH)
    lrspline->refineBasisFunction(prm.elements);
  else
    lrspline->refineElement(prm.elements);

  lrspline->generateIDs();
  this->clear(true);

  return doRefine;
}


Go::BsplineBasis ASMLRSpline::getBezierBasis (int p, double start, double end)
{
  double knot[2*p];
  std::fill(knot,   knot+p,   start);
  std::fill(knot+p, knot+2*p, end);
  return Go::BsplineBasis(p,p,knot);
}


IntVec ASMLRSpline::getFunctionsForElements (const IntVec& elements,
                                             bool globalId) const
{
  IntSet functions; // to get unique function IDs
  this->getFunctionsForElements(functions,elements,globalId);
  return IntVec(functions.begin(), functions.end());
}


void ASMLRSpline::getFunctionsForElements (IntSet& functions,
                                           const IntVec& elements,
                                           bool globalId) const
{
  geomB->generateIDs();
  for (int elmId : elements)
  {
    int iel = globalId ? utl::findIndex(MLGE,1+elmId) : elmId;
    if (iel >= 0 && iel < geomB->nElements())
      for (LR::Basisfunction* b : geomB->getElement(iel)->support())
        functions.insert(globalId ? this->getNodeID(b->getId()+1)-1:b->getId());
  }
}


IntVec ASMLRSpline::getBoundaryCovered (const IntSet& nodes) const
{
  IntSet result;
  int numbEdges = (this->getNoParamDim() == 2) ? 4 : 6;
  for (int edge = 1; edge <= numbEdges; edge++)
  {
    IntVec oneBoundary; // 1-based list of boundary nodes
    this->getBoundaryNodes(edge,oneBoundary,1,1,0,true);
    for (int i : nodes)
      for (int j : oneBoundary)
        if (geomB->getBasisfunction(i)->contains(*geomB->getBasisfunction(j-1)))
          result.insert(j-1);
  }

  return IntVec(result.begin(), result.end());
}


IntVec ASMLRSpline::getOverlappingNodes (const IntSet& nodes, int dir) const
{
  IntSet result;
  for (int i : nodes)
  {
    LR::Basisfunction* b = geomB->getBasisfunction(i);
    for (LR::Element* el : b->support())
      for (LR::Basisfunction* basis : el->support())
      {
        bool support_only_bigger_in_allowed_direction = true;
        for (int j = 0; j < b->nVariate() && support_only_bigger_in_allowed_direction; j++)
        {
          if ((1<<j) & dir) continue; // the function is allowed to grow in the direction j
          if (b->getParmin(j) > basis->getParmin(j) ||
              b->getParmax(j) < basis->getParmax(j))
            support_only_bigger_in_allowed_direction = false;
        }
        if (support_only_bigger_in_allowed_direction)
          result.insert(basis->getId());
      }
  }

  return IntVec(result.begin(), result.end());
}


void ASMLRSpline::Sort (int u, int v, int orient,
                        std::vector<LR::Basisfunction*>& functions)
{
  std::sort(functions.begin(), functions.end(),
            [u,v,orient](const LR::Basisfunction* a, const LR::Basisfunction* b)
            {
              int i,p = a->getOrder(orient < 4 ? v : u);
              int idx = (orient & 4) ? v : u;
              for (i = 0; i <= p; i++)
                if ((*a)[idx][i] != (*b)[idx][i])
                  return (orient & 2) ? (*a)[idx][i] > (*b)[idx][i]
                                      : (*a)[idx][i] < (*b)[idx][i];

              p   = a->getOrder(orient < 4 ? u : v);
              idx = (orient & 4) ? u : v;
              for (i = 0; i <= p; i++)
                if ((*a)[idx][i] != (*b)[idx][i])
                  return (orient & 1) ? (*a)[idx][i] > (*b)[idx][i]
                                      : (*a)[idx][i] < (*b)[idx][i];

              return false;
            });
}


bool ASMLRSpline::transferCntrlPtVars (LR::LRSpline* oldBasis,
                                       const RealArray& oldVars,
                                       RealArray& newVars,
                                       int nGauss, int nf) const
{
  oldBasis->rebuildDimension(nf);
  oldBasis->setControlPoints(const_cast<RealArray&>(oldVars));
  return this->transferCntrlPtVars(oldBasis,newVars,nGauss);
}


std::pair<size_t,double> ASMLRSpline::findClosestNode (const Vec3& X) const
{
  double distance = 0.0;
  size_t inod = 0, iclose = 0;
  for (LR::Basisfunction* b : geomB->getAllBasisfunctions())
  {
    double d = (X-Vec3(&(*b->cp()),nsd)).length();
    if (++inod == 1 || d < distance)
    {
      iclose = inod;
      distance = d;
    }
  }

  return std::make_pair(iclose,distance);
}


Vec3 ASMLRSpline::getElementCenter (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > geomB->nElements())
  {
    std::cerr <<" *** ASMLRSpline::getElementCenter: Element index "<< iel
              <<" out of range [1,"<< geomB->nElements() <<"]."<< std::endl;
    return Vec3();
  }
#endif

  double u0[3] = { 0.0, 0.0, 0.0 };
  LR::Element* elm = geomB->getElement(iel-1);
  for (unsigned char d = 0; d < ndim; d++)
    u0[d] = 0.5*(elm->getParmin(d) + elm->getParmax(d));

  Vec3 XC;
  if (this->evalPoint(iel-1,u0,XC) < 0)
    std::cerr <<" *** ASMLRSpline::getElementCenter: Failed for element "
              << iel << std::endl;
  return XC;
}


bool ASMLRSpline::checkThreadGroups (const IntMat& groups,
                                     const std::vector<const LR::LRSpline*>& bases,
                                     const LR::LRSpline* threadBasis)
{
  bool ok = true;
  for (size_t gId = 1; gId <= groups.size(); gId++)
    for (size_t bId = 1; bId <= bases.size(); bId++)
    {
      IntSet nodes;
      const LR::LRSpline* basis = bases[bId-1];
      for (int elm : groups[gId-1]) {
        RealArray midpoint = threadBasis->getElement(elm)->midpoint();
        int bElm = basis->getElementContaining(midpoint);
        for (const LR::Basisfunction* func : basis->getElement(bElm)->support())
          if (!nodes.insert(func->getId()).second) {
            std::cerr <<" *** ASMLRSpline::checkThreadGroups: Function "
                      << func->getId() <<" on basis "<< bId
                      <<" is present for multiple elements in group "<< gId
                      << std::endl;
            ok = false;
          }
      }
    }

  return ok;
}


void ASMLRSpline::analyzeThreadGroups (const IntMat& groups)
{
  size_t min = std::numeric_limits<size_t>::max() - 1;
  size_t max = 0;
  std::vector<size_t> groupSizes;
  double avg = 0.0;
  for (const IntVec& group : groups) {
    min = std::min(group.size(), min);
    max = std::max(group.size(), max);
    groupSizes.push_back(group.size());
    avg += group.size();
  }
  avg /= groups.size();
  size_t half = groupSizes.size() / 2;
  std::nth_element(groupSizes.begin(), groupSizes.begin() + half, groupSizes.end());
  IFEM::cout << "\n Elements are divided in " << groups.size() << " colors "
             << "(min = " << min
             << ", max = " << max
             << ", avg = " << avg
             << ", med = " << groupSizes[half]
             << ")." << std::endl;
}


bool ASMLRSpline::getParameterDomain (Real2DMat& u, IntVec* corners) const
{
  u.resize(geomB->nVariate(),RealArray(2));
  for (int i = 0; i < geomB->nVariate(); ++i) {
    u[i][0] = geomB->startparam(i);
    u[i][1] = geomB->endparam(i);
  }

  if (corners) {
    if (geomB->nVariate() == 2) {
      for (int J = 0; J < 2; ++J)
        for (int I = 0; I < 2; ++I) {
          std::vector<LR::Basisfunction*> funcs;
          int dir = (I > 0 ? LR::EAST : LR::WEST) |
                    (J > 0 ? LR::NORTH : LR::SOUTH);
          geomB->getEdgeFunctions(funcs, static_cast<LR::parameterEdge>(dir));
          corners->push_back(funcs.front()->getId()+1);
        }
    } else if (geomB->nVariate() == 3) {
      for (int K = 0; K < 2; ++K)
        for (int J = 0; J < 2; ++J)
          for (int I = 0; I < 2; ++I) {
            std::vector<LR::Basisfunction*> funcs;
            int dir = (I > 0 ? LR::EAST  : LR::WEST)  |
                      (J > 0 ? LR::NORTH : LR::SOUTH) |
                      (K > 0 ? LR::TOP   : LR::BOTTOM);
            geomB->getEdgeFunctions(funcs, static_cast<LR::parameterEdge>(dir));
            corners->push_back(funcs.front()->getId()+1);
          }
    }
  }

  return true;
}


void ASMLRSpline::getNoIntPoints (size_t& nPt, size_t& nIPt)
{
  size_t nGp = 1;
  if (nGauss > 0 && nGauss <= 10)
    for (unsigned char d = 0; d < ndim; d++)
      nGp *= nGauss;
  else
  {
    // Use polynomial order to define number of quadrature points
    int ng[3] = { 0, 0, 0 };
    int nG1 = 0;
    this->getOrder(ng[0],ng[1],ng[2]);
    for (unsigned char d = 0; d < ndim && d < 3; d++)
      nG1 = std::max(nG1,ng[d]+nGauss%10);
    for (unsigned char d = 0; d < ndim; d++)
      nGp = nG1*nG1;
  }

  firstIp = nPt;
  nPt += nel*nGp;

  // Count additional interface quadrature points
  size_t nInterface = MLGE.size() - nel;
  if (nInterface > 0 && nInterface != nel && nGauss > 0 && nGauss <= 10)
    nIPt += nInterface*nGp/nGauss;
}
