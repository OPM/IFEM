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
  if (v.size() == (size_t)basis->nBasisFunctions())
    nf = 1; // This is a scalar field
  else if (v.size() != (size_t)(basis->nBasisFunctions()*nf)) {
    std::cerr <<" *** LR::extendControlPoints: Invalid vector size "
              << v.size() <<", nBasis = "<< basis->nBasisFunctions()
              << std::endl;
    return 0;
  }

  RealArray cpts, cpp;
  for (LR::Basisfunction* b : basis->getAllBasisfunctions()) {
    int id = b->getId();
    b->getControlPoint(cpp);
    cpts.insert(cpts.end(), cpp.begin(), cpp.end());
    cpts.insert(cpts.end(), v.begin()+id*nf+ofs, v.begin()+(id+1)*nf+ofs);
  }
  basis->rebuildDimension(basis->dimension()+nf);
  basis->setControlPoints(cpts);
  return nf;
}


void LR::contractControlPoints (LRSpline* basis, Vector& v, int nf, int ofs)
{
  RealArray cpts, cpp;
  for (LR::Basisfunction* b : basis->getAllBasisfunctions()) {
    int id = b->getId();
    b->getControlPoint(cpp);
    for (int i = 0; i < nf; i++)
      v[id*nf+ofs+i] = cpp[basis->dimension()-nf+i];
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
    std::cerr <<" *** getLRGaussPointParameters: Element index "<< iel
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
    for (int i = 0; i < 2; i++)
      threadGroups[i].clear();
    IntMat& answer = threadGroups[0];

    IntVec status(nElement,0); // status vector for elements:

    std::vector<std::set<int>> additionals;
    if (!addConstraints.empty()) {
      additionals.resize(nElement);
      for (auto e : lr->getAllElements()) {
        for (const LRSpline* lr2 : addConstraints) {
           int elB = lr2->getElementContaining(e->midpoint());
           for (auto b2 : lr2->getElement(elB)->support())
             for (auto el3 : b2->support()) {
               std::vector<double> midpoint = el3->midpoint();
               std::vector<std::vector<double>> points;
               if (lr2->nElements() != lr->nElements()) {
                 std::vector<double> diff(midpoint.size());
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

               for (const std::vector<double>& vec : points)
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
      for (auto e : lr->getAllElements()) {
        int i = e->getId();
        if (status[i] == 0) {
          status[i] = nColors+1;
          thisColor.push_back(i);
          fixedElements++;
          for (auto b : e->support()) // for all basisfunctions with support here
            for (auto el2 : b->support()) {// for all elements this function supports
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
      answer.push_back(thisColor);
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
  geo = nullptr;
}


ASMLRSpline::ASMLRSpline (const ASMLRSpline& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  geo = patch.geo;
}


bool ASMLRSpline::refine (const LR::RefineData& prm, Vectors& sol)
{
  PROFILE2("ASMLRSpline::refine()");

  if (!geo)
    return false;
  else if (shareFE && !prm.refShare)
  {
    nnod = geo->nBasisFunctions();
    return true;
  }
  else if (prm.errors.empty() && prm.elements.empty())
    return true;

  IntVec nf(sol.size());
  for (size_t j = 0; j < sol.size(); j++)
    if (!sol[j].empty())
      if (!(nf[j] = LR::extendControlPoints(geo,sol[j],this->getNoFields(1))))
        return false;

  if (!this->doRefine(prm,geo))
    return false;

  nnod = geo->nBasisFunctions();
  for (int i = sol.size()-1; i >= 0; i--)
    if (!sol[i].empty()) {
      sol[i].resize(nf[i]*geo->nBasisFunctions());
      LR::contractControlPoints(geo,sol[i],nf[i]);
    }

  IFEM::cout <<"Refined mesh: "<< geo->nElements() <<" elements "
             << geo->nBasisFunctions() <<" nodes."<< std::endl;

  bool linIndepTest = prm.options.size() > 3 ? prm.options[3] != 0 : false;
  if (linIndepTest)
  {
    std::cout <<"Testing for linear independence by overloading"<< std::endl;
    bool isLinIndep = geo->isLinearIndepByOverloading(false);
    if (!isLinIndep) {
      std::cout <<"Inconclusive..."<< std::endl;
#ifdef HAS_BOOST
      std::cout <<"Testing for linear independence by full tensor expansion"<< std::endl;
      isLinIndep = geo->isLinearIndepByMappingMatrix(false);
#endif
    }
    if (isLinIndep)
      std::cout <<"...Passed."<< std::endl;
    else {
      std::cout <<"FAILED!!!"<< std::endl;
      exit(228);
    }
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
  geo->generateIDs();
  for (int elmId : elements)
  {
    int iel = globalId ? utl::findIndex(MLGE,1+elmId) : elmId;
    if (iel >= 0 && iel < geo->nElements())
      for (LR::Basisfunction* b : geo->getElement(iel)->support())
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
    for (const int i : nodes)
      for (const int j : oneBoundary)
        if (geo->getBasisfunction(i)->contains(*geo->getBasisfunction(j-1)))
          result.insert(j-1);
  }

  return IntVec(result.begin(), result.end());
}


IntVec ASMLRSpline::getOverlappingNodes (const IntSet& nodes, int dir) const
{
  IntSet result;
  for (const int i : nodes)
  {
    LR::Basisfunction *b = geo->getBasisfunction(i);
    for (auto el : b->support()) // for all elements where *b has support
      for (auto basis : el->support()) // for all functions on this element
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
  for (LR::Basisfunction* b : geo->getAllBasisfunctions())
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


bool ASMLRSpline::checkThreadGroups (const IntMat& groups,
                                     const std::vector<const LR::LRSpline*> bases,
                                     const LR::LRSpline* threadBasis)
{
  size_t groupId = 1;
  bool ok = true;
  for (const auto& elms : groups) {
    size_t basisId = 1;
    for (const LR::LRSpline* basis : bases) {
      std::set<int> nodes;
      for (int elm : elms) {
        std::vector<double> midpoint = threadBasis->getElement(elm)->midpoint();
        int bElm = basis->getElementContaining(midpoint);
        for (const LR::Basisfunction* func : basis->getElement(bElm)->support()) {
          if (nodes.find(func->getId()) != nodes.end()) {
            std::cerr << " ** Warning: Function " << func->getId() << " on basis " << basisId << " present for multiple elements in group " << groupId << std::endl;
            ok = false;
          }
          nodes.insert(func->getId());
        }
      }
      ++basisId;
    }
    ++groupId;
  }

  return ok;
}


void ASMLRSpline::analyzeThreadGroups(const IntMat& groups)
{
  size_t min = std::numeric_limits<size_t>::max() - 1;
  size_t max = 0;
  std::vector<size_t> groupSizes;
  double avg = 0;
  for (const auto& group : groups) {
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
             << ").\n";

}
