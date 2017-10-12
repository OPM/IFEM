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

#include "ASMunstruct.h"
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Basisfunction.h"
#include "Profiler.h"
#include "ThreadGroups.h"
#include <fstream>
#include <numeric>
#include <set>

#ifdef USE_OPENMP
#include <omp.h>
#endif


int LR::extendControlPoints (LR::LRSpline* basis, const Vector& v,
                             int nf, int ofs)
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
  for (auto it = basis->basisBegin(); it != basis->basisEnd(); ++it) {
    int id = (*it)->getId();
    (*it)->getControlPoint(cpp);
    cpts.insert(cpts.end(), cpp.begin(), cpp.end());
    cpts.insert(cpts.end(), v.begin()+id*nf+ofs, v.begin()+(id+1)*nf+ofs);
  }
  basis->rebuildDimension(basis->dimension()+nf);
  basis->setControlPoints(cpts);
  return nf;
}


void LR::contractControlPoints (LR::LRSpline* basis, Vector& v,
                                int nf, int ofs)
{
  RealArray cpts, cpp;
  for (auto it = basis->basisBegin(); it != basis->basisEnd(); ++it) {
    int id = (*it)->getId();
    (*it)->getControlPoint(cpp);
    for (int i = 0; i < nf; i++)
      v[id*nf+ofs+i] = cpp[basis->dimension()-nf+i];
    cpts.insert(cpts.end(), cpp.begin(), cpp.end()-nf);
  }
  basis->rebuildDimension(basis->dimension()-nf);
  basis->setControlPoints(cpts);
}


void LR::getGaussPointParameters (const LR::LRSpline* lrspline, RealArray& uGP,
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


void LR::generateThreadGroups (ThreadGroups& threadGroups, const LR::LRSpline* lr)
{
  int nElement = lr->nElements();
  int nt = 1;
#ifdef USE_OPENMP
  nt = omp_get_max_threads();
#endif
  if (nt == 1) {
    threadGroups[0].resize(1);
    threadGroups[0][0].resize(nElement);
    std::iota(threadGroups[0][0].begin(), threadGroups[0][0].end(), 0);
  } else {
    std::vector<int> status(nElement, 0); // status vector for elements: -1 is unusable for current color, 0 is available, any other is assigned color
    std::vector<std::vector<int> > answer;
    int fixedElements = 0;
    int nColors       = 0;
    while(fixedElements < nElement) {
      // reset un-assigned element tags
      for (int i=0; i<nElement; i++)
        if (status[i]<0)
          status[i] = 0;
      std::vector<int> thisColor;

      // look for available elements
      for (auto e : lr->getAllElements() ) {
        int i = e->getId();
        if (status[i] == 0) {
          status[i] = nColors+1;
          thisColor.push_back(i);
          fixedElements++;
          for (auto b : e->support()) // for all basisfunctions with support here
            for (auto el2 : b->support()) {// for all element this function supports
              int j = el2->getId();
              if (status[j] == 0)  // if not assigned a color yet
                status[j] = -1; // set as unavailable (with current color)
            }
        }
      }
      answer.push_back(thisColor);
      nColors++;
    }

    threadGroups[0] = answer;
  }
}


ASMunstruct::~ASMunstruct ()
{
}


bool ASMunstruct::refine (const LR::RefineData& prm,
                          Vectors& sol, const char* fName)
{
  PROFILE2("ASMunstruct::refine()");

  if (!geo)
    return false;
  else if (shareFE && !prm.refShare)
  {
    nnod = geo->nBasisFunctions();
    return true;
  }

  // to pick up if LR splines get stuck while doing refinement,
  // print entry and exit point of this function

  double beta         = prm.options.size() > 0 ? prm.options[0]/100.0 : 0.10;
  int    multiplicity = prm.options.size() > 1 ? prm.options[1]       : 1;
  if (multiplicity > 1)
    for (int d = 0; d < geo->nVariate(); d++) {
      int p = geo->order(d) - 1;
      if (multiplicity > p) multiplicity = p;
    }

  enum refinementStrategy strat = LR_FULLSPAN;
  if (prm.options.size() > 2)
    switch (prm.options[2]) {
    case 1: strat = LR_MINSPAN; break;
    case 2: strat = LR_STRUCTURED_MESH; break;
    }

  bool linIndepTest   = prm.options.size() > 3 ? prm.options[3] != 0 : false;
  int  maxTjoints     = prm.options.size() > 4 ? prm.options[4]      : -1;
  int  maxAspectRatio = prm.options.size() > 5 ? prm.options[5]      : -1;
  bool closeGaps      = prm.options.size() > 6 ? prm.options[6] != 0 : false;

  char doRefine = 0;
  if (!prm.errors.empty())
    doRefine = 'E'; // Refine based on error indicators
  else if (!prm.elements.empty())
    doRefine = 'I'; // Refine the specified elements

  if (doRefine) {

    std::vector<int> nf(sol.size());
    for (size_t j = 0; j < sol.size(); j++)
      if (!(nf[j] = LR::extendControlPoints(geo,sol[j],this->getNoFields(1))))
        return false;

    // set refinement parameters
    if (maxTjoints > 0)
      geo->setMaxTjoints(maxTjoints);
    if (maxAspectRatio > 0)
      geo->setMaxAspectRatio((double)maxAspectRatio);
    geo->setCloseGaps(closeGaps);
    geo->setRefMultiplicity(multiplicity);
    geo->setRefStrat(strat);

    // do actual refinement
    if (doRefine == 'E')
      geo->refineByDimensionIncrease(prm.errors,beta);
    else if (strat == LR_STRUCTURED_MESH)
      geo->refineBasisFunction(prm.elements);
    else
      geo->refineElement(prm.elements);

    geo->generateIDs();
    nnod = geo->nBasisFunctions();

    for (int i = sol.size()-1; i >= 0; i--) {
      sol[i].resize(nf[i]*geo->nBasisFunctions());
      LR::contractControlPoints(geo,sol[i],nf[i]);
    }
  }

  if (fName)
  {
    char fullFileName[256];

    strcpy(fullFileName, "lrspline_");
    strcat(fullFileName, fName);
    std::ofstream lrOut(fullFileName);
    lrOut << *geo;
    lrOut.close();

    LR::LRSplineSurface* lr = dynamic_cast<LR::LRSplineSurface*>(geo);
    if (lr) {
      // open files for writing
      strcpy(fullFileName, "param_");
      strcat(fullFileName, fName);
      std::ofstream paramMeshFile(fullFileName);

      strcpy(fullFileName, "physical_");
      strcat(fullFileName, fName);
      std::ofstream physicalMeshFile(fullFileName);

      strcpy(fullFileName, "param_dot_");
      strcat(fullFileName, fName);
      std::ofstream paramDotMeshFile(fullFileName);

      strcpy(fullFileName, "physical_dot_");
      strcat(fullFileName, fName);
      std::ofstream physicalDotMeshFile(fullFileName);

      lr->writePostscriptMesh(paramMeshFile);
      lr->writePostscriptElements(physicalMeshFile);
      lr->writePostscriptFunctionSpace(paramDotMeshFile);
      lr->writePostscriptMeshWithControlPoints(physicalDotMeshFile);

      // close all files
      paramMeshFile.close();
      physicalMeshFile.close();
      paramDotMeshFile.close();
      physicalDotMeshFile.close();
    }
  }

  if (doRefine)
    std::cout <<"Refined mesh: "<< geo->nElements() <<" elements "
              << geo->nBasisFunctions() <<" nodes."<< std::endl;

  if (linIndepTest)
  {
    std::cout <<"Testing for linear independence by overloading "<< std::endl;
    bool isLinIndep = geo->isLinearIndepByOverloading(false);
    if (!isLinIndep) {
      std::cout <<"Inconclusive..."<< std::endl;
#ifdef HAS_BOOST
      std::cout <<"Testing for linear independence by full tensor expansion "<< std::endl;
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


Go::BsplineBasis ASMunstruct::getBezierBasis (int p, double start, double end)
{
  double knot[2*p];
  std::fill(knot,   knot+p,   start);
  std::fill(knot+p, knot+2*p, end);
  return Go::BsplineBasis(p,p,knot);
}


IntVec ASMunstruct::getFunctionsForElements (const IntVec& elements)
{
  geo->generateIDs();
  IntSet functions; // to get unique function IDs
  for (const int iel : elements)
    for (LR::Basisfunction* b : geo->getElement(iel)->support())
      functions.insert(b->getId());

  return IntVec(functions.begin(), functions.end());
}


IntVec ASMunstruct::getBoundaryNodesCovered (const IntSet& nodes) const
{
  IntSet result;
  int numbEdges = (this->getNoParamDim() == 2) ? 4 : 6;
  for (int edge = 1; edge <= numbEdges; edge++)
  {
    IntVec oneBoundary;
    this->getBoundaryNodes(edge, oneBoundary, 1, 1, 0, true); // this returns a 1-indexed list
    for (const int i : nodes)
      for (const int j : oneBoundary)
        if (geo->getBasisfunction(i)->contains(*geo->getBasisfunction(j-1)))
          result.insert(j-1);
  }

  return IntVec(result.begin(), result.end());
}


IntVec ASMunstruct::getOverlappingNodes (const IntSet& nodes, int dir) const
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


void ASMunstruct::Sort(int u, int v, int orient,
                       std::vector<LR::Basisfunction*>& functions)
{
  std::sort(functions.begin(), functions.end(),
            [u,v,orient](const LR::Basisfunction* a, const LR::Basisfunction* b)
            {
              int p1 = a->getOrder(u);
              int p2 = a->getOrder(v);

              int idx1 = orient & 4 ? v : u;
              int idx2 = orient & 4 ? u : v;

              for (int i = 0; i < 1 + (orient < 4 ? p2 : p1); ++i)
                if ((*a)[idx1][i] != (*b)[idx1][i])
                  return orient & 2 ? (*a)[idx1][i] > (*b)[idx1][i]
                                    : (*a)[idx1][i] < (*b)[idx1][i];
              for (int i = 0; i < 1 + (orient < 4 ? p1 : p2); ++i)
                if ((*a)[idx2][i] != (*b)[idx2][i])
                  return orient & 1 ? (*a)[idx2][i] > (*b)[idx2][i]
                                    : (*a)[idx2][i] < (*b)[idx2][i];

              return false;
            });
}
