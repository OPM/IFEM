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
#include "Profiler.h"
#include <fstream>
#include <set>


void LR::extendControlPoints (LR::LRSpline* basis, Vector& v, int nf, int ofs)
{
  std::vector<double> cpts;
  for (auto it = basis->basisBegin(); it != basis->basisEnd(); ++it) {
    int id = (*it)->getId();
    std::vector<double> cpp;
    (*it)->getControlPoint(cpp);
    cpts.insert(cpts.end(), cpp.begin(), cpp.end());
    cpts.insert(cpts.end(), v.begin()+id*nf+ofs, v.begin()+(id+1)*nf+ofs);
  }
  basis->rebuildDimension(basis->dimension()+nf);
  basis->setControlPoints(cpts);
}


void LR::contractControlPoints (LR::LRSpline* basis, Vector& v, int nf, int ofs)
{
  std::vector<double> cpts;
  for (auto it = basis->basisBegin(); it != basis->basisEnd(); ++it) {
    int id = (*it)->getId();
    std::vector<double> cpp;
    (*it)->getControlPoint(cpp);
    for (int i=1; i <= nf; ++i)
      v(id*nf+ofs+i) = cpp[basis->dimension()-nf+i-1];
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


ASMunstruct::~ASMunstruct ()
{
  if (geo) delete geo;
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

    for (size_t j = 0; j < sol.size(); j++)
      LR::extendControlPoints(geo, sol[j], this->getNoFields(1));

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
      sol[i].resize(this->getNoFields(1)*geo->nBasisFunctions());
      LR::contractControlPoints(geo, sol[i], this->getNoFields(1));
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


Go::BsplineBasis ASMunstruct::getBezierBasis (int p)
{
  double knot[2*p];
  std::fill(knot,   knot+p,  -1.0);
  std::fill(knot+p, knot+2*p, 1.0);
  return Go::BsplineBasis(p,p,knot);
}


IntVec ASMunstruct::getFunctionsForElements (const IntVec& elements)
{
  geo->generateIDs();
  std::set<int> functions; // to get unique function IDs
  for (const int iel : elements)
    for (LR::Basisfunction* b : geo->getElement(iel)->support())
      functions.insert(b->getId());

  IntVec result(functions.size());
  std::copy(functions.begin(), functions.end(), result.begin());

  return result;
}
