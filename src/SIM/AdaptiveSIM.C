// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.C
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for linear static FEM simulators.
//!
//==============================================================================

#include "AdaptiveSIM.h"
#include "SIMoutput.h"
#include "SIMenums.h"
#include "ASMunstruct.h"
#include "ASMmxBase.h"
#include "IntegrandBase.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <sstream>
#include <cstdio>


AdaptiveSIM::AdaptiveSIM (SIMoutput& sim, bool sa) : SIMadmin(sim), model(sim)
{
  alone  = sa;
  geoBlk = nBlock = 0;

  // Default grid adaptation parameters
  storeMesh    = false;
  linIndepTest = false;
  beta         = 10.0;
  errTol       = 1.0;
  rCond        = 0.0;
  condLimit    = 1.0e12;
  maxStep      = 10;
  maxDOFs      = 1000000;
  scheme       = 0;     // fullspan
  symmetry     = 1;     // no symmetry
  knot_mult    = 1;     // maximum regularity (continuity)
  threshold    = NONE;
  adaptor      = 0;
  adNorm       = 0;
  maxTjoints   = -1;
  maxAspRatio  = -1.0;
  closeGaps    = false;

  solution.resize(1);
}


AdaptiveSIM::~AdaptiveSIM ()
{
  for (const char* pfx : prefix)
    free(const_cast<char*>(pfx));
}


bool AdaptiveSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"adaptive"))
    return alone ? model.parse(elem) : true;

  const char* value = nullptr;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if ((value = utl::getValue(child,"maxstep")))
      maxStep = atoi(value);
    else if ((value = utl::getValue(child,"maxdof")))
      maxDOFs = atoi(value);
    else if ((value = utl::getValue(child,"errtol")))
      errTol = atof(value);
    else if ((value = utl::getValue(child,"maxcondition")))
      condLimit = atof(value);
    else if ((value = utl::getValue(child,"symmetry")))
      symmetry = atoi(value);
    else if ((value = utl::getValue(child,"knot_mult")))
      knot_mult = atoi(value);
    else if (!strcasecmp(child->Value(), "store_eps_mesh"))
      storeMesh = true; // no need for value here
    else if (!strcasecmp(child->Value(), "test_linear_independence"))
      linIndepTest = true; // no need for value here
    else if ((value = utl::getValue(child,"scheme"))) {
      if (!strcasecmp(value,"fullspan"))
        scheme = 0;
      else if (!strcasecmp(value,"minspan"))
        scheme = 1;
      else if (!strcasecmp(value,"isotropic_function"))
        scheme = 2;
      else
        std::cerr <<"  ** AdaptiveSIM::parse: Unknown refinement scheme \""
                  << value <<"\" (ignored)"<< std::endl;
      utl::getAttribute(child, "maxTjoints", maxTjoints);
      utl::getAttribute(child, "maxAspectRatio", maxAspRatio);
      utl::getAttribute(child, "closeGaps", closeGaps);
      IFEM::cout <<"\tRefinement scheme: "<< scheme << std::endl;
    }
    else if ((value = utl::getValue(child,"use_norm")))
    {
      adaptor = atoi(value);
      utl::getAttribute(child, "index", adNorm);
    }
    else if ((value = utl::getValue(child,"use_sub_norm")))
      adNorm = atoi(value);
    else if ((value = utl::getValue(child,"beta"))) {
      beta = atof(value);
      std::string type;
      utl::getAttribute(child, "type", type, true);
      if (type.compare("threshold") == 0 || type.compare("threashold") == 0)
        threshold = MAXIMUM;
      else if (type.compare("maximum") == 0)
        threshold = MAXIMUM;
      else if (type.compare("average") == 0)
        threshold = AVERAGE;
      else if (type.compare("minimum") == 0)
        threshold = MINIMUM;
      else if (type.compare("truebeta") == 0)
        threshold = TRUE_BETA;
      else if (type.compare("dorfel") == 0)
        threshold = DORFEL;
      IFEM::cout <<"\tRefinement percentage: "<< beta <<" type="<< threshold
                 << std::endl;
    }

  return true;
}


bool AdaptiveSIM::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ADAPTIVE",8))
  {
    std::istringstream cline(utl::readLine(is));

    cline >> beta >> errTol;
    if (cline.fail() || cline.bad()) return false;

    int itmp;
    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxStep = itmp;

    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxDOFs = itmp;

    cline >> itmp; // read knotline multiplicity
    if (!cline.fail() && !cline.bad())
      knot_mult = itmp;

    cline >> itmp; // read refinement scheme
    // (0=fullspan, 1=minspan, 2=structured mesh)
    if (!cline.fail() && !cline.bad())
      scheme = itmp;

    cline >> itmp; // read symmetry
    // (0=none, <n> always requests a multiplum of <n> elements for refinement)
    if (!cline.fail() && !cline.bad())
      symmetry = itmp;
  }
  else
    return model.parse(keyWord,is);

  return true;
}


void AdaptiveSIM::setAdaptationNorm (size_t normGroup, size_t normIdx)
{
  adaptor = normGroup;
  adNorm  = normIdx;
}


bool AdaptiveSIM::initAdaptor (size_t normGroup)
{
  if (normGroup > 0)
    adaptor = normGroup; // Override value from input file
  else if (adaptor == 0 && !model.haveAnaSol() && !opt.project.empty())
    adaptor = 1; // Set default to first projection if no analytical solution

  SIMoptions::ProjectionMap::const_iterator pit = opt.project.begin();
  for (size_t j = 1; pit != opt.project.end(); ++pit, j++)
    if (j == adaptor) break;

  IFEM::cout <<"\n\n >>> Starting adaptive simulation based on";
  if (pit != opt.project.end())
  {
    IFEM::cout <<"\n     "<< pit->second
               <<" error estimates (norm group "<< adaptor
               <<") <<<"<< std::endl;
    if (!adNorm) adNorm = 2; // Use norm 2, unless specified
  }
  else if (model.haveAnaSol())
  {
    adaptor = 0; // Assuming the exact errors are stored in the first group
    if (!adNorm) adNorm = 4; // Use norm 4, unless specified
    IFEM::cout <<" exact errors ";
    if (adNorm != 4)
      IFEM::cout <<"(norm index " << adNorm << ") ";
    IFEM::cout <<"<<<" << std::endl;
  }
  else
  {
    IFEM::cout <<" - nothing, can not do that! <<<"<< std::endl;
    return false;
  }

  projs.resize(opt.project.size());
  if (opt.format >= 0)
  {
    prefix.reserve(opt.project.size());
    for (pit = opt.project.begin(); pit != opt.project.end(); ++pit)
      prefix.push_back(strdup(pit->second.c_str()));
  }

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, int iStep, bool withRF,
                             std::streamsize precision)
{
  model.getProcessAdm().cout <<"\nAdaptive step "<< iStep << std::endl;
  if (iStep > 1)
  {
    SIMoptions oldOpt(opt);
    // Re-generate the FE model after the refinement
    model.clearProperties();
    // Caution: If we are using XML-input and have specified old command-line
    // options in order to override simulation options read from the input file,
    // those options will not be overridden here, so please don't do that..
    if (!model.read(inputfile) || !model.preprocess())
      return false;
    opt = oldOpt;
  }
  else if (storeMesh)
    // Output the initial grid to eps-file
    model.refine(LR::RefineData(),"mesh_001.eps");

  // Assemble the linear FE equation system
  model.setMode(SIM::STATIC,true);
  model.initSystem(opt.solver,1,model.getNoRHS(),0,withRF);
  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem())
    return false;

  // Solve the linear system of equations
  int printSol = 1;
  solution.resize(model.getNoRHS());
  for (size_t i = 0; i < solution.size(); i++)
    if (!model.solveSystem(solution[i],printSol,&rCond,"displacement",i==0,i))
      return false;
    else if (solution.size() > 2)
      printSol = 0; // Print summary only for the first two solutions

  eNorm.clear();
  gNorm.clear();

  // Project the secondary solution onto the splines basis
  size_t idx = 0;
  model.setMode(SIM::RECOVERY);
  for (const auto& prj : opt.project)
    if (prj.first <= SIMoptions::NONE)
      projs[idx++].clear(); // No projection for this norm group
    else if (!model.project(projs[idx++],solution.front(),prj.first))
      return false;

  if (msgLevel > 1 && !projs.empty())
    model.getProcessAdm().cout << std::endl;

  // Evaluate solution norms
  model.setQuadratureRule(opt.nGauss[1]);
  if (!model.solutionNorms(solution.front(),projs,eNorm,gNorm))
    return false;

  return model.dumpResults(solution.front(),0.0,
                           model.getProcessAdm().cout,true,precision);
}


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> DblIdx;


bool AdaptiveSIM::adaptMesh (int iStep, std::streamsize outPrec)
{
  ASMbase* thePatch = model.getPatch(1);
  if (!thePatch)
    return false; // No patches in the model, nothing to do here

  if (iStep < 2)
    return true; // No refinement in the first adaptive cycle

  std::streamsize oldPrec = outPrec > 0 ? IFEM::cout.precision(outPrec) : 0;
  this->printNorms();
  if (outPrec > 0) IFEM::cout.precision(oldPrec);

  if (adaptor >= gNorm.size() || adaptor >= eNorm.rows() || eNorm.cols() == 0)
    return false; // Refinement norm index out of range, or no element errors

  // Cannot adapt on exact errors without an exact solution
  if (adaptor == 0 && !model.haveAnaSol())
    return false; // Cannot adapt on exact errors without an exact solution

  // Define the reference norm
  double refNorm = 0.01*model.getReferenceNorm(gNorm,adaptor);
  if (refNorm < -epsZ) {
    std::cerr <<" *** AdaptiveSIM::adaptMesh: Negative reference norm."
              <<" Check orientation of your model."<< std::endl;
    return false; // Negative reference norm, modelling error?
  }
  else if (refNorm < epsZ)
    return false; // Zero reference norm, probably no load on the model

  // Check if further refinement is required
  if (iStep > maxStep || model.getNoDOFs() > (size_t)maxDOFs)
    return false; // Refinement cycle or model size limit reached
  else if (gNorm[adaptor](adNorm) < errTol*refNorm)
    return false; // Discretization error tolerance reached
  else if (1.0/rCond > condLimit)
  {
    IFEM::cout <<"\n  ** Terminating the adaptive cycles due to instability."
               <<"\n     The last condition number "<< 1.0/rCond
               <<" is higher than the limit "<< condLimit << std::endl;
    return false; // Condition number limit reached
  }

  // Calculate row index in eNorm of the error norm to adapt based on
  size_t i, eRow = adNorm;
  NormBase* norm = model.getNormIntegrand();
  for (i = 0; i < adaptor; i++)
    eRow += norm->getNoFields(i+1);
  delete norm;

  LR::RefineData prm;
  prm.options.reserve(8);
  prm.options.push_back(beta);
  prm.options.push_back(knot_mult);
  prm.options.push_back(scheme);
  prm.options.push_back(linIndepTest);
  prm.options.push_back(maxTjoints);
  prm.options.push_back(floor(maxAspRatio));
  prm.options.push_back(closeGaps);
  prm.options.push_back(threshold == TRUE_BETA ? 1 : 0);

  if (threshold == TRUE_BETA)
  {
    if (model.getNoPatches() > 1) {
      std::cerr <<" *** AdaptiveSIM::adaptMesh: True beta refinement"
                <<" is not available for multi-patch models."<< std::endl;
      return false;
    }
    IFEM::cout <<"\nRefining by increasing solution space by "<< beta
               <<" percent."<< std::endl;
    prm.errors.resize(thePatch->getNoRefineElms());
    static_cast<ASMunstruct*>(thePatch)->remapErrors(prm.errors,
                                                     eNorm.getRow(eRow), true);
    if (!storeMesh)
      return model.refine(prm);

    char fname[13];
    sprintf(fname,"mesh_%03d.eps",iStep);
    return model.refine(prm,fname);
  }

  std::vector<DblIdx> errors;
  if (scheme == 2) // use errors per function
  {
    if (thePatch->getNoRefineNodes() == thePatch->getNoNodes(1))
      errors.resize(model.getNoNodes(),DblIdx(0.0,0));
    else if (model.getNoPatches() == 1)
      errors.resize(thePatch->getNoRefineNodes(),DblIdx(0.0,0));
    else
    {
      std::cerr <<" *** AdaptiveSIM::adaptMesh: Multi-patch refinement"
                <<" is not available for mixed models."<< std::endl;
      return false;
    }

    for (i = 0; i < errors.size(); i++)
      errors[i].second = i;

    for (ASMbase* patch : model.getFEModel())
    {
      // extract element norms for this patch
      Vector locNorm(patch->getNoElms());
      for (i = 1; i <= patch->getNoElms(); ++i)
        locNorm(i) = eNorm(eRow, patch->getElmID(i));

      // remap from geometry basis to refinement basis
      Vector locErr(patch->getNoProjectionNodes());
      static_cast<ASMunstruct*>(patch)->remapErrors(locErr, locNorm);

      // insert into global error array
      for (i = 0; i < locErr.size(); ++i)
        if (model.getNoPatches() > 1)
          errors[patch->getNodeID(i+1)-1].first += locErr[i];
        else
          errors[i].first += locErr[i];
    }
  }
  else { // use errors per element
    if (model.getNoPatches() > 1) // not supported for multi-patch models
    {
      std::cerr <<" *** AdaptiveSIM::adaptMesh: Multi-patch refinement"
                <<" is available for isotropic_function only."<< std::endl;
      return false;
    }
    for (i = 0; i < eNorm.cols(); i++)
      errors.push_back(DblIdx(eNorm(eRow,1+i),i));
  }

  // Sort the elements in the sequence of decreasing errors
  std::sort(errors.begin(),errors.end(),std::greater<DblIdx>());

  // find the list of refinable elements/basisfunctions
  // the variable 'toBeRefined' contains one of the following:
  //   - list of elements to be refined (if fullspan or minspan)
  //   - list of basisfunctions to be refined (if structured mesh)
  double limit, sumErr = 0.0, curErr = 0.0;
  switch (threshold) {
  case MAXIMUM: // beta percent of max error (less than 100%)
    limit = errors.front().first * beta/100.0;
    break;
  case AVERAGE: // beta percent of avg error (typical 100%)
    for (const DblIdx& error : errors)
      sumErr += error.first;
    limit = sumErr/errors.size() * beta/100.0;
    break;
  case MINIMUM: // beta percent of min error (more than 100%)
    limit = errors.back().first  * beta/100.0;
    break;
  case DORFEL:
    limit = errors.back().first;
    for (const DblIdx& error : errors)
      sumErr += error.first;
    sumErr *= beta/100.0;
    for (const DblIdx& error : errors)
      if (curErr < sumErr)
        curErr += error.first;
      else
      {
        limit = error.first;
        break;
      }
    break;
  default:
    limit = 0.0;
  }

  size_t refineSize;
  if (threshold == NONE)
    refineSize = ceil(errors.size()*beta/100.0);
  else
    refineSize = std::upper_bound(errors.begin(), errors.end(), DblIdx(limit,0),
                                  std::greater_equal<DblIdx>())
               - errors.begin();

  IFEM::cout <<"\nRefining "<< refineSize
             << (scheme < 2 ? " elements" : " basis functions")
             <<" with errors in range ["<< errors[refineSize-1].first
             <<","<< errors.front().first <<"] ";

  switch (threshold) {
  case NONE:
    IFEM::cout << beta <<"% of all "
               << (scheme < 2 ? "elements" : "basis functions");
    break;
  case MAXIMUM:
    IFEM::cout << beta <<"% of max error ("<< limit << ")";
    break;
  case AVERAGE:
    IFEM::cout << beta <<"% of average error ("<< limit <<")";
    break;
  case MINIMUM:
    IFEM::cout << beta <<"% of min error ("<< limit <<")";
    break;
  case DORFEL:
    IFEM::cout << beta <<"% of total error ("<< limit <<")";
    break;
  default:
    break;
  }
  IFEM::cout << std::endl;
/*
  if (symmetry > 0) // Make refineSize a multiplum of 'symmetry'
    refineSize += (symmetry-refineSize%symmetry);
*/

  if (refineSize < 1 || refineSize > errors.size())
    return false;

  prm.elements.reserve(refineSize);
  for (i = 0; i < refineSize; i++)
    prm.elements.push_back(errors[i].second);

  // Now refine the mesh
  if (!storeMesh)
    return model.refine(prm);

  char fname[13];
  sprintf(fname,"mesh_%03d.eps",iStep);
  return model.refine(prm,fname);
}


void AdaptiveSIM::printNorms (size_t w) const
{
  model.printNorms(gNorm,w);

  // Print out the norm used for mesh adaptation
  size_t i, eRow = adNorm;
  if (adaptor > 0 && adaptor < gNorm.size())
  {
    NormBase* norm = model.getNormIntegrand();

    const Vector& aNorm = gNorm[adaptor];
    IFEM::cout <<"\nError estimate"
               << utl::adjustRight(w-14,norm->getName(adaptor+1,adNorm))
               << aNorm(adNorm)
               <<"\nRelative error (%) : "
               << 100.0*aNorm(adNorm)/model.getReferenceNorm(gNorm,adaptor);
    if (model.haveAnaSol() && adaptor != 0)
      IFEM::cout <<"\nEffectivity index  : "
                 << model.getEffectivityIndex(gNorm,adaptor,adNorm);
    IFEM::cout <<"\n"<< std::endl;

    // Calculate the row-index for the corresponding element norms
    for (i = 0; i < adaptor; i++)
      eRow += norm->getNoFields(i+1);

    delete norm;
  }

  if (eNorm.rows() < eRow || eNorm.cols() < 1)
    return;

  // Compute some additional error measures
  double avg_norm = eNorm(eRow,1);
  double min_err  = avg_norm;
  double max_err  = avg_norm;
  for (i = 2; i <= eNorm.cols(); i++)
  {
    avg_norm += eNorm(eRow,i);
    if (min_err > eNorm(eRow,i))
      min_err = eNorm(eRow,i);
    else if (max_err < eNorm(eRow,i))
      max_err = eNorm(eRow,i);
  }
  avg_norm /= eNorm.cols();

  double RMS_norm = 0.0;
  for (i = 1; i <= eNorm.cols(); i++)
    RMS_norm += pow(eNorm(eRow,i)-avg_norm,2.0);
  RMS_norm = sqrt(RMS_norm/eNorm.cols())/avg_norm;

  IFEM::cout <<  "Root mean square (RMS) of error      : "<< RMS_norm
             <<"\nMin element error                    : "<< min_err
             <<"\nMax element error                    : "<< max_err
             <<"\nAverage element error                : "<< avg_norm
             << std::endl;
}


bool AdaptiveSIM::writeGlv (const char* infile, int iStep)
{
  if (opt.format < 0) return true;

  // Write VTF-file with model geometry
  if (!model.writeGlvG(geoBlk, iStep == 1 ? infile : nullptr))
    return false;

  // Write boundary tractions, if any
  if (!model.writeGlvT(iStep,geoBlk,nBlock))
    return false;

  // Write Dirichlet boundary conditions
  if (!model.writeGlvBC(nBlock,iStep))
    return false;

  // Write solution fields
  if (!model.writeGlvS(solution.front(),iStep,nBlock))
    return false;

  // Write projected solution fields
  SIMoptions::ProjectionMap::const_iterator pit = opt.project.begin();
  for (size_t i = 0; i < projs.size(); i++, ++pit)
    if (!model.writeGlvP(projs[i],iStep,nBlock,100+10*i,pit->second.c_str()))
      return false;

  // Write element norms
  if (!model.writeGlvN(eNorm,iStep,nBlock,prefix.data()))
    return false;

  // Write state information
  return model.writeGlvStep(iStep,iStep,1);
}
