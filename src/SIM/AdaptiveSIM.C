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
      adaptor = atoi(value);
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


bool AdaptiveSIM::initAdaptor (size_t indxProj)
{
  if (indxProj > 0)
    adaptor = indxProj; // override value from XML input

  SIMoptions::ProjectionMap::const_iterator pit = opt.project.begin();
  for (size_t j = 1; pit != opt.project.end(); pit++, j++)
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
    IFEM::cout <<" exact errors <<<"<< std::endl;
    adaptor = 0; // Assuming the exact errors are stored in the first group
    if (!adNorm) adNorm = 4; // Use norm 4, unless specified
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
      prefix.push_back(pit->second.c_str());
  }

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, int iStep, bool withRF)
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
  if (!model.solveMatrixSystem(solution,1))
    return false;

  eNorm.clear();
  gNorm.clear();

  // Project the secondary solution onto the splines basis
  model.setMode(SIM::RECOVERY);
  SIMoptions::ProjectionMap::const_iterator pit = opt.project.begin();
  for (size_t i = 0; pit != opt.project.end(); i++, pit++)
    if (!model.project(projs[i],solution.front(),pit->first))
      return false;

  if (msgLevel > 1 && !projs.empty())
    model.getProcessAdm().cout << std::endl;

  // Evaluate solution norms
  model.setQuadratureRule(opt.nGauss[1]);
  if (!model.solutionNorms(solution,projs,eNorm,gNorm))
    return false;

  return model.dumpResults(solution.front(),0.0,
                           model.getProcessAdm().cout,true,6);
}


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> DblIdx;


bool AdaptiveSIM::adaptMesh (int iStep)
{
  if (iStep < 2)
    return true;

  this->printNorms();

  if (adaptor >= gNorm.size() || adaptor >= eNorm.rows())
    return false;

  // Cannot adapt on exact errors without an exact solution
  if (adaptor == 0 && !model.haveAnaSol())
    return false;

  // Define the reference norm
  double refNorm = 0.01*model.getReferenceNorm(gNorm,adaptor);
  if (refNorm <= 0.0) return false;

  // Check if further refinement is required
  if (iStep > maxStep || model.getNoDOFs() > (size_t)maxDOFs) return false;
  if (eNorm.cols() < 1 || gNorm[adaptor](adNorm) < errTol*refNorm) return false;

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
    IFEM::cout <<"\nRefining by increasing solution space by "<< beta
               <<" percent."<< std::endl;
    prm.errors = eNorm.getRow(eRow);
    if (!storeMesh)
      return model.refine(prm);

    char fname[13];
    sprintf(fname,"mesh_%03d.eps",iStep);
    return model.refine(prm,fname);
  }

  std::vector<DblIdx> errors;
  if (scheme == 2)
  {
    // Sum up the total error over all supported elements for each function
    ASMbase* patch = model.getPatch(1);
    if (!patch) return false;
    IntMat::const_iterator eit;
    IntVec::const_iterator nit;
    for (i = 0; i < patch->getNoNodes(); i++) // Loop over basis functions
      errors.push_back(DblIdx(0.0,i));
    for (i = 1, eit = patch->begin_elm(); eit < patch->end_elm(); eit++, i++)
      for (nit = eit->begin(); nit < eit->end(); nit++)
        errors[*nit].first += eNorm(eRow,i);
  }
  else
    for (i = 0; i < eNorm.cols(); i++)
      errors.push_back(DblIdx(eNorm(eRow,1+i),i));

  // Sort the elements in the sequence of decreasing errors
  std::sort(errors.begin(),errors.end(),std::greater<DblIdx>());

  // find the list of refinable elements/basisfunctions
  // the variable 'toBeRefined' contains one of the following:
  //   - list of elements to be refined (if fullspan or minspan)
  //   - list of basisfunctions to be refined (if structured mesh)
  size_t refineSize;
  double limit, sumErr = 0.0;
  switch (threshold) {
  case MAXIMUM: // beta percent of max error (less than 100%)
    limit = errors.front().first * beta/100.0;
    break;
  case AVERAGE: // beta percent of avg error (typical 100%)
    for (size_t i = 0; i < errors.size(); i++)
      sumErr += errors[i].first;
    limit = sumErr/errors.size() * beta/100.0;
    break;
  case MINIMUM: // beta percent of min error (more than 100%)
    limit = errors.back().first  * beta/100.0;
    break;
  default:
    limit = 0.0;
  }
  if (threshold == NONE)
    refineSize = ceil(errors.size()*beta/100.0);
  else
    refineSize = std::upper_bound(errors.begin(), errors.end(),
                                  DblIdx(limit,0), std::greater_equal<DblIdx>())
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

    const Vector& fNorm = gNorm.front();
    const Vector& aNorm = gNorm[adaptor];
    IFEM::cout <<"\nError estimate"
               << utl::adjustRight(w-14,norm->getName(adaptor+1,adNorm))
               << aNorm(adNorm)
               <<"\nRelative error (%) : "
               << 100.0*aNorm(adNorm)/model.getReferenceNorm(gNorm,adaptor);
    if (model.haveAnaSol() && fNorm.size() > 3 && adNorm == 2)
      IFEM::cout <<"\nEffectivity index  : "<< aNorm(2)/fNorm(4);
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
  for (size_t i = 0; i < projs.size(); i++, pit++)
    if (!model.writeGlvP(projs[i],iStep,nBlock,100+10*i,pit->second.c_str()))
      return false;

  // Write element norms
  if (!model.writeGlvN(eNorm,iStep,nBlock,prefix.data()))
    return false;

  // Write state information
  return model.writeGlvStep(iStep,iStep,1);
}
