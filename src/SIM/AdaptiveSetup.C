// $Id$
//==============================================================================
//!
//! \file AdaptiveSetup.C
//!
//! \date May 7 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution setup for linear and nonlinear FEM simulators.
//!
//==============================================================================

#include "AdaptiveSetup.h"
#include "SIMoutput.h"
#include "ASMbase.h"
#include "ASMunstruct.h"
#include "Functions.h"
#include "IntegrandBase.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"
#ifdef HAS_LRSPLINE
#include <LRSpline/Basisfunction.h>
#endif
#include <fstream>
#include <sstream>
#include <cstdio>


AdaptiveSetup::AdaptiveSetup (SIMoutput& sim, bool sa) : model(sim), alone(sa)
{
  // Default grid adaptation parameters
  linIndep   = false;
  beta       = 10.0;
  betaFunc   = nullptr;
  errTol     = 1.0;
  rCond      = 1.0;
  condLimit  = 1.0e12;
  maxStep    = 10;
  maxDOFs    = 1000000;
  scheme     = FULLSPAN;
  knot_mult  = 1;
  threshold  = NONE;
  adaptor    = 0;
  adNorm     = 0;
  eRow       = 0;
  maxTjoints = -1;
  maxAspect  = -1.0;
  closeGaps  = false;
  symmEps    = 1.0e-6;
  storeMesh  = 1;
}


AdaptiveSetup::~AdaptiveSetup()
{
  delete betaFunc;
}


bool AdaptiveSetup::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"adaptive"))
    return alone ? model.parse(elem) : true;

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement()) {
    const char* value = utl::getValue(child,"maxstep");
    if (value)
      maxStep = atoi(value);
    else if ((value = utl::getValue(child,"maxdof")))
      maxDOFs = atoi(value);
    else if ((value = utl::getValue(child,"errtol")))
      errTol = atof(value);
    else if ((value = utl::getValue(child,"maxcondition")))
      condLimit = atof(value);
    else if ((value = utl::getValue(child,"knot_mult")))
      knot_mult = atoi(value);
    else if ((value = utl::getValue(child,"store_mesh"))) {
      mshPrefix = value;
      utl::getAttribute(child,"type",storeMesh);
    }
    else if (!strcasecmp(child->Value(),"store_eps_mesh")) {
      mshPrefix = "mesh";
      storeMesh = 30; // all postscript mesh files and lr-files
    }
    else if ((value = utl::getValue(child,"store_errors")))
      errPrefix = value;
    else if (!strcasecmp(child->Value(),"store_errors"))
      errPrefix = "error";
    else if (!strcasecmp(child->Value(),"test_linear_independence"))
      linIndep = true;
    else if ((value = utl::getValue(child,"scheme"))) {
      if (!strcasecmp(value,"fullspan"))
        scheme = FULLSPAN;
      else if (!strcasecmp(value,"minspan"))
        scheme = MINSPAN;
      else if (!strcasecmp(value,"isotropic_function"))
        scheme = ISOTROPIC_FUNCTION;
      else
        std::cerr <<"  ** AdaptiveSetup::parse: Unknown refinement scheme \""
                  << value <<"\" (ignored)"<< std::endl;
      utl::getAttribute(child,"maxTjoints",maxTjoints);
      utl::getAttribute(child,"maxAspectRatio",maxAspect);
      utl::getAttribute(child,"closeGaps",closeGaps);
      IFEM::cout <<"\tRefinement scheme: "<< scheme << std::endl;
    }
    else if ((value = utl::getValue(child,"use_norm"))) {
      adaptor = atoi(value);
      utl::getAttribute(child,"index",adNorm);
    }
    else if ((value = utl::getValue(child,"same_func_tolerance"))) {
#ifdef LRSPLINE_HAS_FUNC_TOLERANCE
      LR::Basisfunction::sameFuncTolerance = atof(value);
      IFEM::cout << "\tLRSpline: setting coinciding function tolerance to "
                 << LR::Basisfunction::sameFuncTolerance << std::endl;
#else
      IFEM::cout << "  ** LRSpline: setting coinciding function tolerance not supported" << std::endl;
#endif
    } else if ((value = utl::getValue(child,"use_sub_norm")))
      adNorm = atoi(value);
    else if ((value = utl::getValue(child,"beta"))) {
      std::string functype;
      utl::getAttribute(child,"functype",functype);
      if (functype == "expression") {
        IFEM::cout << "\tRefinement function: ";
        betaFunc = utl::parseTimeFunc(value,functype,1.0);
      } else
        beta = atof(value);
      std::string type;
      utl::getAttribute(child,"type",type,true);
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
      else if (type.compare("symmetrized") == 0 ||
               type.compare("symmetric") == 0) {
        threshold = SYMMETRIZED;
        utl::getAttribute(child,"eps",symmEps);
      }
      IFEM::cout <<"\tRefinement percentage: "<< (betaFunc ? (*betaFunc)(1) : beta)
                                              <<" type="<< threshold;
      if (threshold == SYMMETRIZED)
        IFEM::cout <<" (eps = "<< symmEps <<")";
      IFEM::cout << std::endl;
    }
  }

  return true;
}


bool AdaptiveSetup::parse (char* keyWord, std::istream& is)
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
      scheme = (RefScheme)itmp;
  }
  else
    return model.parse(keyWord,is);

  return true;
}


bool AdaptiveSetup::initPrm (size_t normGroup)
{
  if (normGroup > 0)
    adaptor = normGroup; // Override value from input file
  else if (adaptor == 0 && !model.haveAnaSol() && !model.opt.project.empty())
    adaptor = 1; // Set default to first projection if no analytical solution

  size_t iGroup = 0;
  IFEM::cout <<"\n\n >>> Starting adaptive simulation based on";
  for (const SIMoptions::ProjectionMap::value_type& prj : model.opt.project)
    if (++iGroup == adaptor)
    {
      IFEM::cout <<"\n     "<< prj.second
                 <<" error estimates (norm group "<< iGroup
                 <<") <<<"<< std::endl;
      break;
    }

  if (iGroup == adaptor && adaptor > 0)
  {
    if (!adNorm) adNorm = 2; // Use norm 2, unless specified
  }
  else if (model.haveAnaSol())
  {
    adaptor = 0; // Assuming the exact errors are stored in the first group
    if (!adNorm) adNorm = 4; // Use norm 4, unless specified
    IFEM::cout <<" exact errors ";
    if (adNorm != 4)
      IFEM::cout <<"(norm index "<< adNorm <<") ";
    IFEM::cout <<"<<<"<< std::endl;
  }
  else
  {
    IFEM::cout <<" - nothing, can not do that! <<<"<< std::endl;
    return false;
  }

  // Calculate row index in eNorm of the quantity to base the adaptation on
  eRow = adNorm;
  NormBase* norm = model.getNormIntegrand();
  for (iGroup = 1; iGroup <= adaptor; iGroup++)
    eRow += norm->getNoFields(iGroup);
  delete norm;

  return true;
}


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> DblIdx;


int AdaptiveSetup::calcRefinement (LR::RefineData& prm, int iStep,
                                   const Vectors& gNorm,
                                   const Vector& refIn) const
{
  prm.clear();

  ASMbase* thePatch = model.getPatch(1);
  if (!thePatch)
    return 0; // No patches in the model, nothing to do here
  else if (adaptor >= gNorm.size() || refIn.empty())
    return 0; // Refinement norm index out of range, or no element errors
  else if (adaptor == 0 && !model.haveAnaSol())
    return 0; // Cannot adapt on exact errors without an exact solution

  // Define the reference norm
  double refNorm = 0.01*model.getReferenceNorm(gNorm,adaptor);
  if (refNorm < -epsZ)
  {
    std::cerr <<" *** AdaptiveSetup::calcRefinement: Negative reference norm."
              <<" Check orientation of your model."<< std::endl;
    return -1; // Negative reference norm, modelling error?
  }
  else if (refNorm < epsZ)
    return 0; // Zero reference norm, probably no load on the model

  // Check if further refinement is required
  if (iStep > maxStep)
  {
    IFEM::cout << "\n   * Stopping the adaptive cycles as max steps " << maxStep << " was reached." << std::endl;
    return 0;
  }
  else if (model.getNoDOFs() > (size_t)maxDOFs)
  {
    IFEM::cout << "\n   * Stopping the adaptive cycles as max DOFs " << maxDOFs <<" was reached." << std::endl;
    return 0; // Refinement cycle or model size limit reached
  }
  else if (gNorm[adaptor](adNorm) < errTol*refNorm)
  {
    IFEM::cout << "\n   * Stopping the adaptive cycles as error tolerance " << errTol << " was reached." << std::endl;
    return 0; // Discretization error tolerance reached
  }
  else if (1.0/rCond > condLimit)
  {
    IFEM::cout <<"\n  ** Terminating the adaptive cycles due to instability."
               <<"\n     The last condition number "<< 1.0/rCond
               <<" is higher than the limit "<< condLimit << std::endl;
    return 0; // Condition number limit reached
  }

  prm.options.reserve(8);
  prm.options.push_back(betaFunc ? (*betaFunc)(iStep-1) : beta);
  prm.options.push_back(knot_mult);
  prm.options.push_back(scheme);
  prm.options.push_back(linIndep);
  prm.options.push_back(maxTjoints);
  prm.options.push_back(floor(maxAspect));
  prm.options.push_back(closeGaps);
  prm.options.push_back(threshold == TRUE_BETA ? 1 : 0);

  if (threshold == TRUE_BETA)
  {
    if (model.getNoPatches() > 1)
    {
      std::cerr <<" *** AdaptiveSetup::calcRefinement: True beta refinement"
                <<" is not available for multi-patch models."<< std::endl;
      return -2;
    }

    IFEM::cout <<"\nRefining by increasing solution space by "<< prm.options[0]
               <<" percent."<< std::endl;
    prm.errors.resize(thePatch->getNoRefineElms());
    dynamic_cast<ASMunstruct*>(thePatch)->remapErrors(prm.errors,refIn,true);
    return prm.errors.size();
  }

  size_t i, refineSize;
  std::vector<DblIdx> error;

  if (scheme == ISOTROPIC_FUNCTION) // use errors per function
  {
    if (thePatch->getNoRefineNodes() == thePatch->getNoNodes(1))
      error.resize(model.getNoNodes(1),DblIdx(0.0,0));
    else if (model.getNoPatches() == 1)
      error.resize(thePatch->getNoRefineNodes(),DblIdx(0.0,0));
    else
    {
      std::cerr <<" *** AdaptiveSetup::calcRefinement: Multi-patch refinement"
                <<" is not available for mixed models."<< std::endl;
      return -3;
    }

    for (i = 0; i < error.size(); i++)
      error[i].second = i;

    for (ASMbase* patch : model.getFEModel())
    {
      // Extract element norms for this patch
      RealArray locNorm(patch->getNoElms());
      for (i = 1; i <= locNorm.size(); i++)
        locNorm[i-1] = refIn(patch->getElmID(i));

      // Remap from geometry basis to refinement basis
      RealArray locErr(patch->getNoRefineNodes());
      dynamic_cast<ASMunstruct*>(patch)->remapErrors(locErr,locNorm);

      // Insert into global error array
      for (i = 0; i < locErr.size(); i++)
        if (model.getNoPatches() > 1)
          error[patch->getNodeID(i+1)-1].first += locErr[i];
        else
          error[i].first += locErr[i];
    }
  }

  else // use errors per element
  {
    if (model.getNoPatches() > 1) // not supported for multi-patch models
    {
      std::cerr <<" *** AdaptiveSetup::adaptMesh: Multi-patch refinement"
                <<" is available for isotropic_function only."<< std::endl;
      return -4;
    }

    for (i = 0; i < refIn.size(); i++)
      error.push_back(DblIdx(refIn(1+i),i));
  }

  // Sort the elements in the sequence of decreasing errors
  std::sort(error.begin(),error.end(),std::greater<DblIdx>());

  // Find the list of refinable elements or basis functions.
  // The variable prm.elements will contain one of the following:
  // - list of elements to be refined (if fullspan or minspan)
  // - list of basis functions to be refined (if structured mesh)
  double limit, sumErr = 0.0, curErr = 0.0;
  switch (threshold) {
  case MAXIMUM: // beta percent of max error (less than 100%)
    limit = error.front().first * prm.options[0]*0.01;
    break;
  case AVERAGE: // beta percent of avg error (typical 100%)
    for (const DblIdx& e : error) sumErr += e.first;
    limit = (sumErr/error.size()) * prm.options[0]*0.01;
    break;
  case MINIMUM: // beta percent of min error (more than 100%)
    limit = error.back().first * prm.options[0]*0.01;
    break;
  case DORFEL:
    limit = error.back().first;
    for (const DblIdx& e : error) sumErr += e.first;
    sumErr *= prm.options[0]*0.01;
    for (const DblIdx& e : error)
      if (curErr < sumErr)
        curErr += e.first;
      else
      {
        limit = e.first;
        break;
      }
    break;
  default:
    limit = 0.0;
  }

  if (threshold == NONE || threshold == SYMMETRIZED)
    refineSize = ceil(error.size()*prm.options[0]/100.0);
  else
    refineSize = std::upper_bound(error.begin(), error.end(), DblIdx(limit,0),
                                  std::greater_equal<DblIdx>()) - error.begin();

  if (threshold == SYMMETRIZED)
  {
    double fErr = error[refineSize-1].first;
    double fEps = fabs(fErr)*symmEps;
    while (refineSize < error.size() &&
           fabs(error[refineSize].first-fErr) < fEps)
      ++refineSize;
  }

  if (!errPrefix.empty())
  {
    char suffix[9];
    sprintf(suffix,"_%03d.txt",iStep-1);
    std::ofstream of(errPrefix+suffix);
    of.precision(16);
    of <<"# Index   Error\n";
    for (i = 0; i < error.size(); i++)
    {
      of << error[i].second <<" "<< error[i].first << std::endl;
      if (1+i == refineSize)
        of <<"\n# ------------- Below here not refined ---------#\n\n";
    }
  }

  const char* str = scheme < ISOTROPIC_FUNCTION ? "elements":"basis functions";
  IFEM::cout <<"\nRefining "<< refineSize <<" "<< str
             <<" with errors in range ["<< error[refineSize-1].first
             <<","<< error.front().first <<"] ";

  switch (threshold) {
  case NONE:
    IFEM::cout << prm.options[0] <<"% of all "<< str;
    break;
  case SYMMETRIZED:
    IFEM::cout << 100.0*refineSize/error.size() <<"% of all "<< str;
    break;
  case MAXIMUM:
    IFEM::cout << prm.options[0] <<"% of max error ("<< limit <<")";
    break;
  case AVERAGE:
    IFEM::cout << prm.options[0] <<"% of average error ("<< limit <<")";
    break;
  case MINIMUM:
    IFEM::cout << prm.options[0] <<"% of min error ("<< limit <<")";
    break;
  case DORFEL:
    IFEM::cout << prm.options[0] <<"% of total error ("<< limit <<")";
    break;
  default:
    break;
  }
  IFEM::cout << std::endl;

  prm.elements.reserve(refineSize);
  for (i = 0; i < refineSize; i++)
    prm.elements.push_back(error[i].second);

  return refineSize;
}


void AdaptiveSetup::printNorms (const Vectors& gNorm, const Vectors& dNorm,
                                const Matrix& eNorm, size_t w, bool printModelNorms) const
{
  if (printModelNorms)
    model.printNorms(gNorm,w);

  // Print out the norm used for mesh adaptation
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
    if (adaptor < dNorm.size())
    {
      const Vector& bNorm = dNorm[adaptor];
      IFEM::cout <<"\nError estimate, dual solution"
                 << utl::adjustRight(w-29,"(z)") << bNorm(adNorm)
                 <<"\nRelative error (%) : "
                 << 100.0*bNorm(adNorm)/hypot(dNorm.front()(1),bNorm(2))
                 <<"\nError estimate "
                 << utl::adjustRight(w-15,"E(u)*E(z), E(v)=a(v^r-v^h,v^r-v^h)")
                 << dNorm.front()(2);
    }
    IFEM::cout <<"\n"<< std::endl;

    delete norm;
  }

  if (eNorm.rows() < eRow || eNorm.cols() < 1)
    return;

  // Compute some additional error measures
  double avg_norm = eNorm(eRow,1);
  double min_err  = avg_norm;
  double max_err  = avg_norm;
  for (size_t i = 2; i <= eNorm.cols(); i++)
  {
    avg_norm += eNorm(eRow,i);
    if (min_err > eNorm(eRow,i))
      min_err = eNorm(eRow,i);
    else if (max_err < eNorm(eRow,i))
      max_err = eNorm(eRow,i);
  }
  avg_norm /= eNorm.cols();

  double RMS_norm = 0.0;
  for (size_t c = 1; c <= eNorm.cols(); c++)
    RMS_norm += pow(eNorm(eRow,c)-avg_norm,2.0);
  RMS_norm = sqrt(RMS_norm/eNorm.cols())/avg_norm;

  IFEM::cout <<  "Root mean square (RMS) of error      : "<< RMS_norm
             <<"\nMin element error                    : "<< min_err
             <<"\nMax element error                    : "<< max_err
             <<"\nAverage element error                : "<< avg_norm
             << std::endl;
}


bool AdaptiveSetup::writeMesh (int iStep) const
{
  if (mshPrefix.empty() || storeMesh < 1)
    return true;

  if (storeMesh%2)
  {
    // Write all patches to the same file
    char suffix[8];
    sprintf(suffix,"_%03d.lr",iStep);
    std::ofstream os(mshPrefix+suffix);
    model.dumpGeometry(os);
  }

  if (storeMesh > 2)
    for (ASMbase* pch : model.getFEModel())
    {
      char patchFile[256];
      if (model.getNoPatches() > 1)
        sprintf(patchFile,"%zu_%s_%03d",pch->idx+1,mshPrefix.c_str(),iStep);
      else
        sprintf(patchFile,"%s_%03d",mshPrefix.c_str(),iStep);
      dynamic_cast<ASMunstruct*>(pch)->storeMesh(patchFile,storeMesh/2);
    }

  return true;
}
