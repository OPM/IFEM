// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.C
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for linear isogeometric FEM simulators.
//!
//==============================================================================

#include "AdaptiveSIM.h"
#ifdef HAS_LRSPLINE
#include "ASMunstruct.h"
#else
#include "ASMbase.h"
#endif
#include "IntegrandBase.h"
#include "SIMbase.h"
#include "SIMenums.h"
#include "SystemMatrix.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <sstream>
#include <cstdio>


AdaptiveSIM::AdaptiveSIM (SIMbase* sim) : model(sim)
{
  // Default grid adaptation parameters
  storeMesh      = false;
  linIndepTest   = false;
  beta           = 10.0;
  errTol         = 1.0;
  maxStep        = 10;
  maxDOFs        = 1000000;
  scheme         = 0; // fullspan
  symmetry       = 1; // no symmetry
  knot_mult      = 1; // maximum regularity (continuity)
  trueBeta       = false; // beta measured in dimension increase
  threashold     = false; // beta generates a threshold (err/max{err})
  adaptor        = 0;
  maxTjoints     = -1;
  maxAspectRatio = -1.0;
  closeGaps      = false;
}


AdaptiveSIM::~AdaptiveSIM ()
{
  if (model) delete model;
}


bool AdaptiveSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"adaptive"))
    return model->parse(elem);

  const char* value = 0;
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
      utl::getAttribute(child, "maxAspectRatio", maxAspectRatio);
      if (child->Attribute("closeGaps"))
        closeGaps = true;
    }
    else if ((value = utl::getValue(child,"use_norm")))
      adaptor = atoi(value);
    else if ((value = utl::getValue(child,"beta"))) {
      beta = atof(value);
      std::string type;
      utl::getAttribute(child, "type", type, true);
      if (type.compare("threshold") == 0 || type.compare("threashold") == 0)
        threashold = true;
      else if (type.compare("truebeta") == 0)
        trueBeta = true;
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
    // (0=refine all, 1=minimum span, 2=isotropic)
    if (!cline.fail() && !cline.bad())
      scheme = itmp;

    cline >> itmp; // read symmetry
    // (0=none, <n> always requests a multiplum of <n> elements for refinement)
    if (!cline.fail() && !cline.bad())
      symmetry = itmp;
  }
  else
    return model->parse(keyWord,is);

  return true;
}


bool AdaptiveSIM::initAdaptor (size_t indxProj, size_t nNormProj)
{
  if (indxProj > 0)
    adaptor = indxProj; // override value from XML input

  SIMoptions::ProjectionMap::const_iterator pit = model->opt.project.begin();
  for (size_t j = 1; pit != model->opt.project.end(); pit++, j++)
    if (j == adaptor) break;

  std::cout <<"\n\n >>> Starting adaptive simulation based on";
  if (pit != model->opt.project.end())
    std::cout <<"\n     "<< pit->second
              <<" error estimates (norm group "<< adaptor <<") <<<"<< std::endl;
  else if (model->haveAnaSol())
  {
    std::cout <<" exact errors <<<"<< std::endl;
    adaptor = 0; // Assuming the exact errors are stored in the first group
  }
  else
  {
    std::cout <<" - nothing, can not do that!"<< std::endl;
    return false;
  }

  if (model->opt.format >= 0)
    prefix.reserve(model->opt.project.size()+1);

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, int iStep)
{
  std::cout <<"\nAdaptive step "<< iStep << std::endl;
  if (iStep > 1)
  {
    SIMoptions opt = model->opt;
    // Re-generate the FE model after the refinement
    model->clearProperties();
    // Caution: If we are using XML-input and have specified old command-line
    // options in order to override simulation options read from the input file,
    // those options will not be overridden here, so please don't do that..
    if (!model->read(inputfile) || !model->preprocess())
      return false;
    model->opt = opt;
  }
  else if (storeMesh)
    // Output the initial grid to eps-file
    model->refine(std::vector<int>(),std::vector<int>(),"mesh_001.eps");

  // Assemble the linear FE equation system
  model->setMode(SIM::STATIC,true);
  model->initSystem(iStep == 1 ? SystemMatrix::DENSE : model->opt.solver, 1, 1);
  model->setQuadratureRule(model->opt.nGauss[0],true);
  if (!model->assembleSystem())
    return false;

  // Solve the linear system of equations
  if (!model->solveSystem(linsol,1))
    return false;

  eNorm.clear();
  gNorm.clear();

  // Project the secondary solution onto the splines basis
  model->setMode(SIM::RECOVERY);
  SIMoptions::ProjectionMap::const_iterator pit;
  int p=0;
  for (pit = model->opt.project.begin(); pit != model->opt.project.end(); pit++)
  {
    Matrix ssol;
    if (!model->project(ssol,linsol,pit->first))
      return false;

    projs[p++] = ssol;
    if (iStep == 1 && model->opt.format >= 0)
      prefix.push_back(pit->second.c_str());
  }
  if (iStep == 1 && model->opt.format >= 0)
    prefix.push_back(NULL);

  if (msgLevel > 1 && !projs.empty())
    std::cout << std::endl;

  // Evaluate solution norms
  model->setQuadratureRule(model->opt.nGauss[1]);
  return (model->solutionNorms(Vectors(1,linsol),projs,eNorm,gNorm) &&
	  model->dumpResults(linsol,0.0,std::cout,true,6));
}


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> IndexDouble;


bool AdaptiveSIM::adaptMesh (int iStep)
{
  this->printNorms(std::cout);

  if (adaptor >= gNorm.size() || adaptor >= eNorm.rows())
    return false;

  // Cannot adapt on exact errors without an exact solution
  if (adaptor == 0 && !model->haveAnaSol())
    return false;

  // Define the reference norm
  double refNorm;
  const Vector& fNorm = gNorm.front();
  const Vector& aNorm = gNorm[adaptor];
  if (adaptor == 0)
    refNorm = fNorm(3); // Using the analytical solution, |u|_ref = |u|
  else // |u|_ref = sqrt( |u^h|^2 + |e^*|^2 )
    refNorm = sqrt(fNorm(1)*fNorm(1) + aNorm(2)*aNorm(2));

  // Use norm 4 when adapting on exact errors, otherwise use norm 2
  int adNorm = adaptor == 0 ? 4 : 2;

  std::cout << " 100.0*aNorm(adNorm) < errTol*refNormA == \""
            <<   100.0*aNorm(adNorm) << " < " <<  errTol*refNorm << "\"\n";
  std::cout << "model->getNoDOFs() > (size_t)maxDOFs) == \""
            <<  model->getNoDOFs() << " > " <<  (size_t)maxDOFs << "\"\n";

  // Check if further refinement is required
  if (iStep > maxStep || model->getNoDOFs() > (size_t)maxDOFs) return false;
  if (eNorm.cols() < 1 || 100.0*aNorm(adNorm) < errTol*refNorm) return false;

  std::cout << "Checkpoint #1\n";

  // Calculate row index in eNorm of the error norm to adapt based on
  size_t i, eRow = adNorm;
  NormBase* norm = model->getNormIntegrand();
  for (i = 0; i < adaptor; i++)
    eRow += norm->getNoFields(i+1);
  delete norm;

  std::vector<int> toBeRefined, options;
  std::vector<IndexDouble> errors;

  options.reserve(8);
  options.push_back(beta);
  options.push_back(knot_mult);
  options.push_back(scheme);
  options.push_back(linIndepTest);
  options.push_back(maxTjoints);
  options.push_back(floor(maxAspectRatio));
  options.push_back(closeGaps);
  options.push_back(trueBeta);

  std::cout << "Checkpoint #1.5\n";
#ifdef HAS_LRSPLINE
  if (trueBeta) {
    std::cout <<"\nRefining by increasing solution space by "<< beta <<" percent\n";
    ASMunstruct* patch = static_cast<ASMunstruct*>(model->getFEModel().front());
    if (!storeMesh)
      return patch->refine(eNorm.getRow(eRow),options);
    std::cout << "Checkpoint #2\n";

    char fname[13];
    sprintf(fname,"mesh_%03d.eps",iStep);
    return patch->refine(eNorm.getRow(eRow),options,fname);
    std::cout << "Checkpoint #3\n";
  }
#endif

  if (scheme == 2)
  {
    // Sum up the total error over all supported elements for each function
    ASMbase* patch = model->getFEModel().front();
    IntMat::const_iterator eit;
    IntVec::const_iterator nit;
    for (i = 0; i < patch->getNoNodes(); i++) // and by "Nodes", we mean basis functions
      errors.push_back(IndexDouble(0.0,i));
    for (i = 1, eit = patch->begin_elm(); eit < patch->end_elm(); eit++, i++)
      for (nit = eit->begin(); nit < eit->end(); nit++)
        errors[*nit].first += eNorm(eRow,i);
  }
  else
    for (i = 0; i < eNorm.cols(); i++)
      errors.push_back(IndexDouble(eNorm(eRow,1+i),i));

  // sort errors
  std::sort(errors.begin(),errors.end(),std::greater<IndexDouble>());

  // find the list of refinable elements/basisfunctions
  // the variable 'toBeRefined' contains one of the following:
  //   - list of elements to be refined (if fullspan or minspan)
  //   - list of basisfunctions to be refined (if structured mesh)
  //   - nothing if trueBeta (function returned above)
  size_t refineSize;
  if (threashold) {
    std::vector<IndexDouble>::const_iterator it;
    it = std::upper_bound(errors.begin(), errors.end(),
                          IndexDouble(errors.front().first * beta/100.0,0),
                          std::greater<IndexDouble>());
    refineSize = it - errors.begin();
  }
  else
    refineSize = ceil(errors.size()*beta/100.0);

  std::cout <<"\nRefining "<< refineSize
            << (scheme < 2 ? " elements" : " basis functions")
            << (threashold ? " (threshold of " : " (") << beta <<"\%)"
            <<" with errors in range ["<< errors[refineSize-1].first
            <<","<< errors.front().first <<"]"<< std::endl;

/*
  if (symmetry > 0) // Make refineSize a multiplum of 'symmetry' in case of symmetric problems
    refineSize += (symmetry-refineSize%symmetry);
*/

  if (refineSize < 1 || refineSize > errors.size()) return false;

  std::cout << "Checkpoint #4\n";

  toBeRefined.reserve(refineSize);
  for (i = 0; i < refineSize; i++)
    toBeRefined.push_back(errors[i].second);

  // Now refine the mesh
  if (!storeMesh)
    return model->refine(toBeRefined,options);
  std::cout << "Checkpoint #5\n";

  char fname[13];
  sprintf(fname,"mesh_%03d.eps",iStep);
  return model->refine(toBeRefined,options,fname);
}


std::ostream& AdaptiveSIM::printNorms (std::ostream& os) const
{
  model->printNorms(gNorm,os);

  // TODO: This needs further work to enable print for all recovered solutions.
  // As for now we only print the norm used for the mesh adaption.
  size_t i, eRow = 4;
  if (adaptor > 0 && adaptor < gNorm.size())
  {
    NormBase* norm = model->getNormIntegrand();

    const Vector& fNorm = gNorm.front();
    const Vector& aNorm = gNorm[adaptor];

    // Define the reference norm, |u|_ref = sqrt( |u^h|^2 + |e^*|^2 )
    double refNorm = sqrt(fNorm(1)*fNorm(1) + aNorm(2)*aNorm(2));

    os <<"Error estimate "<< norm->getName(adaptor,2) <<": "<< aNorm(2)
       <<"\nRelative error (%) : "<< 100.0*aNorm(2)/refNorm;
    if (model->haveAnaSol() && aNorm.size() > 2)
      os <<"Projective error "<< norm->getName(adaptor,3) <<": "
                              << aNorm(3) << std::endl;
    if (model->haveAnaSol() && fNorm.size() > 3)
      os <<"Effectivity index  : "<< aNorm(2)/fNorm(4);

    for (i = 0, eRow = 2; i < adaptor; i++)
      eRow += norm->getNoFields(i+1);

    delete norm;
  }

  if (eNorm.rows() < eRow || eNorm.cols() < 1)
    return os << std::endl;

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

  os <<"\nRoot mean square (RMS) of error      : "<< RMS_norm
     <<"\nMin element error                    : "<< min_err
     <<"\nMax element error                    : "<< max_err
     <<"\nAverage element error                : "<< avg_norm;
  return os << std::endl;
}


bool AdaptiveSIM::writeGlv (const char* infile, int iStep, int& nBlock,
			    size_t nNormProj)
{
  if (model->opt.format < 0) return true;

  // Write VTF-file with model geometry
  if (!model->writeGlvG(nBlock, iStep == 1 ? infile : 0))
    return false;

  // Write boundary tractions, if any
  if (!model->writeGlvT(iStep,nBlock))
    return false;

  // Write Dirichlet boundary conditions
  if (!model->writeGlvBC(nBlock,iStep))
    return false;

  // Write solution fields
  if (!model->writeGlvS(linsol,iStep,nBlock))
    return false;

  // Write projected solution fields
  SIMoptions::ProjectionMap::const_iterator pit = model->opt.project.begin();
  for (size_t i = 0; i < projs.size(); i++, pit++)
    if (!model->writeGlvP(projs[i],iStep,nBlock,100+10*i,pit->second.c_str()))
      return false;

  // Write element norms
  if (model->haveAnaSol()) nNormProj += 2;
  if (!model->writeGlvN(eNorm,iStep,nBlock,&prefix.front()))
    return false;

  // Write state information
  return model->writeGlvStep(iStep,iStep,1);
}


void AdaptiveSIM::setupProjections()
{
  projs.resize(model->opt.project.size());
}


void AdaptiveSIM::setOptions(SIMoptions& opt2)
{
  model->opt = opt = opt2;
}
